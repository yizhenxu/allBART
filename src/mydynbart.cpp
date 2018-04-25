#include <iostream>
#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include "funs.h"
#include "tree.h"
#include "info.h"


extern "C" {

  void mydynbart(double  *pX, double *testpX, double *mu, double *psigest,
              int *pn, double *y,
              int *pn_cov, int *pnu, double *plambda,
              int *testn, int *testnsub,
              int *pndraws,
              int *pburn,
              int *pntrees,
              double *pkfac,
              double *ppswap,
              double *ppbd,
              double *ppb,
              double *palpha,
              double *pbeta,
              int *pnc,
              int *pminobsnode,
              double *psigmasample,
              double *vec_test,
              double *vec_train,
              int *binaryX){
    //    *pX is n_samp x n_cov  matrix
    //    *y is continuous

    dinfo di; dinfo dip;
    di.n_samp = *pn; di.n_cov = *pn_cov; di.n_dim = 0;
    if(*testn){
      dip.n_samp = *testnsub; dip.n_cov = *pn_cov; dip.n_dim = 0; dip.y=0;
    }

    int nn = *pn; /* nsub */
    int nu = *pnu; /* df in sigma prior */
    double lambda = *plambda; /* param in sigma prior */
    size_t minobsnode = *pminobsnode;

    std::vector<std::vector<double> > XMat; /* The train.data of dimensions nsub x ncov*/
    XMat.resize(di.n_samp);
    for(size_t j=0; j < di.n_samp; j++){
      XMat[j].resize(di.n_cov);
    }
    int itemp = 0;
    for(size_t i=0; i < di.n_samp; i++){
      for(size_t k=0; k< di.n_cov; k++){
        XMat[i][k] = pX[itemp++];
      }
    }

    size_t nd = *pndraws; //number of mcmc iterations

    std::vector<std::vector<std::vector<double> > > testXMat; /* The test.data of dimensions nsub_test x ncov*/
    testXMat.resize(nd);
    for(size_t j=0; j < nd; j++){
      testXMat[j].resize(dip.n_samp);
    }

    for(size_t j=0; j < nd; j++){
      for(size_t i=0; i < dip.n_samp; i++){
        testXMat[j][i].resize(dip.n_cov);
      }
    }
    itemp = 0;
    for(size_t j=0; j < nd; j++){
      for(size_t i=0; i < dip.n_samp; i++){
        for(size_t k=0; k< dip.n_cov; k++){
          testXMat[j][i][k] = testpX[itemp++];
        }
      }
    }


    xinfo xi; /* The cutpoint matrix (ncov x nc) */
    int nc=*pnc; // number of equally spaced cutpoints between min and max.
    getcutpoints(nc, (int)di.n_cov, (int)di.n_samp, XMat, xi);

    std::vector<double> allfit; /* The sum of fit of all trees for train.data (nsub x 1)*/
    allfit.resize(di.n_samp);

    std::vector<double> ppredmeanvec; /* The sum of fit of all trees for test.data (nsub_test x 1)*/
    if(*testn){
      ppredmeanvec.resize(dip.n_samp);
    }

    /* fit of current tree */
    std::vector<double> ftemp; /* for allfit  (nsub x 1)*/
    ftemp.resize(di.n_samp);

    std::vector<double> fpredtemp; /* for ppredmeanvec (nsub_test x 1) */
    if(*testn){
      //temporary fit vector to compute prediction
      fpredtemp.resize(dip.n_samp);
    }

    //partial residuals from backfitting
    std::vector<double> rtemp;
    rtemp.resize(di.n_samp);


    // priors and parameters
    size_t burn = *pburn; //number of mcmc iterations called burn-in
    //size_t nd = *pndraws; //number of mcmc iterations
    size_t m=*pntrees;
    double kfac=*pkfac;

    pinfo pi;
    pi.pswap = *ppswap; //prob of swap move, default 0.1
    pi.pbd=*ppbd; //prob of birth/death move, default 0.5
    pi.pb=*ppb; //prob of birth given  birth/death, default 0.25

    pi.alpha=*palpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
    pi.beta=*pbeta; //
    pi.tau=(0.5)/(kfac*sqrt((double)m)); //this is different from categorical outcome
    pi.sigma= *psigest;

    //initialize tree

    std::vector<tree>  t; /* ntree vector of trees */
    t.resize(m);
    for(size_t i=0; i<m; i++){
      t[i].setm(0.00);
    }

    for(size_t i=0;i<di.n_samp;i++) {
      allfit[i] = 0.0;
    }

    double ss = 0.0;

    //MCMC

    time_t tp;
    int time1 = time(&tp);

    /* Initialize counters for outputs psigmasample, vec_test and vec_train */
    int sigdrawcounter = 0;
    int countvectest = 0;
    int countvectrain = 0;


    for(size_t loop=0;loop<(nd+burn);loop++) { /* Start posterior draws */

    GetRNGstate();

      if(loop%100==0) Rprintf("\n iteration: %d of %d \n",loop, nd+burn);

      /* Step 1 */
      /* See tree sampling theory.doc for explanation */
      for(size_t ntree = 0 ; ntree <m; ntree++){
        fit(t[ntree], XMat, di, xi, ftemp);
        for(size_t i=0;i<di.n_samp;i++) {
          allfit[i] -= ftemp[i];
          rtemp[i] = y[i] - allfit[i];
        }

        di.y = &rtemp[0];
        bd(XMat, t[ntree], xi, di, pi, minobsnode, binaryX);

        fit(t[ntree], XMat, di, xi, ftemp);
        for(size_t i=0;i<di.n_samp;i++) {
          allfit[i] += ftemp[i];
        }

      }//ntree


      //done sampling (T,M)

      /* Step 2 */
      for(size_t i=0;i<di.n_samp;i++) {
        rtemp[i] = y[i] - allfit[i];
      }

      ss = 0.0;
      for(size_t i=0;i<di.n_samp;i++) {
        ss += rtemp[i]* rtemp[i];
      }

      int nupost = nu+nn;
      double nlpost = nu*lambda + ss;
      pi.sigma = sqrt(nlpost/rchisq((double)nupost));

      //std::cout << "loop = "<<loop<<"; SSE = "<<ss<<"; sigma = " << pi.sigma<<";\n";

      if(loop>=burn){
        psigmasample[sigdrawcounter++] = pi.sigma;

        for(size_t k = 0; k <di.n_samp; k++){
          vec_train[countvectrain] = allfit[k];
          countvectrain++;
        }//end prediction for train

        if(*testn) {

          for(size_t k=0; k<dip.n_samp; k++){
            ppredmeanvec[k] = 0.0;
          }

          for(size_t j=0;j<m;j++) {
            fit(t[j], testXMat[loop-burn], dip, xi, fpredtemp);
            for(size_t k=0;k<dip.n_samp;k++) ppredmeanvec[k] += fpredtemp[k];
          }

          for(size_t k = 0; k <dip.n_samp; k++){
            vec_test[countvectest] = ppredmeanvec[k];
            countvectest++;
          }//end prediction for test

        }//end if test

      }//end prediction for current loop

      PutRNGstate();

    } //end of loop

    int time2 = time(&tp);
    Rprintf("time for mcmc loop %d secs", time2-time1);

    for(size_t i=0; i<m; i++){
      t[i].tonull();
    }//delete trees

  }
};
