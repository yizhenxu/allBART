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

  void mypbart(double  *pX, double *testpX, double *mu,
               int *pn, double *y,
               int *pn_cov,
               int *testn,
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
               double *vec_test,
               double *vec_train,
               int *binaryX){
    //    *pX is n_samp x n_cov  matrix
    //    *y is binary

    dinfo di; dinfo dip;
    di.n_samp = *pn; di.n_cov = *pn_cov; di.n_dim = 0;
    if(*testn){
      dip.n_samp = *testn; dip.n_cov = *pn_cov; dip.n_dim = 0; dip.y=0;
    }

    size_t minobsnode = *pminobsnode;

    std::vector<std::vector<double> >  XMat; /* The train.data of dimensions nsub x ncov*/
    int itemp = 0;
    for(size_t i=0; i < di.n_samp; i++){
      for(size_t k=0; k< di.n_cov; k++){
        XMat[i][k] = pX[itemp++];
      }
    }



    std::vector<std::vector<double> >   testXMat; /* The test.data of dimensions nsub_test x ncov*/
    if(*testn){
      itemp = 0;
      for(size_t i=0; i < dip.n_samp; i++){
        for(size_t k=0; k< dip.n_cov; k++){
          testXMat[i][k] = testpX[itemp++];
        }
      }
    }


    xinfo xi; /* The cutpoint matrix (ncov x nc) */
    int nc=*pnc; // number of equally spaced cutpoints between min and max.
    getcutpoints(nc, (int)di.n_cov, (int)di.n_samp, XMat, xi);

    std::vector<double> allfit; /* The sum of fit of all trees for train.data (nsub x 1)*/
    allfit.resize(di.n_samp);

    std::vector<double> w; /* The latent variable for train.data (nsub x 1)*/
    w.resize(di.n_samp);

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
    size_t nd = *pndraws; //number of mcmc iterations
    size_t m=*pntrees;
    double kfac=*pkfac;

    pinfo pi;
    pi.pswap = *ppswap; //prob of swap move, default 0.1
    pi.pbd=*ppbd; //prob of birth/death move, default 0.5
    pi.pb=*ppb; //prob of birth given  birth/death, default 0.25

    pi.alpha=*palpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
    pi.beta=*pbeta; //
    pi.tau=3.0/(kfac*sqrt((double)m)); // categorical outcome
    pi.sigma= 1; //error standard deviation of the latent W

    //initialize tree

    std::vector<tree>  t; /* ntree vector of trees */
    t.resize(m);
    for(size_t i=0; i<m; i++){
      t[i].setm(0.00);
    }

    for(size_t i=0;i<di.n_samp;i++) {
      allfit[i] = 0.0;
      w[i] = 2*(y[i]) - 1.0; //y is 0 and 1, w is -1 and 1
    }

    double u,Z; //for updating latent w

    //MCMC


    //cout << "\nMCMC:\n";
    time_t tp;
    int time1 = time(&tp);

    /* Initialize counters for outputs vec_test and vec_train */
    int countvectest = 0;
    int countvectrain = 0;


    for(size_t loop=0;loop<(nd+burn);loop++) { /* Start posterior draws */


    if(loop%100==0) Rprintf("iteration: %d of %d \n",loop, nd+burn);

    /* Step 1 sample trees*/
    /* See tree sampling theory.doc for explanation */
    for(size_t ntree = 0 ; ntree <m; ntree++){
      fit(t[ntree], XMat, di, xi, ftemp);
      for(size_t i=0;i<di.n_samp;i++) {
        allfit[i] -= ftemp[i];
        rtemp[i] = w[i] - allfit[i];
      }

      di.y = &rtemp[0];
      bd(XMat, t[ntree], xi, di, pi, minobsnode, binaryX);
      fit(t[ntree], XMat, di, xi, ftemp);
      for(size_t i=0;i<di.n_samp;i++) {
        allfit[i] += ftemp[i];
      }

    }//ntree

    //done sampling (T,M)

    /* Step 2 update latent variable w*/
    for(size_t i=0;i<di.n_samp;i++) {
      u = unif_rand();
      if(y[i] > 0) {
        Z = qnorm((1.0-u)*pnorm(-allfit[i],0.0,1.0,1,0) + u,0.0,1.0,1,0);
      } else {
        Z = -qnorm((1.0-u)*pnorm(allfit[i],0.0,1.0,1,0) + u,0.0,1.0,1,0);
      }
      w[i] = allfit[i] + Z;
    }

    if(loop>=burn){
      for(size_t k = 0; k <di.n_samp; k++){
        vec_train[countvectrain] = allfit[k];
        countvectrain++;
      }//end prediction for train

      if(*testn) {

        for(size_t k=0; k<dip.n_samp; k++){
          ppredmeanvec[k] = 0.0;
        }

        for(size_t j=0;j<m;j++) {
          fit(t[j], testXMat, dip, xi, fpredtemp);
          for(size_t k=0;k<dip.n_samp;k++) ppredmeanvec[k] += fpredtemp[k];
        }

        for(size_t k = 0; k <dip.n_samp; k++){
          vec_test[countvectest] = ppredmeanvec[k];
          countvectest++;
        }//end prediction for test

      }//end if test

    }//end prediction for current loop

    } //end of loop

    int time2 = time(&tp);
    Rprintf("time for mcmc loop %d secs", time2-time1);

    for(size_t i=0; i<m; i++){
      t[i].tonull();
    }//delete trees
  }
};
