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

  void mydynpbart(double  *pX, double *testpX, double *mu,
               int *pn, double *y,
               int *pn_cov,
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
               double *vec_test,
               double *vec_train,
               int *binaryX,
               int *diagnostics,
               double *percA,
               double *numNodes,
               double *numLeaves,
               double *treeDepth,
               double *incProp){
    //    *pX is n_samp x n_cov  matrix
    //    *y is binary

    dinfo di; dinfo dip;
    di.n_samp = *pn; di.n_cov = *pn_cov; di.n_dim = 0;
    if(*testn){
      dip.n_samp = *testnsub; dip.n_cov = *pn_cov; dip.n_dim = 0; dip.y=0;
    }

    size_t minobsnode = *pminobsnode;

    std::vector<std::vector<double> >  XMat; /* The train.data of dimensions nsub x ncov*/
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
    //size_t nd = *pndraws; //number of mcmc iterations
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

    int dgn = *diagnostics;
    //percA[i] = average over all percAtmp for draw i
    std::vector<double> percAtmp;
    percAtmp.resize(m);
    //numNodes[i] = average over all numNodestmp for draw i
    std::vector<double> numNodestmp;
    numNodestmp.resize(m);
    //numLeaves[i] = average over all numLeavestmp for draw i
    std::vector<double> numLeavestmp;
    numLeavestmp.resize(m);
    //treeDepth[i] = average over all treeDepthtmp for draw i
    std::vector<double> treeDepthtmp;
    treeDepthtmp.resize(m);
    //nLtD[i] = number of leaves at tree depth of the ith tree for the current round of trai
    std::vector<double> nLtDtmp;
    nLtDtmp.resize(m);

    for(size_t i=0; i<m; i++){
      percAtmp[i] = 0.0;
      numNodestmp[i] = 1.0;
      numLeavestmp[i] = 1.0;
      treeDepthtmp[i] = 0.0;
      nLtDtmp[i] = 1.0;
    }

    //incProptmp[j] = count occurence of xj in splitting rules for the current round of draw
    //incProp[j] = average of incProptmp[j] over all draws
    std::vector<double> incProptmp;
    incProptmp.resize(di.n_cov);

    for(size_t i=0; i<di.n_cov; i++){
      incProptmp[i] = 0.0;
    }

    //MCMC


    //cout << "\nMCMC:\n";
    time_t tp;
    int time1 = time(&tp);

    /* Initialize counters for outputs vec_test and vec_train */
    int countvectest = 0;
    int countvectrain = 0;


    for(size_t loop=0;loop<(nd+burn);loop++) { /* Start posterior draws */
    GetRNGstate();
      if(loop%100==0) Rprintf("\n iteration: %d of %d \n",loop, nd+burn);

      /* Step 1 sample trees*/
      /* See tree sampling theory.doc for explanation */
      for(size_t ntree = 0 ; ntree <m; ntree++){
        fit(t[ntree], XMat, di, xi, ftemp);
        for(size_t i=0;i<di.n_samp;i++) {
          allfit[i] -= ftemp[i];
          rtemp[i] = w[i] - allfit[i];
        }

        di.y = &rtemp[0];

        if(dgn){
          bd1(XMat, t[ntree], xi, di, pi, minobsnode, binaryX, &nLtDtmp[ntree], &percAtmp[ntree], &numNodestmp[ntree], &numLeavestmp[ntree], &treeDepthtmp[ntree], incProptmp);
        } else {
          bd(XMat, t[ntree], xi, di, pi, minobsnode, binaryX);
        }

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

        if(dgn){
          for(size_t i=0;i<m;i++) {
            percA[loop-burn] += percAtmp[i];
            numNodes[loop-burn] += numNodestmp[i];
            numLeaves[loop-burn] += numLeavestmp[i];
            treeDepth[loop-burn] += treeDepthtmp[i];
          }
          percA[loop-burn] /= m;
          numNodes[loop-burn] /= m;
          numLeaves[loop-burn] /= m;
          treeDepth[loop-burn] /= m;

          double numsplits = 0;
          for(size_t i=0; i<di.n_cov; i++){
            //std::cout<< "loop"<<loop<<" "<<incProptmp[i] <<" \n";
            numsplits += incProptmp[i];
          }
          for(size_t i=0; i<di.n_cov; i++){
            incProp[i] += incProptmp[i]/numsplits;//proportion of all splitting rules that uses xi
          }
        }

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

    if(dgn){
      for(size_t i=0; i<di.n_cov; i++){
        incProp[i] /= nd;
      }
    }

    int time2 = time(&tp);
    Rprintf("time for mcmc loop %d secs", time2-time1);


    for(size_t i=0; i<m; i++){
      t[i].tonull();
    }//delete trees
  }
};
