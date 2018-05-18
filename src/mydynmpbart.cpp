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

  void mydynmpbart(double *w, double  *pX, double *testpX, double *mu, double *sigmai, double *V,
                   int *pn, int *pn_dim, int *y, int *pn_cov, int *pnu,
                   int *testn,
                   int *testnsub,
                   int *pndraws,
                   int *pburn,
                   int *pntrees,
                   int *pSigDr,
                   double *pkfac,
                   double *ppswap,
                   double *ppbd,
                   double *ppb,
                   double *palpha,
                   double *pbeta,
                   int *pnc,
                   int *psavesigma,
                   int *pminobsnode,
                   double *psigmasample,
                   int *vec_class_pred_test,
                   int *vec_class_pred_train,
                   int *binaryX,
                   int *diagnostics,
                   double *percA,
                   double *numNodes,
                   double *numLeaves,
                   double *treeDepth,
                   double *incProp){
    // w is the starting value of latents

    //    *w is n_samp x n_dim vector
    //    *pX is (n_samp x n_dim) x n_cov  matrix. eg. row 0 corresponds to first choic, row 1 corresponds to 2nd choice
    //     *y is multinomial 1,..., (n_dim + 1)
    //  *sigmai is (n_dim) x (n_dim)


    dinfo di; dinfo dip;
    di.n_samp = *pn; di.n_cov = *pn_cov; di.n_dim = *pn_dim;
    if(*testn){
      dip.n_samp = *testnsub; dip.n_cov = *pn_cov; dip.n_dim = *pn_dim; dip.y=0;
    }


    //Initialize latents
    int nn = *pn; /* nsub */
    int nndim = *pn_dim; /* nlatent */
    int nu = *pnu;
    int savesigma= *psavesigma;
    size_t minobsnode = *pminobsnode;


    double maxy; /* The reference level */
    maxy = R_NegInf;
    for(size_t i=0; i<di.n_samp;i++){
      if(maxy<y[i]) maxy = y[i];
    }

    ////
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
    ////

    std::vector<std::vector<double> > allfit; /* The sum of fit of all trees for train.data (nlatent x nsub)*/
    std::vector<std::vector<double> > wtilde; /* (nlatent x nsub) */

    std::vector<std::vector<double> > ppredmeanvec; /* The sum of fit of all trees for test.data (nlatent x nsub_test)*/
    if(*testn){

      ppredmeanvec.resize(dip.n_dim);
      for(size_t k=0; k<dip.n_dim; k++){
        ppredmeanvec[k].resize(dip.n_samp);
      }
    }
    /* fit of current tree */
    std::vector<std::vector<double> > ftemp; /* for allfit  (nlatent x nsub)*/
    ftemp.resize(di.n_dim);
    std::vector<double> fpredtemp; /* for ppredmeanvec (nsub_test vector) */
    if(*testn){
      //temporary fit vector to compute prediction
      fpredtemp.resize(dip.n_samp);
    }
    std::vector<double> mvnsample,mvnmean; /* temp storage for mvn samples (nlatent vector) */
    mvnsample.resize(di.n_dim);
    mvnmean.resize(di.n_dim);
    allfit.resize(di.n_dim);
    wtilde.resize(di.n_dim);



    for(size_t k=0; k<di.n_dim; k++){
      ftemp[k].resize(di.n_samp);
      allfit[k].resize(di.n_samp);
      wtilde[k].resize(di.n_samp);
    }


    std::vector<std::vector<double> > r; /* pseudoresponse (nlatent x nsub) */
    std::vector<std::vector<double> > rtemp;
    rtemp.resize(di.n_dim);
    r.resize(di.n_dim);
    for(size_t k=0; k<di.n_dim; k++){
      r[k].resize(di.n_samp);
      rtemp[k].resize(di.n_samp);
    }


    std::vector<double> condsig;
    condsig.resize(di.n_dim);

    // priors and parameters
    size_t burn = *pburn; //number of mcmc iterations called burn-in
    //size_t nd = *pndraws; //number of mcmc iterations
    size_t m=*pntrees;
    size_t nSigDr = *pSigDr;
    double kfac=*pkfac;

    pinfo pi;
    pi.pswap = *ppswap; //prob of swap move, default 0.1
    pi.pbd=*ppbd; //prob of birth/death move, default 0.5
    pi.pb=*ppb; //prob of birth given  birth/death, default 0.25

    pi.alpha=*palpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
    pi.beta=*pbeta; //
    pi.tau=(3.0)/(kfac*sqrt((double)m));
    pi.sigma=1.0;

    //storage for ouput
    //in sample fit

    double max_temp = 0.0;// used to find class membership
    int pclass = 0;// used to find class membership

    //initialize tree

    std::vector<std::vector<tree> > t; /* ntree x nlatent matrix of trees */
    t.resize(m);
    for(size_t i=0; i<m; i++){
      t[i].resize(di.n_dim);
    }
    for(size_t i=0; i<m; i++){
      for(size_t k=0;k< di.n_dim;k++) t[i][k].setm(0.00);
    }

    for(size_t k=0; k<di.n_dim; k++){
      for(size_t i=0;i<di.n_samp;i++) {
        allfit[k][i] = 0.0;
        wtilde[k][i] = 0.0;
      }
    }

    std::vector<std::vector<double> > mtemp1;
    std::vector<std::vector<double> > WishSample,WishSampleTilde,WishSampleInv, SigmaTmp, SigmaTmpInv;
    WishSampleInv.resize(di.n_dim);
    WishSample.resize(di.n_dim);
    SigmaTmp.resize(di.n_dim);
    SigmaTmpInv.resize(di.n_dim);

    WishSampleTilde.resize(di.n_dim);

    WishSampleInv.resize(di.n_dim);
    mtemp1.resize(di.n_dim);
    for(size_t j=0;j<di.n_dim;j++){
      WishSample[j].resize(di.n_dim);
      WishSampleTilde[j].resize(di.n_dim);
      mtemp1[j].resize(di.n_dim);
      WishSampleInv[j].resize(di.n_dim);
      SigmaTmp[j].resize(di.n_dim);
      SigmaTmpInv[j].resize(di.n_dim);

    }

    /* Initialize Zi(r) in Algorithm 3.2, Page 20 of Jiao & van Dyk 2015 */
    std::vector<std::vector<double> > zir;
    zir.resize(di.n_dim);
    for(size_t j=0;j<di.n_dim;j++){
      zir[j].resize(di.n_samp);
    }


    std::vector<std::vector<double> > WishMat1, WishMat1Inv, WishSampleTildeInv;

    WishMat1.resize(di.n_dim);
    WishSampleTildeInv.resize(di.n_dim);
    for(size_t j=0;j<di.n_dim;j++){
      WishMat1[j].resize(di.n_dim);
      WishSampleTildeInv[j].resize(di.n_dim);
    }

    //parameters for correction from Jiao & van Dyk 2015
    size_t check_temp;
    size_t max_class_zir = 0;
    double max_zir;
    size_t itercnt;


    double alpha2, alpha2old, ss;
    int sigdrawcounter = 0;

    int dgn = *diagnostics;
    //percA[k][i] = percA[k*m + i] = average over all percAtmp for draw i latent k
    std::vector<std::vector<double> > percAtmp;
    percAtmp.resize(di.n_dim);
    //numNodes[k][i] =  numNodes[k*m + i] = average over all numNodestmp for draw i latent k
    std::vector<std::vector<double> > numNodestmp;
    numNodestmp.resize(di.n_dim);
    //numLeaves[k][i] = numLeaves[k*m + i] = average over all numLeavestmp for draw i latent k
    std::vector<std::vector<double> > numLeavestmp;
    numLeavestmp.resize(di.n_dim);
    //treeDepth[k][i] = treeDepth[k*m + i] = average over all treeDepthtmp for draw i latent k
    std::vector<std::vector<double> > treeDepthtmp;
    treeDepthtmp.resize(di.n_dim);
    //nLtD[k][i] = nLtD[k*m + i] = number of leaves at tree depth of the ith tree for the current round of draw of latent k
    std::vector<std::vector<double> > nLtDtmp;
    nLtDtmp.resize(di.n_dim);


    for(size_t i=0;i<di.n_dim;i++) {
      percAtmp[i].resize(m);
      numNodestmp[i].resize(m);
      numLeavestmp[i].resize(m);
      treeDepthtmp[i].resize(m);
      nLtDtmp[i].resize(m);
    }

    for(size_t k=0; k<di.n_dim; k++){
      for(size_t i=0;i<m;i++) {
        percAtmp[k][i] = 0.0;
        numNodestmp[k][i] = 1.0;
        numLeavestmp[k][i] = 1.0;
        treeDepthtmp[k][i] = 0.0;
        nLtDtmp[k][i] = 1.0;
      }
    }

    //incProptmp[j] = count occurence of xj in splitting rules for the current round of draw
    //incProp[j] = average of incProptmp[j] over all draws
    std::vector<std::vector<double> > incProptmp;
    incProptmp.resize(di.n_dim);

    for(size_t i=0; i<di.n_dim; i++){
      incProptmp[i].resize(di.n_cov);
    }

    for(size_t k=0; k<di.n_dim; k++){
      for(size_t i=0;i<di.n_cov;i++) {
        incProptmp[k][i] = 0.0;
      }
    }

    size_t countdgn = 0;

    //MCMC


    //cout << "\nMCMC:\n";
    time_t tp;
    int time1 = time(&tp);

    /* Initialize counters for outputs vec_class_pred_test and vec_class_pred_train */
    int countvectest = 0;
    int countvectrain = 0;


    for(size_t loop=0;loop<(nd+burn);loop++) { /* Start posterior draws */


    if(loop%100==0) Rprintf("\n iteration: %d of %d \n",loop, nd+burn);

    /* Step 1 (a) */
    draww(w, mu, sigmai, &nn,&nndim,y);
    if(loop==0) draww(w, mu, sigmai, &nn,&nndim,y); /* regenerate w for the initial draw */

    /* Step 1 (b) */
    /* mtemp1 = V x inverse(Sigma) */
    ss=0;
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++) {
        mtemp1[j][k]=0;
      }
    }

    for(size_t i=0;i<di.n_dim;i++){
      for(size_t j=0;j<di.n_dim;j++){
        for(size_t k=0;k<di.n_dim;k++){
          mtemp1[j][k]+=V[j*di.n_dim + i]*sigmai[i*di.n_dim + k];
        }
      }
    }
    /* ss = trace(V x inverse(Sigma)) */
    for(size_t j=0;j<di.n_dim;j++) ss+=mtemp1[j][j];
    /* alpha^2 = trace(V x inverse(Sigma)) / rchisq */alpha2=ss/(double)rchisq((double)nu*di.n_dim);
    /* Step 1 (c) */
    for(size_t k=0; k<di.n_dim; k++){
      for(size_t i=0;i<di.n_samp;i++) {
        wtilde[k][i] = sqrt(alpha2) * (w[i*di.n_dim + k]);
      }
    }

    //done sampling alpha2, w

    /* Step 2 */
    /* See tree sampling theory.doc for explanation */
    for(size_t ntree = 0 ; ntree <m; ntree++){
      for(size_t k=0; k<di.n_dim; k++){
        fit(t[ntree][k], XMat, di, xi, ftemp[k]);
        for(size_t i=0;i<di.n_samp;i++) {
          allfit[k][i] -= ftemp[k][i];
          rtemp[k][i] = wtilde[k][i] - allfit[k][i];
        }
      }


      //get pseudo response
      getpseudoresponse(di, ftemp, rtemp, sigmai, r,condsig);
      //condsig[k] is sqrt psi_k



      for(size_t k=0; k<di.n_dim; k++){
        di.y = &r[k][0];
        pi.sigma = sqrt(alpha2) * condsig[k]; //sqrt psi_k tilde

        if(dgn){
          bd1(XMat, t[ntree][k], xi, di, pi, minobsnode, binaryX, &nLtDtmp[k][ntree], &percAtmp[k][ntree], &numNodestmp[k][ntree], &numLeavestmp[k][ntree], &treeDepthtmp[k][ntree], incProptmp[k]);
        } else {
          bd(XMat, t[ntree][k], xi, di, pi, minobsnode, binaryX);
        }

        fit(t[ntree][k], XMat, di, xi, ftemp[k]);
        for(size_t i=0;i<di.n_samp;i++) {
          allfit[k][i] += ftemp[k][i]	;
        }
      }

    }//ntree


    //done sampling (T,M)

    /* Step 3 */

    /* WishMat1 = V+ sum_i Zi*Zi^T*/
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        WishMat1[j][k]=0;
      }
    }

    for(size_t i=0;i<di.n_samp;i++){
      for(size_t j=0;j<di.n_dim;j++){
        for(size_t k=0;k<di.n_dim;k++){
          WishMat1[j][k] += (wtilde[j][i]-allfit[j][i])* (wtilde[k][i] - allfit[k][i]);
        }
      }
    }

    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        WishMat1[j][k] += V[k*di.n_dim + j] ;
      }
    }


    dinv(WishMat1 ,di.n_dim,WishMat1Inv);

    /* keep the old alpha2 from step 1 sampling (w,alpha)*/
    alpha2old = alpha2;

    //correction from Jiao & van Dyk 2015
    check_temp = 0;
    itercnt = 1;

    while(check_temp != di.n_samp && itercnt<= nSigDr){

      /* Step 3 (a) */
      // generate inverse Sigma
      rWish(WishSampleTildeInv, WishMat1Inv, (int)(nu+di.n_samp),(int)di.n_dim);
      // invert to get Sigma
      dinv(WishSampleTildeInv ,di.n_dim,WishSampleTilde);

      // now check condition 10

      /* Step 3 (b) */
      // definition of zir needs new alpha2 based on Sigma, alpha2 = trace(Sigma)/p, y has p+1 levels
      alpha2 = 0;
      for(size_t j=0; j< di.n_dim; j++) alpha2 += (WishSampleTilde[j][j])/double(di.n_dim);

      /* Step 3 (c) */
      // difine zir
      for(size_t i=0;i<di.n_samp;i++){
        for(size_t j=0;j<di.n_dim;j++){
          zir[j][i] = wtilde[j][i]- (1 - sqrt(alpha2/alpha2old))* allfit[j][i]; //see pseudocode for explanation
        }
      }
      // condition 10 should exist for every training sample
      check_temp = 0;// count samples that satisfies
      for(size_t i=0;i<di.n_samp;i++){

        max_zir = R_NegInf;
        for(size_t j=0;j<di.n_dim;j++){
          if(zir[j][i] > max_zir){
            max_zir = zir[j][i];
            max_class_zir = j+1;
          }
        }
        if(max_zir <= 0){
          max_class_zir = (size_t)maxy;
        }
        if((size_t)y[i] == max_class_zir){
          check_temp++;
        }
      }

      itercnt++;
    }//while

    if(itercnt>= nSigDr){
      std::cout << "\n iteration on Sigma reached upper limit\n";
    }
    //correction end

    /* Step 3 (e) and (f) */
    for(size_t i=0; i<di.n_samp; i++){
      for(size_t k=0; k < di.n_dim; k++){

        mu[i*di.n_dim + k] = allfit[k][i]/sqrt(alpha2); //divide allfit this to transform
        w[i*di.n_dim +k] = allfit[k][i]/sqrt(alpha2old) + (wtilde[k][i]-allfit[k][i]) /sqrt(alpha2) ;
      }
    }

    /* Step 3 (d) */
    for(size_t j=0;j<di.n_dim;j++){
      for(size_t k=0;k<di.n_dim;k++){
        sigmai[j*di.n_dim + k] = WishSampleTildeInv[j][k]*alpha2;
        SigmaTmpInv[j][k] = WishSampleTildeInv[j][k]*alpha2;
        if( (savesigma==1) && (loop>=burn) ){
          psigmasample[sigdrawcounter++] = WishSampleTilde[j][k]/alpha2;
        }
      }
    }

    /* Make Predictions */
    if(loop>=burn){

      if(dgn){
        for(size_t k=0;k<di.n_dim;k++){
          for(size_t i=0;i<m;i++) {
            percA[countdgn] += percAtmp[k][i];
            numNodes[countdgn] += numNodestmp[k][i];
            numLeaves[countdgn] += numLeavestmp[k][i];
            treeDepth[countdgn] += treeDepthtmp[k][i];
          }
          percA[countdgn] /= m;
          numNodes[countdgn] /= m;
          numLeaves[countdgn] /= m;
          treeDepth[countdgn] /= m;

          countdgn++;

          double numsplits = 0;
          for(size_t i=0; i<di.n_cov; i++){
            //std::cout<< "loop"<<loop<<" "<<incProptmp[i] <<" \n";
            numsplits += incProptmp[k][i];
          }
          for(size_t i=0; i<di.n_cov; i++){
            incProp[k*di.n_cov + i] += incProptmp[k][i]/numsplits;//proportion of all splitting rules that uses xi
          }
        }//k

      }//dgn

      dinv(SigmaTmpInv ,di.n_dim,SigmaTmp);
      for(size_t k = 0; k <di.n_samp; k++){
        max_temp = R_NegInf;
        for(size_t l=0; l<di.n_dim; l++){
          mvnmean[l] = mu[k*di.n_dim + l];
        }

        rMVN(mvnsample, mvnmean, SigmaTmp,di.n_dim);

        for(int l = 0 ; l < *pn_dim; l++){
          if(mvnsample[l] > max_temp){
            max_temp = mvnsample[l];
            pclass = l+1;
          }
        }
        if(max_temp <=0) {
          pclass = (int)maxy;
        }
        vec_class_pred_train[countvectrain] = pclass;
        countvectrain++;
        //cout << "pclass: " << pclass << endl;
      }//end prediction for train

      if(*testn) {

        for(size_t k=0; k<dip.n_dim; k++){
          for(size_t i=0;i<dip.n_samp;i++) {
            ppredmeanvec[k][i] = 0.0;

          }
        }

        for(size_t l = 0; l < dip.n_dim; l++){
          for(size_t j=0;j<m;j++) {
            fit(t[j][l], testXMat[loop-burn], dip, xi, fpredtemp);
            for(size_t k=0;k<dip.n_samp;k++) ppredmeanvec[l][k] += fpredtemp[k];
          }
        }

        for(size_t k = 0; k <dip.n_samp; k++){
          max_temp = R_NegInf;

          for(size_t l=0; l<di.n_dim; l++){
            mvnmean[l] = ppredmeanvec[l][k]/sqrt(alpha2) ;
          }

          rMVN(mvnsample, mvnmean, SigmaTmp,di.n_dim);

          for(int l = 0 ; l < *pn_dim; l++){
            if(mvnmean[l] > max_temp){
              max_temp = mvnmean[l];
              pclass = l+1;
            }
          }
          if(max_temp <=0) {
            pclass = (int)maxy;
          }
          vec_class_pred_test[countvectest] = pclass;
          countvectest++;
          //cout << "pclass: " << pclass << endl;
        }//end prediction for test

      }//end if test

    }//end prediction for current loop



    } //end of loop

    if(dgn){
      for(size_t i=0; i<(di.n_dim*di.n_cov); i++){
        incProp[i] /= nd;
      }
    }

    int time2 = time(&tp);
    Rprintf("time for mcmc loop %d secs", time2-time1);

    for(size_t i=0; i<m; i++){
      for(size_t k=0;k< di.n_dim;k++) t[i][k].tonull();
    }//delete trees
  }
};
