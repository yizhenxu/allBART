#'Multinomial Probit Bayesian Additive Regression Trees
#'
#'Multinomial probit modeling using Bayesian Additive Regression Trees,
#'@param formula response ~ demographic covariates. Demographic covariates are the covariates that are not specific to levels of the multinomial response,
#'@param train.data Training Data in wide format (for details on wide format, see documentation in R package \pkg{mlogit}). Names of alternative specific covariates (the covariates that are specific to levels of the multinomial response) should in the format of “Xa.y”, where Xa is a variable name, and y is a level of the multinomial response, and “.” is the separator. For examle, if “price” is an alternative specific covariate, and the response takes values 2, 3, and 4, then the train.data should contain “price.2”, “price.3”, and “price.4”,
#'@param test.data Test Data in wide format, typically without the response,
#'@param base order index of the reference level of the multinomial response. For example, if the response takes values 2, 3, and 4, then base = 2 sets response value 3 as the reference. Default is the highest class,
#'@param varying The indeces of the variables in the train.data that are alternative specific. The length of varying should be a multiple of the number of levels in the multinomial response,
#'@param sep The seperator of the variable name and the alternative level in the alternative specific covariates. The default separator is dot (.).
#'@param Prior List of Priors for MPBART: e.g., Prior = list(nu=p+2,  V= diag(p - 1), ntrees=200,  kfac=2.0, pswap=0,  pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, priorindep = FALSE,  minobsnode = 10).
#'The components of Prior are
#' \itemize{
#'\item nu : The covariance matrix of latent variables is assumed to have prior distribution Inv-Wish(nu, V). nu is the degree of freedom and nu > (nlatent - 1).
#'\item V : The positive definite scale matrix in the Inverse-Wishart prior of the covariance matrix.
#'\item ntrees : The total number of trees in each round of BART fitting.
#'\item kfac : A tuning parameter that satisfies mu - kfac * sigma = ymin and mu + kfac * sigma = ymax, where mu and sigma are the mean and std of the Gaussian prior distribution of the sum of fit of all trees.
#'\item pswap : The prior probability of swap move in simulating trees; default 0, there should be pswap no larger than 1-pbd.
#'\item pbd : The prior probability of birth/death move in simulating trees; default 1.
#'\item pb : The prior probability of birth given birth/death; default 0.5.
#'\item alpha : The prior probability of a bottom node splits is alpha/(1+d)^beta, d is depth of node.
#'\item beta : see alpha.
#'\item nc : The number of equally spaced cutpoints between min and max for each covariate.
#'\item priorindep
#'\item minobsnode : The minimum number of observations in bottom nodes for birth in simulating trees.
#'}
#'@param Mcmc List of MCMC starting values, burn-in ...: e.g.,     list(sigma0 = diag(p - 1), burn = 100, ndraws = 1000, nSigDr = 50, keep_sigma_draws=FALSE)
#'The components of Mcmc are
#' \itemize{
#'\item sigma0 : The starting value of the covariance matrix of latent variables.
#'\item nSigDr: User-specified upper limit to repeated draws of the covariance variance matrix of latent variables in each round of posterior draw when condition 10 in Jiao and van Dyk 2015 is not satisfied. Default 50.
#'}
#'@return samp_train ndraws x n posterior matrix of the training data outcome,
#'@return samp_test ndraws x testn posterior matrix of the test data outcome,
#'@return sigmasample posterior samples of the latent variable covariance matrix.
#'@import bayesm mlbench mlogit cvTools stats
#'@export
#'@useDynLib allBART
MPBART_call  <- function(formula, data, base = NULL,test.data = NULL,
                         Prior = NULL, Mcmc = NULL,
                         varying = NULL, sep = '.')
{

  callT <- match.call(expand.dots = TRUE)

  formula <- mFormula(formula)


  response.name <- paste(deparse(attr(formula, "lhs")[[1L]]))
  m <- match(c("formula", "data"), names(callT), 0L)
  mf <- callT
  mf <- mf[c(1L, m)]

  mf$formula <- formula


  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf)
  Y <- as.factor(Y)
  lev <- levels(Y)
  p <- length(lev)
  if (p < 3){
    stop("The number of alternatives should be at least 3.")
  }

  if (!is.null(base))
  {
    base <- base
  } else {
    base <- lev[p]
  }

  if (base <= p) {
    Y <- relevel(Y, ref = base)
    lev <- levels(Y)
  } else {
    stop(paste("Error: `base' does not exist in the response variable."))
  }

  base <- lev[1]
  counts <- table(Y)
  if (any(counts == 0)) {
    warning(paste("group(s)", paste(lev[counts == 0], collapse = " "), "are empty"))
    Y <- factor(Y, levels  = lev[counts > 0])
    lev <- lev[counts > 0]
  }
  Y <- as.numeric(unclass(Y)) - 1
  Y <- ifelse(Y==0, p,Y)

  cat("The base level is: '", lev[1], "'.\n\n", sep="")

  relvelved <- c(lev[2:length(lev)], lev[1])

  Terms <- attr(mf, "terms")

  X <- model.matrix.default(Terms, mf)



  xcolnames <- colnames(X)[-1]

  if(length(xcolnames) == 1 ){
    X <- data.frame(X[,xcolnames])
    names(X) <- xcolnames[1]

  } else {

    X <- X[,xcolnames]
  }
  if (!is.null(test.data)){

    if(length(xcolnames) == 1 ){
      Xtest <- data.frame(test.data[,xcolnames])
      names(Xtest) <- xcolnames[1]
    } else {
      Xtest <- test.data[,xcolnames]
    }

    testXEx = NULL;
    for(i in 1:nrow(Xtest)){
      testXEx = rbind(testXEx, matrix(rep(Xtest[i,], p-1), byrow = TRUE, ncol = ncol(Xtest) ) )
    }


  } else {
    testXEx = 0
  }


  XEx = NULL;
  for(i in 1:nrow(X)){
    XEx = rbind(XEx, matrix(rep(X[i,], p-1), byrow = TRUE, ncol = ncol(X) ) )
  }


  Data = list(p=p,y=Y,X= XEx)
  testData = list(p=p,X= testXEx)


  cat("Table of y values",fill=TRUE)
  print(table(model.response(mf) ))


  n=length(Data$y)

  pm1=p-1

  if (!is.null(test.data)){
    testn <- nrow(testData$X)/(p-1)
  } else {
    testn <- 0
  }


  #reading alternate specific variables
  if (!is.null(varying)) {

    varying.names <-    names(data)[varying]

    alt.names <- NULL
    for(vv in 1:(length(varying))){
      alt.names <- c(alt.names, unlist(strsplit(varying.names[vv], sep, fixed = TRUE))[1L])
    }



    alt.names <- unique(alt.names)


    if(length(alt.names) != length(varying)/p){
      stop("alternative variables names mismatch. Check names of alternative variables.")
    }

    reordered_names <- NULL
    for(nn in 1:(length(varying)/p)){
      reordered_names <- c(reordered_names,  paste0(alt.names[nn], sep, relvelved) )
    }

    XChSp <- createX(p,na=length(varying)/p ,
                     nd=NULL,Xa= data[,reordered_names],Xd=NULL,
                     INT = FALSE,DIFF=TRUE,base=p)

    Data$X <- cbind(Data$X, XChSp)

    if (!is.null(test.data)){
      testXChSp = createX(p,na=length(varying)/p ,
                          nd=NULL,Xa= test.data[,reordered_names],Xd=NULL,
                          INT = FALSE,DIFF=TRUE,base=p)
      testData$X <- cbind(testData$X, testXChSp)
    }

  }


  binaryX = rep(NA,ncol(Data$X))
  for(i in 1:ncol(Data$X)){
    binaryX[i] = 1*(length(unique(Data$X[,i]))==2)
  }


  if(missing(Prior))
  {nu=pm1+3; V=nu*diag(pm1);
  ntrees=200; kfac=2.0;pswap=0;pbd=1.0;pb=0.5;beta = 2.0;alpha = 0.95; nc = 100; priorindep = 0; minobsnode = 10;
  }
  else
  {if(is.null(Prior$nu)) {nu=pm1+3} else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(pm1)} else {V=Prior$V}
    if(is.null(Prior$ntrees)) {ntrees=200} else {ntrees=Prior$ntrees}
    if(is.null(Prior$kfac)) {kfac=2.0} else {kfac=Prior$kfac}
    if(is.null(Prior$pswap)) {pswap=0.0} else {pswap=Prior$pswap}
    if(is.null(Prior$pbd)) {pbd=1.0} else {pbd=Prior$pbd}
    if(is.null(Prior$pb)) {pb=0.5} else {pb=Prior$pb}
    if(is.null(Prior$beta)) {beta = 2.0} else {beta=Prior$beta}
    if(is.null(Prior$alpha)) {alpha = 0.95} else {alpha=Prior$alpha}
    if(is.null(Prior$nc)) {nc=100} else {nc=Prior$nc}
    if(is.null(Prior$priorindep)) {priorindep= FALSE} else {priorindep=Prior$priorindep}
    if(is.null(Prior$minobsnode)) {minobsnode= 10} else {minobsnode=Prior$minobsnode}


  }

  if(is.null(Mcmc$sigma0)) {sigma0=diag(pm1)} else {sigma0=Mcmc$sigma0}

  if(is.null(Mcmc$burn)) {burn=100} else {burn=Mcmc$burn}
  if(is.null(Mcmc$ndraws)) {ndraws=1000} else {ndraws=Mcmc$ndraws}
  if(is.null(Mcmc$nSigDr)) {nSigDr=50} else {nSigDr=Mcmc$nSigDr}
  if(is.null(Mcmc$keep_sigma_draws)) {keep_sigma_draws=FALSE} else {keep_sigma_draws=Mcmc$keep_sigma_draws}

  sigmai = solve(sigma0)
  if( (priorindep ==TRUE) || (keep_sigma_draws==FALSE)){
    sigmasample = as.double(0);
    savesigma = 0;
  } else {
    sigmasample = as.double(rep(sigma0, ndraws));
    savesigma = 1;
  }

  cat("Number of trees: ", ntrees, ".\n\n", sep="")
  cat("Number of draws: ", ndraws, ".\n\n", sep="")
  cat("burn-in: ", burn, "'\n\n", sep="")



  res =   .C('mympbart',w=as.double(rep(0,nrow(Data$X))),
             trainx= as.double(t(Data$X)),
             testx= as.double(t(testData$X)),
             mu = as.double(rep(0,nrow(Data$X))),
             sigmai = as.double(sigmai),
             V = as.double(V),
             n = as.integer(length(Data$y)),
             n_dim = as.integer(ncol(sigmai)),
             y = as.integer(Data$y),
             n_cov = as.integer(ncol(Data$X)),
             nu = as.integer(nu),
             testn = as.integer(testn),
             ndraws = as.integer(ndraws),
             burn = as.integer(burn),
             ntrees = as.integer(ntrees),
             nSigDr = as.integer(nSigDr),
             kfac = as.double(kfac),
             pswap = as.double(pswap),
             pbd = as.double(pbd),
             pb = as.double(pb),
             alpha = as.double(alpha),
             beta =  as.double(beta),
             nc = as.integer(nc),
             savesigma = as.integer(savesigma),
             minobsnode = as.integer(minobsnode),
             sigmasample = sigmasample,
             vec_test = as.integer(rep(0,testn*ndraws)),
             vec_train = as.integer(rep(0,n*ndraws)),
             binaryX = as.integer(binaryX))

  vec_class_train <- matrix(relvelved[res$vec_train], byrow = TRUE, ncol = n)

  if (!is.null(test.data)){

    vec_class_test <- matrix(relvelved[res$vec_test], byrow = TRUE, ncol = testn)

  } else {
    vec_class_test <- NULL
  }

  if(savesigma){
    ret = list(samp_test = vec_class_test,
               samp_train = vec_class_train,
               sigmasample = res$sigmasample);
  } else {
    ret = list(samp_test = vec_class_test,
               samp_train = vec_class_train);
  }
  return(ret)
}
