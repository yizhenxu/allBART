#'Probit Bayesian Additive Regression Trees
#'
#'Bayesian Additive Regression Trees Modeling for Binary Outcome,
#'@param formula response ~ covariates,
#'@param train.data Training Data with the response taking values 0-1,
#'@param test.data Test Data, typically without the response,
#'@param Prior List of Priors for MPBART: e.g., Prior = list(ntrees=200,  kfac=2.0, pswap=0,  pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, priorindep = FALSE,  minobsnode = 10).
#'The components of Prior are
#' \itemize{
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
#'@param Mcmc List of MCMC starting values, burn-in ...: e.g.,     list(burn = 100, ndraws = 1000)
#'@return treefit_train ndraws x n posterior matrix of the training data sum of trees fit,
#'@return treefit_test ndraws x testn posterior matrix of the test data sum of trees fit,
#'@return samp_train ndraws x n posterior matrix of the training data outcome,
#'@return samp_test ndraws x testn posterior matrix of the test data outcome.
#'@import bayesm mlbench mlogit cvTools stats
#'@export
#'@useDynLib allBART
PBART_call  <- function(formula, data, test.data = NULL,
                        Prior = NULL, Mcmc = NULL)
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

  ####################################

  Terms <- attr(mf, "terms")

  X <- model.matrix.default(Terms, mf)

  xcolnames <- colnames(X)[-1]

  if(length(xcolnames) == 1 ){
    X <- data.frame(X[,xcolnames])
    names(X) <- xcolnames[1]

  } else {

    X <- X[,xcolnames]
  }


  if(length(xcolnames) == 1 ){
    Xtest <- data.frame(test.data[,xcolnames])
    names(Xtest) <- xcolnames[1]
  } else {
    Xtest <- test.data[,xcolnames]
  }


  Data = list(y = Y,X = X)
  testData = list(X = Xtest)

  n=length(Data$y)

  if (!is.null(test.data)){
    testn <- nrow(testData$X)
  } else {
    testn <- 0
  }

  binaryX = rep(NA,ncol(X))
  for(i in 1:ncol(X)){
    binaryX[i] = 1*(length(unique(X[,i]))==2)
  }


  if(missing(Prior))
  {ntrees=200; kfac=2.0;pswap=0;pbd=1.0;pb=0.5;beta = 2.0;alpha = 0.95; nc = 100; priorindep = 0; minobsnode = 10;
  }
  else
  { if(is.null(Prior$ntrees)) {ntrees=200} else {ntrees=Prior$ntrees}
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


  if(is.null(Mcmc$burn)) {burn=100} else {burn=Mcmc$burn}
  if(is.null(Mcmc$ndraws)) {ndraws=1000} else {ndraws=Mcmc$ndraws}

  cat("Number of trees: ", ntrees, ".\n\n", sep="")
  cat("Number of draws: ", ndraws, ".\n\n", sep="")
  cat("burn-in: ", burn, ".\n\n", sep="")



  res =   .C('mypbart',
             trainx= as.double(t(Data$X)),
             testx= as.double(t(testData$X)),
             mu = as.double(rep(0,nrow(Data$X))),
             n = as.integer(length(Data$y)),
             y = as.double(Data$y),
             n_cov = as.integer(ncol(Data$X)),
             testn = as.integer(testn),
             ndraws = as.integer(ndraws),
             burn = as.integer(burn),
             ntrees = as.integer(ntrees),
             kfac = as.double(kfac),
             pswap = as.double(pswap),
             pbd = as.double(pbd),
             pb = as.double(pb),
             alpha = as.double(alpha),
             beta =  as.double(beta),
             nc = as.integer(nc),
             minobsnode = as.integer(minobsnode),
             vec_test = as.double(rep(0,testn*ndraws)),
             vec_train = as.double(rep(0,n*ndraws)),
             binaryX = as.integer(binaryX))


  #vec_test and vec_train are the sum of tree fits

  yhat.train = res$vec_train #nsamp X ndraws
  vec_samp_train = rnorm(n*ndraws, yhat.train, 1)
  vec_samp_train = 1*(vec_samp_train > 0)
  vec_samp_train = matrix(vec_samp_train, byrow = TRUE, ncol = n) #ndraws X nsamp

  if (!is.null(test.data)){
    yhat.test = res$vec_test
    vec_samp_test = rnorm(testn*ndraws, yhat.test, 1)
    vec_samp_test = 1*(vec_samp_test > 0)
    vec_samp_test = matrix(vec_samp_test, byrow = TRUE, ncol = testn)

  } else {
    vec_samp_test <- NULL
  }

  ret = list(treefit_test = matrix(yhat.test, byrow = TRUE, ncol = n),
             treefit_train = matrix(yhat.train, byrow = TRUE, ncol = testn),
             samp_test = vec_samp_test,
             samp_train = vec_samp_train);

  return(ret)
}
