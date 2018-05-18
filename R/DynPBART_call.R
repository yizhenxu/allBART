#'Probit Bayesian Additive Regression Trees
#'
#'Bayesian Additive Regression Trees Modeling for Binary Outcome for Dynamic Data,
#'@param formula response ~ covariates,
#'@param data Training data with the response taking values 0-1,
#'@param test.data Test data with number of rows equals nsamp x ndraws, typically without the response. For the ith subject, (i, i+nsamp, ..., i+(ndraws-1)nsamp) rows are its simulated posterior outcomes from the previous simulation,
#'@param Prior List of Priors for MPBART: e.g., Prior = list(ntrees=200,  kfac=2.0, pswap=0,  pbd=1.0, pb=0.5 , beta = 2.0, alpha = 0.95, nc = 100, minobsnode = 10).
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
#'\item minobsnode : The minimum number of observations in bottom nodes for birth in simulating trees.
#'}
#'@param Mcmc List of MCMC starting values, burn-in ...: e.g.,     list(burn = 100, ndraws = 1000)
#'@param diagnostics Returns convergence diagnostics and variable inclusion proportions if True (default),
#'@return treefit_train ndraws x n posterior matrix of the training data sum of trees fit,
#'@return treefit_test ndraws x testn posterior matrix of the test data sum of trees fit,
#'@return samp_train ndraws x n posterior matrix of the training data outcome,
#'@return samp_test ndraws x testn posterior matrix of the test data outcome,
#'@return Percent_Acceptance Percent acceptance of Metropolis-Hastings proposals across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Num_Nodes Average number of tree nodes across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Num_Leaves Average number of leaves across the ntrees number of trees for each posterior draw after burn-in,
#'@return Tree_Depth Average tree depth across the ntrees number of trees for each posterior draw after burn-in,
#'@return Inclusion_Proportions Predictor inclusion frequencies. Smaller value of ntrees (such as 10, 20, 50, 100) is recommended for the purposes of variable selection.
#'@examples
#'##simulate data (example from Friedman MARS paper)
#'f = function(x){
#'  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#'}
#'sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
#'n = 100 #number of observations
#'set.seed(99)
#'###############################################
#'nd = 1000 # ndraws
#'x1 = matrix(runif(n*9),n,9) #10 variables, only first 5 matter
#'x1 = x1[rep(1:n,nd),]
#'x2 = runif(n*nd)
#'x = cbind(x1, x2)
#'Ey = f(x)
#'u=Ey+sigma*rnorm(n*nd)
#'z = (u-min(u))/(max(u)-min(u))
#'y = 1*(z<0.4)+ 2*(z>=0.4 & z<0.6) + 3*(z>=0.6)
#'tdat = data.frame(x,y)
#'###############################################
#'x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
#'Ey = f(x)
#'u=Ey+sigma*rnorm(n)
#'z = (u-min(u))/(max(u)-min(u))
#'y= rbinom(100, 1, z)
#'dat = data.frame(x,y)
#'colnames(tdat) = colnames(dat)
#'###############################################
#'fml = as.formula("y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10")
#'bpfit = DynPBART_call(fml, data = dat, test.data = tdat,
#'                 Prior = list(ntrees = 100,
#'                              kfac = 2,
#'                              pswap = 0.1, pbd = 0.5, pb = 0.25,
#'                              alpha = 0.95, beta = 2.0,
#'                              nc = 100, minobsnode = 10),
#'                 Mcmc = list(burn=100, ndraws = nd))
#'
#'@import bayesm mlbench mlogit cvTools stats
#'@export
#'@useDynLib allBART
DynPBART_call  <- function(formula, data, test.data = NULL,
                        Prior = NULL, Mcmc = NULL, diagnostics = TRUE)
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

  if (!is.null(test.data)){

    if(is.null(Mcmc$ndraws))
      stop(paste("Error: ndraws is required for dynamic model calls."))

    ndraws=Mcmc$ndraws
    testnsub = nrow(test.data) / ndraws

    if(testnsub %% 1 != 0)
      stop(paste("Error: testn does not equal to test.nsub x ndraws."))

    if(length(xcolnames) == 1 ){
      Xtest <- data.frame(test.data[,xcolnames])
      names(Xtest) <- xcolnames[1]
    } else {
      Xtest <- test.data[,xcolnames]
    }
  } else {
    stop(paste("Error: no test data, please use PBART_call instead."))
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
  {ntrees=200; kfac=2.0;pswap=0;pbd=1.0;pb=0.5;beta = 2.0;alpha = 0.95; nc = 100; minobsnode = 10;
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
    if(is.null(Prior$minobsnode)) {minobsnode= 10} else {minobsnode=Prior$minobsnode}
  }


  if(is.null(Mcmc$burn)) {burn=100} else {burn=Mcmc$burn}
  #if(is.null(Mcmc$ndraws)) {ndraws=1000} else {ndraws=Mcmc$ndraws}

  cat("Number of trees: ", ntrees, ".\n\n", sep="")
  cat("Number of draws: ", ndraws, ".\n\n", sep="")
  cat("burn-in: ", burn, ".\n\n", sep="")



  res =   .C('mydynpbart',
             trainx= as.double(t(Data$X)),
             testx= as.double(t(testData$X)),
             mu = as.double(rep(0,nrow(Data$X))),
             n = as.integer(length(Data$y)),
             y = as.double(Data$y),
             n_cov = as.integer(ncol(Data$X)),
             testn = as.integer(testn),
             testnsub = as.integer(testnsub),
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
             vec_test = as.double(rep(0,testn)),
             vec_train = as.double(rep(0,n*ndraws)),
             binaryX = as.integer(binaryX),
             diagnostics = as.integer(diagnostics),
             percA = as.double(rep(0, ndraws)),
             numNodes = as.double(rep(0, ndraws)),
             numLeaves = as.double(rep(0, ndraws)),
             treeDepth = as.double(rep(0, ndraws)),
             incProp = as.double(rep(0, ncol(Data$X))) )

  names(res$incProp) = colnames(Data$X)

  #vec_test and vec_train are the sum of tree fits

  yhat.train = res$vec_train #nsamp X ndraws
  vec_samp_train = rnorm(n*ndraws, yhat.train, 1)
  vec_samp_train = 1*(vec_samp_train > 0)
  vec_samp_train = matrix(vec_samp_train, nrow = n) #nsamp X ndraws

  yhat.train = matrix(yhat.train, nrow = n)

  if (!is.null(test.data)){
    yhat.test = res$vec_test
    vec_samp_test = rnorm(testn, yhat.test, 1)
    vec_samp_test = 1*(vec_samp_test > 0)
    vec_samp_test = matrix(vec_samp_test, nrow = testnsub)
    yhat.test = matrix(yhat.test, nrow = testnsub)

  } else {
    vec_samp_test <- NULL
    yhat.test <- NULL
  }

  if(diagnostics){
    ret = list(treefit_test = yhat.test,
               treefit_train = yhat.train,
               samp_test = vec_samp_test,
               samp_train = vec_samp_train,
               Percent_Acceptance = res$percA,
               Tree_Num_Nodes = res$numNodes,
               Tree_Num_Leaves = res$numLeaves,
               Tree_Depth = res$treeDepth,
               Inclusion_Proportions = res$incProp);
  } else {
    ret = list(treefit_test = yhat.test,
               treefit_train = yhat.train,
               samp_test = vec_samp_test,
               samp_train = vec_samp_train);
  }

  return(ret)
}
