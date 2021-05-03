#' Gaussian Core Vector Regression
#'
#' A method that uses parameters to reweight both samples and features.
#' Ideas based off of Vector Generalized Linear Models. See https://en.wikipedia.org/wiki/Vector_generalized_linear_model
#'
#' @param x Reference dataset of cell by gene matrix
#' @param y Experimental dataset of cell by gene matrix
#' @param covariate_x Experimental covariates of dimension c_x by cell
#' @param covariate_y Reference covariates of dimension c_y by cell
#' @param k_dim Number of k reduced dimensions
#' @param max_iter Maximum number of iterations possible
#' @param min_iter Minimum number of iterations possible
#' @param tol Tolerance value to assess convergence (converges when less than tol)
#' @param log True for log_e(1+...) transform, False for no transform
#' @param factor Initialise with library normalisation to consider total counts per cell library
#' @param select_run Select genes to begin dimension reduction with, otherwise takes the most variable genes
#' @return  A list of two items - a list of test statistics and a list of internal parameters
#' @export
corevec <- function(x, y, covariate_x, covariate_y, k_dim = NULL, max_iter=100, min_iter = 15, tol=1e-5, log=T, factor=T, scale=T, select_run=NULL){

  if (!is.null(select_run) & !is.null(k_dim)){
    print("Please provide one of k_dim or select_run")
    break
  }

  # Prepare for log_e(1+x) transform and library normalisation
  cx <- Matrix::colSums(x)
  cy <- Matrix::colSums(y)

  if (log == T){
    x <- log(1+x)
    y <- log(1+y)
  }
  if (scale==T){
    x <- scale(x)
    y <- scale(y)
  }
  if (factor==T){
    x <- x - Matrix::t(cx%*%MASS::ginv(as.matrix(Matrix::t(cx)%*%cx))%*%Matrix::t(cx)%*%Matrix::t(x))
    y <- y - Matrix::t(cy%*%MASS::ginv(as.matrix(Matrix::t(cy)%*%cy))%*%Matrix::t(cy)%*%Matrix::t(y))
  }

  # Initialise parameters
  param.list <- initialise_param.list(x=x,y=y,covariate_y=covariate_y, covariate_x = covariate_x, alpha.init = NULL, beta.init =  NULL, k_dim = k_dim, select_run = select_run)

  llik<-score<-tol_vec<-beta_vec<-alpha.L_vec<-alpha.L_vec<-c()
  count=0

  # Convert parameters to the Expectation space
  param.list$S.list <- convert_param_to_S(param.list = param.list,alpha_run = T,beta_run = T)

  # Start timing
  a <- Sys.time()
  while (T){

    # Run internal EM algorithm
    param.list <- convert_S_to_param(param.list = param.list, alpha_run = T, beta_run = T)

    # Update meta parameters
    llik <- c(llik,param.list$mllik)
    score <- c(score,param.list$mscore)
    tol_vec <- c(tol_vec,abs(tail(llik,2)[1]-tail(llik,1)[1]))

    # Check convergence
    if (count>min_iter){
      if (tail(llik,2)[2]<tail(llik,2)[1]){
        break
      }
      if (count>max_iter | tail(tol_vec,1)<tol){
        break
      }
    }

    # Print updates
    if (count %% 5 == 0){
      b <- Sys.time()
      print(paste("Iterations ", count," in ", round(b-a,3)," with tol at",round(tail(tol_vec,1),20)))
    }

    count = count + 1

  }

  # Statistics for interpretation and analysis
  statistics <- extract_statistics(param.list = param.list)

  statistics$iter = count
  statistics$max_iter = max_iter
  statistics$tol_vec = tol_vec
  statistics$llik = llik
  statistics$score = score

  return(list(statistics = statistics,
              internal_parameters = param.list))
}


extract_statistics <- function(param.list){

  x <- (param.list$x)
  y <- (param.list$y)
  covariate_y <- (param.list$covariate_y)
  covariate_x <- (param.list$covariate_x)
  y_hat_alpha_beta <- (param.list$y_hat_alpha_beta)
  alpha.L <- (param.list$alpha.L)
  beta <- (param.list$beta)
  delta <- (param.list$delta)
  gamma = (param.list$gamma)

  return(list(
    x = as.matrix((x)),
    y = as.matrix((y)),
    covariate_x = as.matrix(covariate_x),
    covariate_y = as.matrix(covariate_y),
    fitted = as.matrix((y_hat_alpha_beta)),
    residuals = as.matrix((as.matrix((y)) - (y_hat_alpha_beta))),
    alpha.L = as.matrix(alpha.L),
    beta = as.matrix(beta),
    delta = as.matrix(delta),
    gamma = as.matrix(gamma)
  )
  )
}

initialise_param.list <- function(x,y,covariate_y,covariate_x,alpha.init=NULL,beta.init=NULL, select_run=F, k_dim = 10){
  # Prepare y
  yT <- Matrix::t(y)
  y <- Matrix::Matrix(y)

  # Prepare x
  xT <- Matrix::t(x)
  x <- Matrix::Matrix(x)

  if (is.null(covariate_y)){
    covariate_y = as.vector(rep(0,dim(y)[2]))
  }
  if (is.vector(covariate_y)){
    n = 1
    p = length(covariate_y)
  }else{
    n = dim(covariate_y)[1]
    p = dim(covariate_y)[2]
  }

  covariate_y <- Matrix::Matrix(matrix(covariate_y,nrow = n,ncol=p))

  if (is.null(covariate_x)){
    covariate_x = as.vector(rep(0,dim(x)[2]))
  }
  if (is.vector(covariate_x)){
    n = 1
    p = length(covariate_x)
  }else{
    n = dim(covariate_x)[1]
    p = dim(covariate_x)[2]
  }

  covariate_x <- Matrix::Matrix(matrix(covariate_x,nrow = n,ncol=p))

  if (is.null(beta.init)){
    beta = MASS::ginv(as.matrix(Matrix::t(x)%*%(x)))%*%Matrix::t((x))%*%(y)
  } else {
    beta = beta.init
  }

  if (is.null(alpha.init)){
    alpha = (y)%*%Matrix::t(beta)%*%MASS::ginv(as.matrix(beta%*%Matrix::t(beta)))%*%xT%*%MASS::ginv(as.matrix(x%*%xT))
  } else {
    alpha = alpha.init
  }


  if (!is.null(select_run)){
    alpha.L <- alpha[which(row.names(x)%in%c(select_run)),]
  } else {

    alpha.L <- alpha[order(apply(alpha,1,var),decreasing = T)[1:k_dim],]
  }

  delta = (alpha.L%*%y - alpha.L%*%x%*%beta)%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))

  gamma = (alpha.L%*%y - alpha.L%*%x%*%beta - delta%*%covariate_y)%*%Matrix::t(covariate_x%*%beta)%*%MASS::ginv(as.matrix((covariate_x%*%beta)%*%Matrix::t(covariate_x%*%beta)))

  param.list <- list(x = x,
                     y = y,
                     alpha.L = alpha.L,
                     alpha.L = alpha.L,
                     beta = beta,
                     delta = delta,
                     gamma = gamma,
                     covariate_x = covariate_x,
                     covariate_y = covariate_y)

  return(param.list)

}

convert_param_to_S <- function(param.list,alpha_run=T,beta_run=T){

  # Prepare main parameters
  x <- param.list$x
  y <- param.list$y
  covariate_y <- param.list$covariate_y
  alpha.L <- param.list$alpha.L
  alpha.L <- param.list$alpha.L
  beta <- param.list$beta
  delta <- param.list$delta

  main_beta <- c()
  main_alpha <- c()

  # E-step of EM algorithm for beta parameter
  if (beta_run==T){

    alphax <- alpha.L%*%x

    t.alphax_y <- Matrix::t(alphax)%*%alpha.L%*%y

    t.alphax_alphax <- Matrix::t(alphax)%*%alphax


    main_alpha <- list(
      alphax = alphax,
      t.alphax_y = t.alphax_y,
      t.alphax_alphax = t.alphax_alphax
    )
  }

  # E-step of EM algorithm for L parameter
  if (alpha_run==T){

    xbeta <- x%*%beta

    y_t.xbeta <- alpha.L%*%y%*%Matrix::t(xbeta)

    xbeta_t.xbeta <- xbeta%*%Matrix::t(xbeta)

    main_beta <- list(
      xbeta = xbeta,
      y_t.xbeta = y_t.xbeta,
      xbeta_t.xbeta = xbeta_t.xbeta
    )
  }

  return(c(main_alpha,main_beta))
}

convert_S_to_param <- function(param.list, alpha_run=T, beta_run=T){

  # Prepare main parameters
  S.list <- param.list$S.list

  y <- param.list$y
  x <- param.list$x
  beta <- param.list$beta
  alpha.L <- param.list$alpha.L
  alpha.L <- param.list$alpha.L
  covariate_x <- param.list$covariate_x
  covariate_y <- param.list$covariate_y
  delta <- param.list$delta
  gamma <- param.list$gamma

  if (alpha_run){
    # Calculate parts of L parameter
    S.1_alpha <-
      S.list$y_t.xbeta +
      delta%*%covariate_y%*%Matrix::t(S.list$xbeta) -
      gamma%*%covariate_x%*%beta%*%Matrix::t(S.list$xbeta)

    S.2_alpha <-
      S.list$xbeta_t.xbeta

    # Calculate L parameter
    param.list$alpha.L <- alpha.L <- Matrix::Matrix(S.1_alpha%*%MASS::ginv(as.matrix(S.2_alpha)),sparse=T)

    # Update in expectation space
    S.list <- convert_param_to_S(param.list = param.list,alpha_run = alpha_run, beta_run = beta_run)
  }


  if (beta_run){
    # Calculate parts of beta parameter
    S.1_beta <-
      S.list$t.alphax_y +
      Matrix::t(gamma%*%covariate_x)%*%alpha.L%*%y -
      Matrix::t(S.list$alphax)%*%delta%*%covariate_y -
      Matrix::t(gamma%*%covariate_x)%*%delta%*%covariate_y

    S.2_beta <-
      S.list$t.alphax_alphax +
      Matrix::t(gamma%*%covariate_x)%*%S.list$alphax +
      Matrix::t(S.list$alphax)%*%(gamma%*%covariate_x) +
      Matrix::t(gamma%*%covariate_x)%*%(gamma%*%covariate_x)

    # Calculate beta parameter
    beta <- Matrix::Matrix(MASS::ginv(as.matrix(S.2_beta))%*%S.1_beta)
    param.list$beta <- beta
  }

  delta = (alpha.L%*%y - alpha.L%*%(x)%*%beta - gamma%*%covariate_x%*%beta )%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))

  gamma = (alpha.L%*%y - alpha.L%*%(x)%*%beta - delta%*%covariate_y)%*%Matrix::t(covariate_x%*%beta)%*%MASS::ginv(as.matrix((covariate_x%*%beta)%*%Matrix::t(covariate_x%*%beta)))

  param.list$delta <- delta
  param.list$gamma <- gamma

  # Update in expectation space
  S.list <- convert_param_to_S(param.list,alpha_run = alpha_run,beta_run = beta_run)

  # Calculate fitted value of y; covariance of y
  y_hat_alpha_beta <- as.matrix(x)%*%beta + MASS::ginv(as.matrix(Matrix::t(alpha.L)%*%alpha.L))%*%Matrix::t(alpha.L)%*%(delta%*%covariate_y + gamma%*%covariate_x%*%beta)
  D1 = Matrix::diag(Matrix::diag(Matrix::t((y) - (y_hat_alpha_beta))%*%((y) - (y_hat_alpha_beta))/dim(y)[1]))

  # Calculate likelihood ; mean-squared error
  mllik <- sum(mclust::dmvnorm(y - y_hat_alpha_beta,rep(0,dim(D1)[1]),D1,log = T))
  mscore <- mean(sqrt(as.matrix(alpha.L%*%as.matrix(y - y_hat_alpha_beta))^2))

  return(list(
    S.list = S.list,
    x = x,
    y = y,
    covariate_x = covariate_x,
    covariate_y = covariate_y,
    y_hat_alpha_beta = y_hat_alpha_beta,
    alpha.L = alpha.L,
    beta = beta,
    delta = delta,
    gamma = gamma,
    mllik = mllik,
    mscore = mscore
  ))
}
