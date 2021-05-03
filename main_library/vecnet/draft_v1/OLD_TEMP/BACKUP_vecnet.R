#' Gaussian Vector Network Regression
#'
#' Gaussian Vector Network Regression is a method that uses parameters to reweight both samples and features.
#'
#' @param x Reference dataset of cell by gene matrix
#' @param y Experimental dataset of cell by gene matrix
#' @param max_iter Maximum number of iterations possible
#' @param tol Tolerance value to assess convergence (converges when less than tol)
#' @return  A list of two items - a list of test statistics and a list of internal parameters
#' @export
vecnet <- function(x, y, covariate_y, max_iter=1000, tol=1e-10){

  # Initialise parameters
  param.list <- initialise_param.list(x=x,y=y,covariate_y=covariate_y)
  llik<-score<-tol_vec<-c()
  count=0

  # Convert parameters to the Expectation space
  param.list$S.list <- convert_param_to_S(param.list = param.list,alpha_run = T,beta_run = T)

  # Start timing
  a <- Sys.time()
  while (T){

    # Run internal EM algorithm
    param.list <- convert_S_to_param(param.list = param.list)

    # Update meta parameters
    llik <- c(llik,param.list$mllik)
    score <- c(score,param.list$mscore)
    tol_vec <- c(tol_vec,abs(tail(llik,2)[1]-tail(llik,1)[1]))

    # Check convergence
    if (count>25){
      if (count>max_iter | tail(tol_vec,1)<tol){
        break
      }
    }

    # Print updates
    if (count %% 10 == 0){
      b <- Sys.time()
      print(paste("Iterations ", count," in ", round(b-a,3)," with tol at",round(tail(tol_vec,1),20)))
    }

    count = count + 1
  }

  # Statistics for interpretation and analysis
  statistics <- extract_statistics(param.list = param.list)

  statistics$iter = count
  statistics$max_iter = max_iter
  statistics$tol = tol
  statistics$tol_vec = tol_vec
  statistics$llik = llik
  statistics$score = score

  return(statistics)
}

extract_statistics <- function(param.list){

  x <- as.matrix(param.list$x)
  y <- as.matrix(param.list$y)
  covariate_y <- as.matrix(param.list$covariate_y)
  y_hat_alpha_beta <- as.matrix(param.list$y_hat_alpha_beta)
  alpha <- as.matrix(param.list$alpha)
  beta <- as.matrix(param.list$beta)
  delta <- as.matrix(param.list$delta)
  intercept_alpha <- as.matrix(param.list$intercept_alpha)
  intercept_beta <- as.matrix(param.list$intercept_beta)
  intercept_x <- as.matrix(param.list$intercept_x)


  alpha_inverse <-t(alpha)%*%MASS::ginv(as.matrix(alpha%*%t(alpha)))
  beta_inverse <- MASS::ginv(as.matrix(t(beta)%*%beta))%*%t(beta)
  aligned_imputation <- as.matrix(
    x + alpha_inverse%*%intercept_x%*%beta_inverse + alpha_inverse%*%intercept_alpha + intercept_beta%*%beta_inverse
  )

  return(list(
    x = as.matrix((x)),
    y = as.matrix((y)),
    covariate_y <- as.matrix(covariate_y),
    aligned_imputation = as.matrix(aligned_imputation),
    fitted = as.matrix((y_hat_alpha_beta)),
    residuals = as.matrix((as.matrix((y)) - (y_hat_alpha_beta))),
    intercept_x = as.matrix(intercept_x),
    intercept_alpha = as.matrix(intercept_alpha),
    intercept_beta = as.matrix(intercept_beta),
    alpha = as.matrix(alpha),
    beta = as.matrix(beta),
    delta = as.matrix(delta)
  )
  )

}

initialise_param.list <- function(x,y,covariate_y){
  # Prepare y
  yT <- Matrix::t(y)
  y <- Matrix::Matrix(y)

  # # Prepare x
  xT <- Matrix::t(x)
  x <- Matrix::Matrix(x)

  if (is.vector(covariate_y)){
    n = 1
    p = length(covariate_y)
  }else{
    n = dim(covariate_y)[1]
    p = dim(covariate_y)[2]
  }

  covariate_y <- Matrix::Matrix(matrix(covariate_y,nrow = n,ncol=p))

  # Prepare matrix multiplication of x to save time later on
  internal_x <- xT%*%MASS::ginv(as.matrix(x%*%xT))
  internal_xT <- (MASS::ginv(as.matrix(xT%*%(x))))%*%xT

  # Initialise beta
  beta = internal_xT%*%matrix(1/(dim(x)[1]),nrow=dim(x)[1],ncol=dim(y)[1])%*%(y)

  # Initialise alpha
  alpha = (y)%*%Matrix::t(beta)%*%MASS::ginv(as.matrix(beta%*%Matrix::t(beta)))%*%internal_x

  delta = (y - alpha%*%x%*%beta)%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))

  # Initialise intercepts
  intercept_alpha <- Matrix::Matrix(1,nrow=dim(x)[1],ncol=dim(y)[2])
  intercept_beta <- Matrix::Matrix(1,nrow=dim(y)[1],ncol=dim(x)[2])
  intercept_x <- Matrix::Matrix(1,nrow=dim(x)[1],ncol=dim(x)[2])

  param.list <- list(x = x, y = y, covariate_y = covariate_y, alpha = alpha, beta = beta, delta = delta, intercept_x = intercept_x, intercept_alpha = intercept_alpha, intercept_beta = intercept_beta)

  return(param.list)

}

convert_param_to_S <- function(param.list,alpha_run=T,beta_run=T){

  # Prepare main parameters
  x <- param.list$x
  y <- param.list$y
  covariate_y <- param.list$covariate_y
  alpha <- param.list$alpha
  beta <- param.list$beta
  delta <- param.list$delta
  intercept_x <- param.list$intercept_x
  intercept_alpha <- param.list$intercept_alpha
  intercept_beta <- param.list$intercept_beta

  main_beta <- c()
  main_alpha <- c()

  # E-step of EM algorithm for beta parameter
  if (beta_run==T){

    y <- y - delta%*%covariate_y

    alphax <- alpha%*%x
    alphai <- alpha%*%intercept_x
    alphaa <- alpha%*%intercept_alpha

    t.alphax_y <- Matrix::t(alphax)%*%y
    t.alphai_y <- Matrix::t(alphai)%*%y
    t.alphax_alphaa <- Matrix::t(alphax)%*%alphaa

    t.alphax_alphax <- Matrix::t(alphax)%*%alphax
    t.alphax_alphai_p_alphax_t.alphai <- Matrix::t(alphax)%*%(alphai)+Matrix::t(Matrix::t(alphax)%*%(alphai))
    t.alphax_b_p_alphax_t.b <- Matrix::t(alphax)%*%(intercept_beta)+Matrix::t(Matrix::t(alphax)%*%(intercept_beta))

    main_alpha <- list(
      alphax = alphax,
      alphai = alphai,
      alphaa = alphaa,
      t.alphax_y = t.alphax_y,
      t.alphai_y = t.alphai_y,
      t.alphax_alphaa = t.alphax_alphaa,
      t.alphax_alphax = t.alphax_alphax,
      t.alphax_alphai_p_alphax_t.alphai = t.alphax_alphai_p_alphax_t.alphai,
      t.alphax_b_p_alphax_t.b = t.alphax_b_p_alphax_t.b
    )
  }

  # E-step of EM algorithm for beta parameter
  if (alpha_run==T){

    y <- y - delta%*%covariate_y

    xbeta <- x%*%beta
    ibeta <- intercept_x%*%beta
    bbeta <- intercept_beta%*%beta

    y_t.xbeta <- y%*%Matrix::t(xbeta)
    y_t.ibeta <- y%*%Matrix::t(ibeta)
    bbeta_t.xbeta <- bbeta%*%Matrix::t(xbeta)

    xbeta_t.xbeta <- xbeta%*%Matrix::t(xbeta)
    xbeta_t.ibeta_p_ibeta_t.xbeta <- xbeta%*%Matrix::t(ibeta)+Matrix::t(xbeta%*%Matrix::t(ibeta))
    xbeta_t.a_p_a_t.xbeta <- xbeta%*%Matrix::t(intercept_alpha)+Matrix::t(xbeta%*%Matrix::t(intercept_alpha))

    main_beta <- list(
      xbeta = xbeta,
      ibeta = ibeta,
      bbeta = bbeta,
      y_t.xbeta = y_t.xbeta,
      y_t.ibeta = y_t.ibeta,
      bbeta_t.xbeta = bbeta_t.xbeta,
      xbeta_t.xbeta = xbeta_t.xbeta,
      xbeta_t.ibeta_p_ibeta_t.xbeta = xbeta_t.ibeta_p_ibeta_t.xbeta,
      xbeta_t.a_p_a_t.xbeta = xbeta_t.a_p_a_t.xbeta
    )
  }

  return(c(main_alpha,main_beta))
}

convert_S_to_param <- function(param.list){

  # Prepare main parameters
  S.list <- param.list$S.list

  y <- param.list$y
  x <- param.list$x
  intercept_alpha <- param.list$intercept_alpha
  intercept_beta <- param.list$intercept_beta
  intercept_x <- param.list$intercept_x


  # Calculate parts of alpha parameter
  S.1_alpha <-
    S.list$y_t.xbeta +
    S.list$y_t.ibeta +
    - S.list$bbeta_t.xbeta +
    - (S.list$bbeta)%*%Matrix::t(S.list$ibeta) +
    - (S.list$bbeta)%*%Matrix::t(intercept_alpha) +
    y%*%Matrix::t(intercept_alpha)

  S.2_alpha <-
    S.list$xbeta_t.xbeta +
    S.list$xbeta_t.ibeta_p_ibeta_t.xbeta +
    S.list$xbeta_t.a_p_a_t.xbeta +
    (S.list$ibeta)%*%Matrix::t(S.list$ibeta) +
    intercept_alpha%*%Matrix::t(intercept_alpha) +
    (S.list$ibeta)%*%Matrix::t(intercept_alpha) +
    Matrix::t((S.list$ibeta)%*%Matrix::t(intercept_alpha))

  # Calculate alpha parameter
  alpha <- Matrix::Matrix(S.1_alpha%*%MASS::ginv(as.matrix(S.2_alpha)),sparse=T)
  param.list$alpha <- alpha
  # Update EM algorithm in E-step
  S.list <- convert_param_to_S(param.list = param.list,alpha_run = F,beta_run = T)

  # Calculate parts of beta parameter
  S.1_beta <-
    S.list$t.alphax_y +
    S.list$t.alphai_y +
    - S.list$t.alphax_alphaa +
    - Matrix::t(S.list$alphai)%*%(S.list$alphaa) +
    - Matrix::t(intercept_beta)%*%(S.list$alphaa) +
    Matrix::t(intercept_beta)%*%y


  S.2_beta <-
    S.list$t.alphax_alphax +
    S.list$t.alphax_alphai_p_alphax_t.alphai +
    S.list$t.alphax_b_p_alphax_t.b +
    Matrix::t(S.list$alphai)%*%(S.list$alphai) +
    Matrix::t(intercept_beta)%*%intercept_beta +
    Matrix::t(S.list$alphai)%*%intercept_beta +
    Matrix::t(Matrix::t(S.list$alphai)%*%intercept_beta)

  # Calculate beta parameter
  beta <- Matrix::Matrix(MASS::ginv(as.matrix(S.2_beta))%*%S.1_beta)
  param.list$beta <- beta
  # Update EM algorithm in E-step

  # Initialise intercept parameters
  intercept_x <- param.list$intercept_x
  intercept_alpha <- param.list$intercept_alpha
  intercept_beta <- param.list$intercept_beta

  # Calculate intercept parameters
  intercept_x <- y - alpha%*%intercept_alpha - intercept_beta%*%beta - alpha%*%x%*%beta
  intercept_alpha <- MASS::ginv(as.matrix(Matrix::t(alpha)%*%alpha))%*%Matrix::t(alpha)%*%(y - intercept_x - alpha%*%(x)%*%beta - intercept_beta%*%beta)
  intercept_beta <- (y - intercept_x - alpha%*%(x)%*%beta - alpha%*%intercept_alpha)%*%Matrix::t(beta)%*%MASS::ginv(as.matrix(beta%*%Matrix::t(beta)))

  param.list$intercept_x <- intercept_x
  param.list$intercept_alpha <- intercept_alpha
  param.list$intercept_beta <- intercept_beta

  delta = (y - alpha%*%x%*%beta - intercept_x - alpha%*%intercept_alpha - intercept_beta%*%beta )%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))

  S.list <- convert_param_to_S(param.list,alpha_run = T,beta_run = F)

  # Calculate fitted value of y; covariance of y ; likelihood ; mean-squared error
  y_hat_alpha_beta <- alpha%*%(as.matrix(x))%*%beta+alpha%*%intercept_alpha+intercept_beta%*%beta +intercept_x + delta%*%covariate_y
  D1 = Matrix::diag(Matrix::diag(Matrix::t((y) - (y_hat_alpha_beta))%*%((y) - (y_hat_alpha_beta))/dim(y)[1]))
  mllik <- c(-sum(log(diag(D1)))-sum(Matrix::diag(((y) - (y_hat_alpha_beta))%*%Matrix::diag(1/Matrix::diag(D1))%*%Matrix::t((y) - (y_hat_alpha_beta)))))
  mscore <- mean(sqrt(as.matrix(y - y_hat_alpha_beta)^2))
  return(list(
    S.list = S.list,
    x = x,
    y = y,
    covariate_y = covariate_y,
    intercept_x = intercept_x,
    intercept_alpha = intercept_alpha,
    intercept_beta = intercept_beta,
    y_hat_alpha_beta = y_hat_alpha_beta,
    alpha = alpha,
    beta = beta,
    delta = delta,
    intercept_x = intercept_x,
    intercept_alpha = intercept_alpha,
    intercept_beta = intercept_beta,
    mllik = mllik,
    mscore = mscore
  ))
}
