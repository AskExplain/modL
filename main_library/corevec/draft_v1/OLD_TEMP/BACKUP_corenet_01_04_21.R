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
vecnet <- function(x, y, covariate_x, covariate_y, k_dim = 10, max_iter=1000, min_iter = 35, tol=1e-10, log=T, factor=T){
  
  cx <- Matrix::colSums(x)
  cy <- Matrix::colSums(y)
  
  if (log == T){
    x <- log(1+x)
    y <- log(1+y)
  }
  if (factor==T){
    x <- x - Matrix::t(cx%*%MASS::ginv(as.matrix(Matrix::t(cx)%*%cx))%*%Matrix::t(cx)%*%Matrix::t(x))
    y <- y - Matrix::t(cy%*%MASS::ginv(as.matrix(Matrix::t(cy)%*%cy))%*%Matrix::t(cy)%*%Matrix::t(y))
  }
  
  # Initialise parameters
  param.list <- initialise_param.list(x=x,y=y,covariate_y=covariate_y, covariate_x = covariate_x, alpha.init = NULL, beta.init =  NULL, k_dim = k_dim)
  
  llik<-score<-tol_vec<-beta_vec<-alpha.L_vec<-alpha.K_vec<-c()
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
    beta_vec <- c(beta_vec,param.list$beta[1,2])
    alpha.L_vec <- c(alpha.L_vec,param.list$alpha.L[1,2])
    alpha.K_vec <- c(alpha.K_vec,param.list$alpha.K[1,2])
    
    par(mfcol=c(5,1))
    plot(llik,main=tail(llik,1))
    plot(score,main=tail(score,1))
    plot(beta_vec,main=tail(beta_vec,1))
    plot(alpha.L_vec,main=tail(alpha.L_vec,1))
    plot(alpha.K_vec,main=tail(alpha.K_vec,1))
    
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
  statistics$tol_vec = tol_vec
  statistics$llik = llik
  statistics$score = score
  statistics$beta_vec = beta_vec
  statistics$alpha.L_vec = alpha.L_vec
  statistics$alpha.K_vec = alpha.K_vec
  
  return(list(statistics = statistics,
              internal_parameters = param.list))
}



fit.vecnet <- function(vecnet_model,z,covariate_z_remove=NULL){
  
  z <- apply(z,2,as.numeric)
  
  cell_types_z <- colnames(z)
  cell_types_x <- colnames(vecnet_model$statistics$x)
  cell_types_y <- colnames(vecnet_model$statistics$y)
  
  covariate_x.cell_type <- do.call('cbind',lapply(sort(unique(cell_types_x)),function(cell_id){
    c(cell_types_x == cell_id)+0
  }))
  
  covariate_z.cell_type <- do.call('cbind',lapply(sort(unique(cell_types_z)),function(cell_id){
    c(cell_types_z == cell_id)+0
  }))
  
  unique_cell_types_z <- unique(cell_types_z)
  alpha <- as.matrix(vecnet_model$statistics$alpha)
  beta.x <- as.matrix(vecnet_model$statistics$beta)
  
  fitted <- as.matrix(vecnet_model$statistics$fitted)
  y <- as.matrix(vecnet_model$statistics$y)
  x <- as.matrix(vecnet_model$statistics$x)
  Z_obs <- z
  
  if (!is.null(covariate_z_remove)){
    library_size <- matrix(colSums(z),nrow=dim(z)[2],ncol=1)
    covar_z <- cbind(covariate_z_remove,library_size)
    
    beta.Z_obs <- ((MASS::ginv(as.matrix(t(covar_z)%*%(covar_z)))%*%t(covar_z)%*%t(Z_obs)))
    Z_obs <- Z_obs - t(covar_z %*% beta.Z_obs)
  }
  
  fit.z <- alpha%*%Z_obs
  beta.Z_obs <- ((MASS::ginv(as.matrix(t(fit.z)%*%(fit.z)))%*%t(fit.z)%*%(y)))
  beta.inv <- beta.x%*%t(beta.Z_obs)%*%MASS::ginv(as.matrix((beta.Z_obs)%*%t(beta.Z_obs)))
  
  init_params <- vecnet::initialise_param.list(x = x, y = fit.z, covariate_y = t(covariate_z.cell_type), covariate_x = t(covariate_x.cell_type), alpha.init =  alpha, beta.init = beta.inv)
  init_params$S.list <- vecnet::convert_param_to_S(param.list = init_params, alpha_run = F, beta_run = F)
  init_params <- vecnet::convert_S_to_param(param.list = init_params, alpha_run = F, beta_run = F)
  
  Z_align <- as.matrix(init_params$alpha%*%init_params$x%*%init_params$beta + init_params$delta%*%init_params$covariate_y + init_params$gamma%*%init_params$covariate_x%*%init_params$beta)
  row.names(Z_align) <- row.names(z)
  colnames(Z_align) <- colnames(z)
  
  return(list(Z_adjust = Z_obs, Z_align = Z_align, beta.Z_align = beta.Z_obs))
}



extract_statistics <- function(param.list){
  
  x <- as.matrix(param.list$x)
  y <- as.matrix(param.list$y)
  covariate_y <- as.matrix(param.list$covariate_y)
  covariate_x <- as.matrix(param.list$covariate_x)
  y_hat_alpha_beta <- as.matrix(param.list$y_hat_alpha_beta)
  alpha.L <- as.matrix(param.list$alpha.L)
  alpha.K <- as.matrix(param.list$alpha.K)
  beta <- as.matrix(param.list$beta)
  delta <- as.matrix(param.list$delta)
  gamma = as.matrix(param.list$gamma)
  
  return(list(
    x = as.matrix((x)),
    y = as.matrix((y)),
    covariate_x = as.matrix(covariate_x),
    covariate_y = as.matrix(covariate_y),
    fitted = as.matrix((y_hat_alpha_beta)),
    residuals = as.matrix((as.matrix((y)) - (y_hat_alpha_beta))),
    alpha.L = as.matrix(alpha.L),
    alpha.K = as.matrix(alpha.K),
    beta = as.matrix(beta),
    delta = as.matrix(delta),
    gamma = as.matrix(gamma)
  )
  )
}

initialise_param.list <- function(x,y,covariate_y,covariate_x,alpha.init=NULL,beta.init=NULL, k_dim = 10){
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
  
  alpha.L <- Matrix::t(Matrix::Matrix(RSpectra::eigs(as.matrix(alpha%*%Matrix::t(alpha)),k=k_dim)$vectors))
  # alpha.K <- (alpha.L%*%x%*%beta)%*%Matrix::t(y)%*%MASS::ginv(as.matrix(y%*%Matrix::t(y)))
  alpha.K <- alpha.L
  
  delta = (alpha.K%*%y - alpha.L%*%x%*%beta)%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))
  
  gamma = (alpha.K%*%y - alpha.L%*%x%*%beta - delta%*%covariate_y)%*%Matrix::t(covariate_x%*%beta)%*%MASS::ginv(as.matrix((covariate_x%*%beta)%*%Matrix::t(covariate_x%*%beta)))
  
  # Initialise intercepts
  intercept_x <- Matrix::Matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  intercept_alpha <- Matrix::Matrix(0,nrow=dim(x)[1],ncol=dim(y)[2])
  intercept_beta <- Matrix::Matrix(0,nrow=dim(y)[1],ncol=dim(x)[2])
  
  param.list <- list(x = x,
                     y = y,
                     alpha.L = alpha.L,
                     alpha.K = alpha.K,
                     beta = beta,
                     delta = delta,
                     gamma = gamma,
                     covariate_x = covariate_x,
                     covariate_y = covariate_y,
                     intercept_x = intercept_x,
                     intercept_alpha = intercept_alpha,
                     intercept_beta = intercept_beta)
  
  return(param.list)
  
}

convert_param_to_S <- function(param.list,alpha_run=T,beta_run=T){
  
  # Prepare main parameters
  x <- param.list$x
  y <- param.list$y
  covariate_y <- param.list$covariate_y
  alpha.L <- param.list$alpha.L
  alpha.K <- param.list$alpha.K
  beta <- param.list$beta
  delta <- param.list$delta
  intercept_x <- param.list$intercept_x
  intercept_alpha <- param.list$intercept_alpha
  intercept_beta <- param.list$intercept_beta
  
  main_beta <- c()
  main_alpha <- c()
  
  # E-step of EM algorithm for beta parameter
  if (beta_run==T){
    
    alphax <- alpha.L%*%x
    alphai <- alpha.L%*%intercept_x
    
    t.alphax_y <- Matrix::t(alphax)%*%alpha.K%*%y
    t.alphai_y <- Matrix::t(alphai)%*%alpha.K%*%y
    
    t.alphax_alphax <- Matrix::t(alphax)%*%alphax
    
    t.alphax_alphai_p_alphax_t.alphai <- Matrix::t(alphax)%*%(alphai)+Matrix::t(Matrix::t(alphax)%*%(alphai))
    
    main_alpha <- list(
      alphax = alphax,
      alphai = alphai,
      t.alphax_y = t.alphax_y,
      t.alphai_y = t.alphai_y,
      t.alphax_alphax = t.alphax_alphax,
      t.alphax_alphai_p_alphax_t.alphai = t.alphax_alphai_p_alphax_t.alphai
    )
  }
  
  # E-step of EM algorithm for beta parameter
  if (alpha_run==T){
    
    xbeta <- x%*%beta
    ibeta <- intercept_x%*%beta
    
    y_t.xbeta <- alpha.K%*%y%*%Matrix::t(xbeta)
    y_t.ibeta <- alpha.K%*%y%*%Matrix::t(ibeta)
    
    xbeta_t.xbeta <- xbeta%*%Matrix::t(xbeta)
    xbeta_t.ibeta_p_ibeta_t.xbeta <- xbeta%*%Matrix::t(ibeta)+Matrix::t(xbeta%*%Matrix::t(ibeta))
    
    main_beta <- list(
      xbeta = xbeta,
      ibeta = ibeta,
      y_t.xbeta = y_t.xbeta,
      y_t.ibeta = y_t.ibeta,
      xbeta_t.xbeta = xbeta_t.xbeta,
      xbeta_t.ibeta_p_ibeta_t.xbeta = xbeta_t.ibeta_p_ibeta_t.xbeta
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
  alpha.K <- param.list$alpha.K
  intercept_alpha <- param.list$intercept_alpha
  intercept_beta <- param.list$intercept_beta
  intercept_x <- param.list$intercept_x
  covariate_x <- param.list$covariate_x
  covariate_y <- param.list$covariate_y
  delta <- param.list$delta
  gamma <- param.list$gamma
  
  if (alpha_run){
    # Calculate parts of alpha parameter
    S.1_alpha <-
      S.list$y_t.xbeta +
      S.list$y_t.ibeta -
      delta%*%covariate_y%*%Matrix::t(S.list$xbeta) -
      delta%*%covariate_y%*%Matrix::t(S.list$ibeta) -
      gamma%*%covariate_x%*%beta%*%Matrix::t(S.list$xbeta) -
      gamma%*%covariate_x%*%beta%*%Matrix::t(S.list$ibeta)
    
    
    S.2_alpha <-
      S.list$xbeta_t.xbeta +
      S.list$xbeta_t.ibeta_p_ibeta_t.xbeta +
      (S.list$ibeta)%*%Matrix::t(S.list$ibeta)
    
    # Calculate alpha parameter
    alpha.L <- Matrix::Matrix(S.1_alpha%*%MASS::ginv(as.matrix(S.2_alpha)),sparse=T)
    alpha.K <- param.list$alpha.K <- param.list$alpha.L <- alpha.L
    # Update EM algorithm in E-step
    S.list <- convert_param_to_S(param.list = param.list,alpha_run = alpha_run, beta_run = beta_run)
  }
  
  # alpha.K <- (alpha.L%*%x%*%beta)%*%Matrix::t(y)%*%MASS::ginv(as.matrix(y%*%Matrix::t(y)))
  
  if (beta_run){
    # Calculate parts of beta parameter
    S.1_beta <-
      S.list$t.alphax_y +
      S.list$t.alphai_y +
      Matrix::t(gamma%*%covariate_x)%*%alpha.K%*%y -
      Matrix::t(S.list$alphax)%*%delta%*%covariate_y -
      Matrix::t(S.list$alphai)%*%delta%*%covariate_y -
      Matrix::t(gamma%*%covariate_x)%*%delta%*%covariate_y
    
    
    S.2_beta <-
      S.list$t.alphax_alphax +
      S.list$t.alphax_alphai_p_alphax_t.alphai +
      Matrix::t(S.list$alphai)%*%(S.list$alphai) +
      Matrix::t(gamma%*%covariate_x)%*%S.list$alphax +
      Matrix::t(gamma%*%covariate_x)%*%S.list$alphai +
      Matrix::t(S.list$alphax)%*%(gamma%*%covariate_x) +
      Matrix::t(S.list$alphai)%*%(gamma%*%covariate_x) +
      Matrix::t(gamma%*%covariate_x)%*%(gamma%*%covariate_x)
    
    # Calculate beta parameter
    beta <- Matrix::Matrix(MASS::ginv(as.matrix(S.2_beta))%*%S.1_beta)
    param.list$beta <- beta
  }
  
  
  delta = (alpha.K%*%y - alpha.L%*%(x)%*%beta - gamma%*%covariate_x%*%beta )%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))
  
  gamma = (alpha.K%*%y - alpha.L%*%(x)%*%beta - delta%*%covariate_y)%*%Matrix::t(covariate_x%*%beta)%*%MASS::ginv(as.matrix((covariate_x%*%beta)%*%Matrix::t(covariate_x%*%beta)))
  
  param.list$delta <- delta
  param.list$gamma <- gamma
  
  S.list <- convert_param_to_S(param.list,alpha_run = alpha_run,beta_run = beta_run)
  
  # Calculate fitted value of y; covariance of y ; likelihood ; mean-squared error
  y_hat_alpha_beta <- as.matrix(x)%*%beta + MASS::ginv(as.matrix(Matrix::t(alpha.L)%*%alpha.L))%*%Matrix::t(alpha.L)%*%(delta%*%covariate_y + gamma%*%covariate_x%*%beta)
  D1 = Matrix::diag(Matrix::diag(Matrix::t((y) - (y_hat_alpha_beta))%*%((y) - (y_hat_alpha_beta))/dim(y)[1]))
  mllik <- sum(mclust::dmvnorm(y - y_hat_alpha_beta,rep(0,dim(D1)[1]),D1,log = T))
  mscore <- mean(sqrt(as.matrix(y - y_hat_alpha_beta)^2))
  return(list(
    S.list = S.list,
    x = x,
    y = y,
    covariate_x = covariate_x,
    covariate_y = covariate_y,
    y_hat_alpha_beta = y_hat_alpha_beta,
    alpha.L = alpha.L,
    alpha.K = alpha.K,
    beta = beta,
    delta = delta,
    gamma = gamma,
    intercept_x = intercept_x,
    intercept_alpha = intercept_alpha,
    intercept_beta = intercept_beta,
    mllik = mllik,
    mscore = mscore
  ))
}
