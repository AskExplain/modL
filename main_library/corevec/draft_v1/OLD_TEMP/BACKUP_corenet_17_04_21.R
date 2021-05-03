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
  param.list$B.list <- param.list
  
  # Convert parameters to the Expectation space
  param.list$S.list <- convert_param_to_S(param.list = param.list,alpha_run = T,beta_run = T)
  param.list <- convert_S_to_param(param.list = param.list, alpha_run = T, beta_run = T, subset_run = F)
  
  llik<-score<-tol_vec<-beta_vec<-alpha.L_vec<-c()
  count=0
  
  # Start timing
  a <- Sys.time()
  while (T){
    
    # Run internal EM algorithm
    param.list <- convert_S_to_param(param.list = param.list, alpha_run = T, beta_run = T, subset_run = T)
    
    # Update meta parameters
    llik <- c(llik,param.list$mllik)
    score <- c(score,param.list$mscore)
    tol_vec <- c(tol_vec,abs(tail(llik,2)[1]-tail(llik,1)[1]))
    beta_vec <- c(beta_vec,param.list$beta[1,2])
    alpha.L_vec <- c(alpha.L_vec,param.list$alpha.L[1,2])
    
    par(mfcol=c(6,1))
    plot(llik,main=tail(llik,1))
    plot(score,main=tail(score,1))
    plot(beta_vec,main=tail(beta_vec,1))
    plot(alpha.L_vec,main=tail(alpha.L_vec,1))
    barplot(c(length(param.list$gene.prune),length(param.list$cell_x.prune),length(param.list$cell_y.prune)))
    plot(as.matrix(Matrix::t(param.list$alpha.L%*%param.list$x)),col=as.factor(colnames(param.list$x)))
    
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
  
  return(list(statistics = statistics,
              internal_parameters = param.list))
}



extract_statistics <- function(param.list){
  
  x <- as.matrix(param.list$x)
  y <- as.matrix(param.list$y)
  covariate_y <- as.matrix(param.list$covariate_y)
  covariate_x <- as.matrix(param.list$covariate_x)
  y_hat_alpha_beta <- as.matrix(param.list$y_hat_alpha_beta)
  alpha.L <- as.matrix(param.list$alpha.L)
  beta <- as.matrix(param.list$beta)
  delta <- as.matrix(param.list$delta)
  
  return(list(
    x = as.matrix((x)),
    y = as.matrix((y)),
    covariate_y = as.matrix(covariate_y),
    fitted = as.matrix((y_hat_alpha_beta)),
    residuals = as.matrix((as.matrix((y)) - (y_hat_alpha_beta))),
    alpha.L = as.matrix(alpha.L),
    beta = as.matrix(beta),
    delta = as.matrix(delta)
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
  
  
  
  if (is.null(beta.init)){
    beta = MASS::ginv(as.matrix(Matrix::t(x)%*%(x)))%*%Matrix::t((x))%*%(y)
  } else {
    beta = beta.init
  }
  
  if (is.null(alpha.init)){
    alpha = (y)%*%Matrix::t(beta)%*%MASS::ginv(as.matrix(beta%*%Matrix::t(beta)))%*%xT
  } else {
    alpha = alpha.init
  }
  
  alpha.L <- Matrix::t(Matrix::Matrix(RSpectra::eigs(as.matrix(alpha%*%Matrix::t(alpha)),k=k_dim)$vectors))
  
  delta = (alpha.L%*%y - alpha.L%*%x%*%beta)%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))
  
  param.list <- list(x = x,
                     y = y,
                     alpha.L = alpha.L,
                     beta = beta,
                     delta = delta,
                     covariate_y = covariate_y)
  
  return(param.list)
  
}

convert_param_to_S <- function(param.list,alpha_run=T,beta_run=T){
  
  # Prepare main parameters
  x <- param.list$B.list$x
  y <- param.list$B.list$y
  covariate_y <- param.list$B.list$covariate_y
  alpha.L <- param.list$B.list$alpha.L
  beta <- param.list$B.list$beta
  delta <- param.list$B.list$delta
  
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
  
  # E-step of EM algorithm for beta parameter
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



subset_param.list <- function(param.list){
  
  gene.prune = param.list$gene.prune
  if (length(gene.prune)==0){
    gene.prune = c(param.list$cut_gene.prune)
  } else {
    # Update top 500 genes
    param.list$cut_gene.prune[param.list$cut_gene.prune%in%gene.prune] <- setdiff(param.list$all_gene.prune,param.list$cut_gene.prune)[1:length(gene.prune)]
  }
  cell_x.prune = param.list$cell_x.prune
  if (length(cell_x.prune)==0){
    cell_x.prune = c(param.list$cut_cell_x.prune)
  } else {
    # Update top 100 x.cells
    param.list$cut_cell_x.prune[param.list$cut_cell_x.prune%in%cell_x.prune] <- setdiff(param.list$all_cell_x.prune,param.list$cut_cell_x.prune)[1:length(cell_x.prune)]
    
  }
  cell_y.prune = param.list$cell_y.prune
  if (length(cell_y.prune)==0){
    cell_y.prune = c(param.list$cut_cell_y.prune)
  } else {
    # Update top 100 y.cells
    param.list$cut_cell_y.prune[param.list$cut_cell_y.prune%in%cell_y.prune] <- setdiff(param.list$all_cell_y.prune,param.list$cut_cell_y.prune)[1:length(cell_y.prune)]
  }
  
  gene.prune <- param.list$cut_gene.prune
  cell_x.prune <- param.list$cut_cell_x.prune
  cell_y.prune <- param.list$cut_cell_y.prune
  
  param.list$x <- param.list$x[gene.prune,cell_x.prune]
  param.list$y <- param.list$y[gene.prune,cell_y.prune]
  param.list$covariate_y <- param.list$covariate_y[,cell_y.prune]
  param.list$alpha.L <- param.list$alpha.L[,gene.prune]
  param.list$beta <- param.list$beta[cell_x.prune,cell_y.prune]
  param.list$beta.K <- param.list$beta.K[cell_y.prune,]
  param.list$delta <- param.list$delta[gene.prune,]
  
  return(param.list)
  
}



update_param.list <- function(param.list){
  
  gene.prune = param.list$gene.prune
  if (length(gene.prune)==0){
    gene.prune = c(param.list$cut_gene.prune)
  } else {
    # Update top 500 genes
    param.list$cut_gene.prune[param.list$cut_gene.prune%in%gene.prune] <- setdiff(param.list$all_gene.prune,param.list$cut_gene.prune)[1:length(gene.prune)]
  }
  cell_x.prune = param.list$cell_x.prune
  if (length(cell_x.prune)==0){
    cell_x.prune = c(param.list$cut_cell_x.prune)
  } else {
    # Update top 100 x.cells
    param.list$cut_cell_x.prune[param.list$cut_cell_x.prune%in%cell_x.prune] <- setdiff(param.list$all_cell_x.prune,param.list$cut_cell_x.prune)[1:length(cell_x.prune)]
    
  }
  cell_y.prune = param.list$cell_y.prune
  if (length(cell_y.prune)==0){
    cell_y.prune = c(param.list$cut_cell_y.prune)
  } else {
    # Update top 100 y.cells
    param.list$cut_cell_y.prune[param.list$cut_cell_y.prune%in%cell_y.prune] <- setdiff(param.list$all_cell_y.prune,param.list$cut_cell_y.prune)[1:length(cell_y.prune)]
  }
  
  
  gene.prune <- param.list$cut_gene.prune
  cell_x.prune <- param.list$cut_cell_x.prune
  cell_y.prune <- param.list$cut_cell_y.prune
  
  
  param.list$alpha.L <- Matrix::Matrix(as.matrix(param.list$alpha.L))
  
  # Update parameters
  param.list$alpha.L[,gene.prune] <- Matrix::Matrix(as.matrix(param.list$B.list$alpha.L))
  param.list$beta[cell_x.prune,cell_y.prune] <- param.list$B.list$beta
  param.list$delta[gene.prune,] <- param.list$B.list$delta
  param.list$beta.K[cell_y.prune,] <- param.list$B.list$beta.K
  
  
  return(param.list)
  
}




convert_S_to_param <- function(param.list, alpha_run=T, beta_run=T, subset_run = F){
  
  # Prepare main parameters
  S.list <- param.list$S.list
  
  y <- param.list$y
  x <- param.list$x
  beta <- param.list$beta
  alpha.L <- param.list$alpha.L
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
      
      S.2_alpha <-
        S.list$xbeta_t.xbeta +
        S.list$xbeta_t.ibeta_p_ibeta_t.xbeta +
        (S.list$ibeta)%*%Matrix::t(S.list$ibeta)
      
      # Calculate alpha parameter
      alpha.L <- Matrix::Matrix(S.1_alpha%*%MASS::ginv(as.matrix(S.2_alpha)),sparse=T)
      param.list$alpha.L <- alpha.L
      
      # Update EM algorithm in E-step
      S.list <- convert_param_to_S(param.list = param.list,alpha_run = alpha_run, beta_run = beta_run)
  }
  
  if (beta_run){
    # Calculate parts of beta parameter
    S.1_beta <-
      S.list$t.alphax_y +
      S.list$t.alphai_y +
      Matrix::t(gamma%*%covariate_x)%*%alpha.L%*%y -
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
  
  
  delta = (alpha.L%*%y - alpha.L%*%(x)%*%beta - gamma%*%covariate_x%*%beta )%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))
  
  param.list$delta <- delta
  S.list <- convert_param_to_S(param.list,alpha_run = alpha_run,beta_run = beta_run)
  
  
  # Calculate fitted value of y; covariance of y ; likelihood ; mean-squared error
  param.list$y_hat_alpha_beta <- as.matrix(x)%*%beta + (delta%*%(covariate_y))
  D1 = Matrix::diag(Matrix::diag(((y) - (param.list$y_hat_alpha_beta))%*%Matrix::t((y) - (param.list$y_hat_alpha_beta))/dim(y)[1]))
  D2 = Matrix::diag(Matrix::diag(Matrix::t((y) - (param.list$y_hat_alpha_beta))%*%((y) - (param.list$y_hat_alpha_beta))/dim(y)[1]))
  
  param.list$mllik <- sum(mclust::dmvnorm((y - param.list$y_hat_alpha_beta),rep(0,dim(D2)[1]),D2,log = T))
  param.list$mscore <- mean(sqrt(as.matrix(y - param.list$y_hat_alpha_beta)^2))
  
  if (subset_run == T){
    param.list <- update_param.list(param.list)
  }
  
  cov_alpha.L = Matrix::diag(Matrix::t(alpha.L)%*%alpha.L)
  cov_beta.K = Matrix::diag(MASS::ginv(as.matrix(Matrix::t(x%*%beta+delta%*%covariate_y)%*%Matrix::diag(1/Matrix::diag(1e-300+D1))%*%(x%*%beta+delta%*%covariate_y))))
  cov_beta.beta.K = Matrix::diag(MASS::ginv(as.matrix(Matrix::t(x+delta%*%covariate_y%*%MASS::ginv(as.matrix(Matrix::t(beta)%*%beta))%*%Matrix::t(beta))%*%Matrix::diag(1/Matrix::diag(1e-300+D1))%*%(x+delta%*%covariate_y%*%MASS::ginv(as.matrix(Matrix::t(beta)%*%beta))%*%Matrix::t(beta)))))
  
  alpha.L_check.restart = pt((as.matrix(alpha.L)/sqrt(cov_alpha.L+1e-300)),length(cov_alpha.L)-1,lower.tail = T)
  beta.K_check.restart = pt(as.matrix((beta.K)/sqrt(cov_beta.K+1e-300)),length(cov_beta.K)-1,lower.tail = T)
  beta.beta.K_check.restart = pt(as.matrix(beta%*%beta.K)/sqrt(cov_beta.beta.K+1e-300),length(cov_beta.beta.K)-1,lower.tail = T)
  
  param.list$gene.prune <- unique(c(gene.prune,param.list$cut_gene.prune[unique(c(which(colSums(alpha.L_check.restart<0.05/prod(dim(alpha.L_check.restart)))==k_dim)))]))
  param.list$cell_y.prune <- unique(c(cell_y.prune,param.list$cut_cell_y.prune[unique(c(which(rowSums(beta.K_check.restart<0.05/prod(dim(beta_check.restart)))==1)))]))
  param.list$cell_x.prune <- unique(c(cell_x.prune,param.list$cut_cell_x.prune[unique(c(which(rowSums(beta.beta.K_check.restart<0.05/prod(dim(beta.beta.K_check.restart)))==1)))]))
  
  param.list$B.list$S.list <- S.list
  
  print(length(param.list$gene.prune))
  print(sort(row.names(param.list$x)[param.list$gene.prune]))
  
  
  
  
  return(list(
    S.list = S.list,
    x = x,
    y = y,
    covariate_x = covariate_x,
    covariate_y = covariate_y,
    y_hat_alpha_beta = y_hat_alpha_beta,
    alpha.L = alpha.L,
    alpha.L = alpha.L,
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
