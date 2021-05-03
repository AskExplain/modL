#' #' Gaussian Vector Network Regression
#' #'
#' #' Gaussian Vector Network Regression is a method that uses parameters to reweight both samples and features.
#' #'
#' #' @param x Reference dataset of cell by gene matrix
#' #' @param y Experimental dataset of cell by gene matrix
#' #' @param max_iter Maximum number of iterations possible
#' #' @param tol Tolerance value to assess convergence (converges when less than tol)
#' #' @return  A list of two items - a list of test statistics and a list of internal parameters
#' #' @export
#' vecnet <- function(x, y, covariate_y, k_dim = 10, max_iter=1000, min_iter = 35, tol=1e-10){
#'   
#'   # cx <- Matrix::colSums(x)
#'   # cy <- Matrix::colSums(y)
#'   #
#'   # if (log.transform == T){
#'   #   x <- log(1+x)
#'   #   y <- log(1+y)
#'   # }
#'   # if (factor==T){
#'   #   x <- x - Matrix::t(cx%*%MASS::ginv(as.matrix(Matrix::t(cx)%*%cx))%*%Matrix::t(cx)%*%Matrix::t(x))
#'   #   y <- y - Matrix::t(cy%*%MASS::ginv(as.matrix(Matrix::t(cy)%*%cy))%*%Matrix::t(cy)%*%Matrix::t(y))
#'   # }
#'   
#'   # Initialise parameters
#'   param.list <- initialise_param.list(x = x, y = y, covariate_y = covariate_y, alpha.init = NULL, beta.init = NULL, k_dim  = 10)
#'   param.list <- calculate_scores_and_update_params(param.list = param.list, alpha.L_run = T, beta_run = T, delta_run = T, D_run = T)
#'   
#'   llik<-score<-tol_vec<-beta_vec<-alpha.L_vec<-delta_vec<-D_vec<-c()
#'   count=0
#'   
#'   # Start timing
#'   a <- Sys.time()
#'   while (T){
#'     
#'     # Run internal EM algorithm
#'     param.list <- calculate_scores_and_update_params(param.list = param.list, alpha.L_run = T, beta_run = T, delta_run = T, D_run = T)
#'     
#'     # Update meta parameters
#'     llik <- c(llik,param.list$mllik)
#'     score <- c(score,param.list$mscore)
#'     tol_vec <- c(tol_vec,abs(tail(llik,2)[1]-tail(llik,1)[1]))
#'     beta_vec <- c(beta_vec,param.list$beta[1,2])
#'     alpha.L_vec <- c(alpha.L_vec,param.list$alpha.L[1,2])
#'     delta_vec <- c(delta_vec,param.list$delta[1,2])
#'     D_vec <- c(D_vec,param.list$D[1,1])
#'     
#'     par(mfcol=c(3,2))
#'     plot(llik,main=tail(llik,1))
#'     plot(score,main=tail(score,1))
#'     plot(beta_vec,main=tail(beta_vec,1))
#'     plot(alpha.L_vec,main=tail(alpha.L_vec,1))
#'     plot(delta_vec,main=tail(delta_vec,1))
#'     plot(D_vec,main=tail(D_vec,1))
#'     
#'     # Check convergence
#'     if (count>min_iter){
#'       if (tail(llik,2)[2]<tail(llik,2)[1]){
#'         break
#'       }
#'       if (count>max_iter | tail(tol_vec,1)<tol){
#'         break
#'       }
#'     }
#'     
#'     # Print updates
#'     if (count %% 10 == 0){
#'       b <- Sys.time()
#'       print(paste("Iterations ", count," in ", round(b-a,3)," with tol at",round(tail(tol_vec,1),20)))
#'     }
#'     
#'     count = count + 1
#'   }
#'   
#'   # Statistics for interpretation and analysis
#'   statistics <- extract_statistics(param.list = param.list)
#'   
#'   statistics$iter = count
#'   statistics$max_iter = max_iter
#'   statistics$tol_vec = tol_vec
#'   statistics$llik = llik
#'   statistics$score = score
#'   statistics$beta_vec = beta_vec
#'   statistics$alpha.L_vec = alpha.L_vec
#'   statistics$alpha.K_vec = alpha.K_vec
#'   
#'   return(list(statistics = statistics,
#'               internal_parameters = param.list))
#' }
#' 
#' 
#' 
#' fit.vecnet <- function(vecnet_model,z,covariate_z_remove=NULL){
#'   
#'   z <- apply(z,2,as.numeric)
#'   
#'   cell_types_z <- colnames(z)
#'   cell_types_x <- colnames(vecnet_model$statistics$x)
#'   cell_types_y <- colnames(vecnet_model$statistics$y)
#'   
#'   covariate_x.cell_type <- do.call('cbind',lapply(sort(unique(cell_types_x)),function(cell_id){
#'     c(cell_types_x == cell_id)+0
#'   }))
#'   
#'   covariate_z.cell_type <- do.call('cbind',lapply(sort(unique(cell_types_z)),function(cell_id){
#'     c(cell_types_z == cell_id)+0
#'   }))
#'   
#'   unique_cell_types_z <- unique(cell_types_z)
#'   alpha <- as.matrix(vecnet_model$statistics$alpha)
#'   beta.x <- as.matrix(vecnet_model$statistics$beta)
#'   
#'   fitted <- as.matrix(vecnet_model$statistics$fitted)
#'   y <- as.matrix(vecnet_model$statistics$y)
#'   x <- as.matrix(vecnet_model$statistics$x)
#'   Z_obs <- z
#'   
#'   if (!is.null(covariate_z_remove)){
#'     library_size <- matrix(colSums(z),nrow=dim(z)[2],ncol=1)
#'     covar_z <- cbind(covariate_z_remove,library_size)
#'     
#'     beta.Z_obs <- ((MASS::ginv(as.matrix(t(covar_z)%*%(covar_z)))%*%t(covar_z)%*%t(Z_obs)))
#'     Z_obs <- Z_obs - t(covar_z %*% beta.Z_obs)
#'   }
#'   
#'   fit.z <- alpha%*%Z_obs
#'   beta.Z_obs <- ((MASS::ginv(as.matrix(t(fit.z)%*%(fit.z)))%*%t(fit.z)%*%(y)))
#'   beta.inv <- beta.x%*%t(beta.Z_obs)%*%MASS::ginv(as.matrix((beta.Z_obs)%*%t(beta.Z_obs)))
#'   
#'   init_params <- vecnet::initialise_param.list(x = x, y = fit.z, covariate_y = t(covariate_z.cell_type), covariate_x = t(covariate_x.cell_type), alpha.init =  alpha, beta.init = beta.inv)
#'   init_params$S.list <- vecnet::convert_param_to_S(param.list = init_params, alpha_run = F, beta_run = F)
#'   init_params <- vecnet::convert_S_to_param(param.list = init_params, alpha_run = F, beta_run = F)
#'   
#'   Z_align <- as.matrix(init_params$alpha%*%init_params$x%*%init_params$beta + init_params$delta%*%init_params$covariate_y + init_params$gamma%*%init_params$covariate_x%*%init_params$beta)
#'   row.names(Z_align) <- row.names(z)
#'   colnames(Z_align) <- colnames(z)
#'   
#'   return(list(Z_adjust = Z_obs, Z_align = Z_align, beta.Z_align = beta.Z_obs))
#' }
#' 
#' 
#' 
#' extract_statistics <- function(param.list){
#'   
#'   x <- as.matrix(param.list$x)
#'   y <- as.matrix(param.list$y)
#'   covariate_y <- as.matrix(param.list$covariate_y)
#'   covariate_x <- as.matrix(param.list$covariate_x)
#'   y_hat_alpha_beta <- as.matrix(param.list$y_hat_alpha_beta)
#'   alpha.L <- as.matrix(param.list$alpha.L)
#'   alpha.L <- as.matrix(param.list$alpha.L)
#'   beta <- as.matrix(param.list$beta)
#'   delta <- as.matrix(param.list$delta)
#'   gamma = as.matrix(param.list$gamma)
#'   
#'   return(list(
#'     x = as.matrix((x)),
#'     y = as.matrix((y)),
#'     covariate_x = as.matrix(covariate_x),
#'     covariate_y = as.matrix(covariate_y),
#'     fitted = as.matrix((y_hat_alpha_beta)),
#'     residuals = as.matrix((as.matrix((y)) - (y_hat_alpha_beta))),
#'     alpha.L = as.matrix(alpha.L),
#'     alpha.L = as.matrix(alpha.L),
#'     beta = as.matrix(beta),
#'     delta = as.matrix(delta),
#'     gamma = as.matrix(gamma)
#'   )
#'   )
#' }
#' 
#' initialise_param.list <- function(x,y,covariate_y,alpha.init=NULL,beta.init=NULL, k_dim = 10){
#'   # Prepare y
#'   
#'   yy <- y
#'   y <- log(1+y)
#'   yT <- Matrix::t(y)
#'   y <- Matrix::Matrix(y)
#'   
#'   # Prepare x
#'   xT <- Matrix::t(x)
#'   x <- Matrix::Matrix(x)
#'   
#'   if (is.null(covariate_y)){
#'     covariate_y = as.vector(rep(0,dim(y)[2]))
#'   }
#'   if (is.vector(covariate_y)){
#'     n = 1
#'     p = length(covariate_y)
#'   }else{
#'     n = dim(covariate_y)[1]
#'     p = dim(covariate_y)[2]
#'   }
#'   
#'   covariate_y <- Matrix::Matrix(matrix(covariate_y,nrow = n,ncol=p))
#'   
#'   if (is.null(covariate_x)){
#'     covariate_x = as.vector(rep(0,dim(x)[2]))
#'   }
#'   if (is.vector(covariate_x)){
#'     n = 1
#'     p = length(covariate_x)
#'   }else{
#'     n = dim(covariate_x)[1]
#'     p = dim(covariate_x)[2]
#'   }
#'   
#'   covariate_x <- Matrix::Matrix(matrix(covariate_x,nrow = n,ncol=p))
#'   
#'   if (is.null(beta.init)){
#'     beta = MASS::ginv(as.matrix(Matrix::t(x)%*%(x)))%*%Matrix::t((x))%*%(y)
#'   } else {
#'     beta = beta.init
#'   }
#'   
#'   if (is.null(alpha.init)){
#'     alpha = (y)%*%Matrix::t(beta)%*%MASS::ginv(as.matrix(beta%*%Matrix::t(beta)))%*%xT%*%MASS::ginv(as.matrix(x%*%xT))
#'   } else {
#'     alpha = alpha.init
#'   }
#'   
#'   alpha.L <- Matrix::t(Matrix::Matrix(RSpectra::eigs(as.matrix(alpha%*%Matrix::t(alpha)),k=k_dim)$vectors))
#'   
#'   delta = (alpha.L%*%y - alpha.L%*%x%*%beta)%*%Matrix::t(covariate_y)%*%MASS::ginv(as.matrix(covariate_y%*%Matrix::t(covariate_y)))
#'   
#'   # Calculate fitted value of y; covariance of y ; likelihood ; mean-squared error
#'   y_hat_alpha_beta <- as.matrix(x)%*%beta + MASS::ginv(as.matrix(Matrix::t(alpha.L)%*%alpha.L))%*%Matrix::t(alpha.L)%*%(delta%*%covariate_y)
#'   D = Matrix::diag(Matrix::diag(Matrix::t((y) - (y_hat_alpha_beta))%*%((y) - (y_hat_alpha_beta))/dim(y)[1]))
#'   
#'   param.list <- list(x = x,
#'                      y = yy,
#'                      D = D,
#'                      alpha.L = alpha.L,
#'                      beta = beta,
#'                      delta = delta,
#'                      covariate_y = covariate_y)
#'   
#'   return(param.list)
#'   
#' }
#' 
#' calculate_scores_and_update_params <- function(param.list,alpha.L_run=T,beta_run=T,delta_run=T,D_run=T){
#'   
#'   # Prepare main parameters
#'   x <- param.list$x
#'   y <- param.list$y
#'   D <- param.list$D
#'   covariate_y <- param.list$covariate_y
#'   alpha.L <- param.list$alpha.L
#'   beta <- param.list$beta
#'   delta <- param.list$delta
#'   
#'   
#'   mu <- exp(alpha.L%*%x%*%beta + delta%*%covariate_y)
#'   mu[mu>10e5] <- 10e5
#'   mu[mu<1] <- 1
#'   
#'   if (alpha.L_run){
#'     first.alpha.L <- Matrix::t(y %*% log(D%*%Matrix::t(mu)) - y %*% (log((D%*%Matrix::t(mu)+1))) -
#'                                  x%*%beta%*%Matrix::t(mu-alpha.L%*%y)%*%
#'                                  MASS::ginv(as.matrix((mu)%*%Matrix::t(mu%*%D+1))))
#'     
#'     info.mat.alpha.L <- MASS::ginv(as.matrix(first.alpha.L%*%Matrix::t(first.alpha.L)))
#'     alpha.L <- alpha.L + info.mat.alpha.L%*%first.alpha.L
#'     param.list$alpha.L <- alpha.L
#'     
#'   }
#'   
#'   
#'   mu <- exp(alpha.L%*%x%*%beta + delta%*%covariate_y)
#'   mu[mu>10e5] <- 10e5
#'   mu[mu<1] <- 1
#'   
#'   if (beta_run){
#'     first.beta <- Matrix::t(alpha.L%*%x)%*%(alpha.L%*%x%*%beta + delta%*%covariate_y + alpha.L%*%y)%*%
#'       MASS::ginv(as.matrix(Matrix::t(alpha.L%*%x%*%beta+delta%*%covariate_y)%*%((alpha.L%*%x%*%beta+delta%*%covariate_y)%*%D+1)))
#'     
#'     info.mat.beta <- MASS::ginv(as.matrix(first.beta%*%Matrix::t(first.beta)))
#'     beta <- beta + info.mat.beta%*%first.beta
#'     param.list$beta <- beta
#'     
#'   }
#'   
#'   
#'   mu <- exp(alpha.L%*%x%*%beta + delta%*%covariate_y)
#'   mu[mu>10e5] <- 10e5
#'   mu[mu<1] <- 1
#'   
#'   if (delta_run){
#'     
#'     first.delta <- Matrix::t((covariate_y)%*%Matrix::t(alpha.L%*%x%*%beta + delta%*%covariate_y + alpha.L%*%y)%*%
#'                                MASS::ginv(as.matrix((alpha.L%*%x%*%beta+delta%*%covariate_y)%*%Matrix::t((alpha.L%*%x%*%beta+delta%*%covariate_y)%*%D+1))))
#'     
#'     info.mat.delta <- MASS::ginv(as.matrix(first.delta%*%Matrix::t(first.delta)))
#'     delta <- delta + info.mat.delta%*%first.delta
#'     param.list$delta <- delta
#'   }
#'   
#'   
#'   
#'   mu <- exp(alpha.L%*%x%*%beta + delta%*%covariate_y)
#'   mu[mu>10e5] <- 10e5
#'   mu[mu<1] <- 1
#'   if (D_run){
#'     
#'     first.D <- ((((MASS::ginv(as.matrix(D%*%D)))%*%diag(as.vector(c(rowSums(as.matrix(log(1+MASS::ginv(as.matrix(D))%*%Matrix::t(mu))))) -
#'                                                                     c(Reduce('+',lapply(c(1:dim(y)[1]),function(y_ID){
#'                                                                       do.call('c',lapply(c(1:length(y[y_ID,])),function(YY_ID){
#'                                                                         sum(do.call('c',lapply(c(1:y[y_ID,YY_ID]),function(YYYY){
#'                                                                           1/(YYYY+Matrix::diag(D)[YY_ID])
#'                                                                         })))
#'                                                                       }))
#'                                                                     }))/prod(dim(y)))
#'     )))) +
#'       (Matrix::t(alpha.L%*%y - mu)%*%MASS::ginv(as.matrix(MASS::ginv(as.matrix(D))%*%(1+MASS::ginv(as.matrix(D))%*%Matrix::t(mu)))))
#'     )
#'     
#'     info.mat.D <- MASS::ginv(as.matrix(first.D%*%Matrix::t(first.D)))
#'     D <- D + Matrix::diag(Matrix::diag(info.mat.D%*%first.D))
#'     param.list$D <- D
#'     
#'   }
#'   
#'   
#'   # param.list$mllik <- sum(alpha.L%*%y%*%(log((D%*%Matrix::t(mu)))) - alpha.L%*%y%*%(log((1+D%*%Matrix::t(mu)))) - Matrix::t(MASS::ginv(as.matrix(D))%*%(log((1+D%*%Matrix::t(mu)))))) + sum(do.call('c',lapply(dim(y)[1],function(y_ID){y[y_ID,]+Matrix::diag(D)})))
#'   param.list$mllik <- 0
#'   param.list$mscore <- mean(sqrt(as.matrix(alpha.L%*%y - mu)^2))
#'   
#'   return(param.list)
#' }
