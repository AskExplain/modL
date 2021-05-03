deepgmm <- function(y, k, r, W = 10, layers = 3,
            it = 250, eps = 0.001, init = 'kmeans', init_est = 'factanal',
            seed = NULL, scale = TRUE) {


  if (any(class(y) %in% 'data.frame'))
  	y <- as.matrix(y)

  if (!is.null(seed)) {

    if(!is.numeric(seed)) {
      stop("The value of seed must be an integer")
    }

    set.seed(seed)
  }

  if (scale)
    y <- scale(y)

  if (any(tolower(init) %in% c('kmeans', 'k-means', 'k')))
    init <- 'kmeans'
  if (any(tolower(init) %in% c('random', 'r')))
    init <- 'random'
  if (any(tolower(init) %in% c('hclass', 'h')))
    init <- 'hclass'
  if (any(tolower(init_est) == c('factanal', 'factana', 'fact', 'f')))
    init_est <- 'factanal'


  # check arguments
  # tmp <- valid_args(Y = y, layers = layers, k = k, r = r, it = it,
  #                   eps = eps, init = init)

  all_converged <- c(F,F)
  hh <- 0
  out.lst <- list(out = list(), ps.y.list = list(), ps.y = list(), ps.y.main = list(), lik = list(), bic = list())
  IDS <- c()
  comb_main <- c()
  main_IDS <- rep(T,dim(y)[1])

  while(T){
    hh <- hh + 1

    if (all_converged[1]==F){
      orig_data = y

    } else {
      orig_data = init_data

    }

    lst <- list(w = list(), H = list(), mu = list(), psi = list(),
                psi.inv = list())

    numobs <- nrow(orig_data)
    p <- ncol(orig_data)

    if (all_converged[1]==F){
      r <- c(p, r)
    } else {
      r <- c(p,p-1,p-2)
    }

    print(tail(r,2))
    i_lst <- deepgmm::deepgmm(y = y,layers = 2,k = k,r = r[-1], it = 3)

    lst$w <- i_lst$w
    lst$H <- i_lst$H
    lst$mu <- i_lst$mu
    lst$psi <- i_lst$psi
    lst$psi.inv <- lapply(c(1:length(i_lst$psi)),function(X){
      for (j in 1:k[X]){
        i_lst$psi[[X]][j,,]<-diag(1/diag(i_lst$psi[[X]][j,,]))
      }
      return(i_lst$psi[[X]])
    })


    if (all_converged[1]==T){
      i=1
      lst$w[i] <- list(out.lst$out[[hh-1]]$w[[2]])
      lst$H[i] <- list(out.lst$out[[hh-1]]$H[[2]])
      lst$mu[i] <- list(out.lst$out[[hh-1]]$mu[[2]])
      lst$psi[i] <- list(out.lst$out[[hh-1]]$psi[[2]])
      lst$psi.inv[i] <- list(out.lst$out[[hh-1]]$psi.inv[[2]])

    }

    if (all_converged[1]==F){
      lst$llik <- compute.lik(orig_data, numobs, k, lst$mu, lst$H, lst$psi, lst$w, main_IDS)

    } else {
      lst$llik <- init_llik_c

    }

    out <- deep.sem.alg.2(orig_data, numobs, p, r, k, lst$H, lst$psi,
                          lst$psi.inv, lst$mu, lst$w, lst$llik$ps.y.list,it, eps, main_IDS)

    init_llik_c <- out$ps.y.list
    init_data <- out$Ez
    all_converged[hh] <- T

    comb_main <- cbind(comb_main,out$s)
    main <- as.integer(as.factor(apply(comb_main,1,function(X){paste0(X,collapse="_")})))
    print(table(main))

    main_IDS <- c(main_IDS & main%in%names(table(main))[which(table(main)>W)])
    sum(main_IDS)
    # if (length(names(table(main))[which(table(main)<W)])>0){
    #   IDS <- unique(c(IDS,which(main_IDS==T)))
    #   init_llik_c <- lapply(init_llik_c,function(X){
    #     X[IDS,]<-0
    #     return(X)
    #   })
    #   init_data[IDS,] <- 0
    # }


    if (length(names(table(main))[which(table(main)<W)])!=length(unique(main))){
      all_converged <- c(all_converged,F)
    }

    out.lst$out[[hh]] <- out
    out.lst$ps.y.list[[hh]] <- out$ps.y.list
    out.lst$ps.y[[hh]] <- out$ps.y
    out.lst$main[[hh]] <- out$s
    out.lst$lik[[hh]] <- c(out$lik)
    out.lst$bic[[hh]] <- c(out$bic)


    if (hh>1){
      ps.y.main <- c()
      for (t in c(1:dim(out.lst$ps.y.main)[2])){
        ps.y.main<-cbind(ps.y.main,out.lst$ps.y.main[,t]*out$ps.y)
      }
      out.lst$ps.y.main <- ps.y.main
    } else {
      out.lst$ps.y.main <- out$ps.y
    }


    if (all(!main_IDS) | tail(r,1)==1 | hh >= layers){
      break
    }
  }

  out.lst$comb_main <- comb_main
  out.lst$main <- main
  out.lst$layers <- hh

  return(out.lst)

}
