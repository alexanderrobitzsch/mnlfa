## File Name: mnlfa_proc_data.R
## File Version: 0.224

mnlfa_proc_data <- function(dat, items, item_type, formula_mean, formula_sd,
            parm_trait_init, weights=NULL)
{
    resp <- dat[,items]
    N <- nrow(resp)
    I <- ncol(resp)
    if (length(item_type)==1 ){
        item_type <- rep(item_type, I)
        names(item_type) <- items
    }
    resp_ind <- ! is.na(resp)
    N_item <- colSums(resp_ind)
    resp[ is.na(resp) ] <- 0

    # design matrix predictors trait distribution
    Xdes_mean <- stats::model.matrix( object=formula_mean, data=dat)
    Xdes_sd <- stats::model.matrix( object=formula_sd, data=dat)
    mu <- rep(0,ncol(Xdes_mean))
    n_mu <- length(mu)
    if (n_mu>0){
        names(mu) <- paste0('mu_', colnames(Xdes_mean))
    }
    sigma <- rep(0,ncol(Xdes_sd))
    n_sigma <- length(sigma)
    if (n_sigma>0){
        names(sigma) <- paste0('sigma_', colnames(Xdes_sd))
    }
    if ( is.null(parm_trait_init) ){
        parm_trait <- list(mu=mu, sigma=sigma)
        parm_trait$index$mu <- seq_len(n_mu)
        parm_trait$index$sigma <- n_mu + seq_len( length(sigma) )
        parm_trait$npar <- length( unlist( parm_trait$index ))
    } else {
        parm_trait <- parm_trait_init
    }
    if (is.null(weights)){
        weights <- rep(1,N)
    }
    K <- as.numeric( apply(resp, 2, max, na.rm=TRUE) )

    #-- output
    res <- list(resp=resp, N=N, I=I, resp_ind=resp_ind, item_type=item_type,
                    N_item=N_item, Xdes_mean=Xdes_mean, Xdes_sd=Xdes_sd,
                    parm_trait=parm_trait, weights=weights, K=K)
    return(res)
}
