## File Name: mnlfa.R
## File Version: 0.767

mnlfa <- function( dat, items, weights=NULL, item_type="2PL", formula_int=~1,
    formula_slo=~1, formula_res=~1, formula_mean=~0, formula_sd=~0, theta=NULL,
    parm_list_init=NULL, parm_trait_init=NULL, prior_init=NULL,
    regular_lam=c(0,0,0), regular_alpha=c(0,0,0), regular_type=c("none","none","none"),
    maxit=1000, msteps=4, conv=1E-5, conv_mstep=1E-4, h=1E-4,
    parms_regular_types=NULL, parms_regular_lam=NULL, parms_regular_alpha=NULL,
    parms_iterations=NULL, center_parms=NULL, center_max_iter=6,
    L_max=.07, verbose=TRUE, numdiff=FALSE)
{
    CALL <- match.call()
    s1 <- Sys.time()
    mstep_shortcut <- FALSE

    #- theta design matrix
    res <- mnlfa_proc_theta_design(theta=theta)
    theta <- res$theta
    TP <- res$TP

    #- data processing
    res <- mnlfa_proc_data( dat=dat, items=items, item_type=item_type,
                formula_mean=formula_mean, formula_sd=formula_sd,
                parm_trait_init=parm_trait_init, weights=weights )
    resp <- res$resp
    N <- res$N
    I <- res$I
    resp_ind <- res$resp_ind
    item_type <- res$item_type
    N_item <- res$N_item
    Xdes_mean <- res$Xdes_mean
    Xdes_sd <- res$Xdes_sd
    parm_trait <- res$parm_trait
    weights <- res$weights
    K <- res$K

    #- inits item parameters
    res <- mnlfa_proc_item_parameters( dat=dat, formula_int=formula_int,
                formula_slo=formula_slo, formula_res=formula_res, item_type=item_type,
                parms_regular_types=parms_regular_types,
                parms_regular_lam=parms_regular_lam,
                parms_regular_alpha=parms_regular_alpha, regular_type=regular_type,
                regular_lam=regular_lam, regular_alpha=regular_alpha,
                parms_iterations=parms_iterations,
                parm_list_init=parm_list_init, K=K)
    parm_list <- res$parm_list
    parm_Xdes <- res$parm_Xdes
    parms_regular_types <- res$parms_regular_types
    parms_regular_lam <- res$parms_regular_lam
    parms_regular_alpha <- res$parms_regular_alpha
    parms_iterations <- res$parms_iterations
    parm_index <- res$parm_index

    #* inits prior distribution
    if ( is.null(prior_init) ){
        prior <- mnlfa_proc_inits_prior(theta=theta, N=N)
    } else {
        prior <- prior_init
    }

    center_group_parms <- TRUE

    disp <- '----------------------------------------------\n'
    iter <- 1
    eps0 <- 1E-50
    iterate <- TRUE

    #----------------------------------
    #--- begin EM algorithm
    while (iterate){

        zz0 <- Sys.time()
        parm_list0 <- parm_list
        parm_trait0 <- parm_trait

        res <- mnlfa_compute_likelihood( resp=resp, resp_ind=resp_ind,
                            parm_Xdes=parm_Xdes, parm_list=parm_list, prior=prior,
                            item_type=item_type, theta=theta)
        like <- res$like
        post <- like / rowSums(like)
        loglike_case <- log( rowSums(like) + eps0 )
        ll <- sum( weights*loglike_case )
        dev <- -2*ll
        if (verbose){
            cat(disp)
            cat('Iteration', iter, '\n')
            cat( 'Deviance','=', round(dev, 6) )
            utils::flush.console()
            cat('\nM-step item parameters ')
        }

        #* M-step item parameters
        parms_values <- list()
        parms_regularized <- list()
        parms_estimated <- list()
        for (ii in 1L:I){
            y <- resp[,ii]
            y_resp <- resp_ind[,ii]
            parms <- c(parm_list[[ii]]$b, parm_list[[ii]]$a, parm_list[[ii]]$psi,
                            parm_list[[ii]]$K)
            b_index <- parm_index[[ii]]$b_index
            a_index <- parm_index[[ii]]$a_index
            psi_index <- parm_index[[ii]]$psi_index
            K_index <- parm_index[[ii]]$K_index
            Xdes_int <- parm_Xdes[[ii]]$Xdes_int
            Xdes_slo <- parm_Xdes[[ii]]$Xdes_slo
            Xdes_res <- parm_Xdes[[ii]]$Xdes_res
            item_type_ii <- item_type[ii]
            res <- mnlfa_mstep_item( y=y, y_resp=y_resp, theta=theta, parms=parms,
                        Xdes_int=Xdes_int, Xdes_slo=Xdes_slo, Xdes_res=Xdes_res,
                        post=post, b_index=b_index, a_index=a_index, psi_index=psi_index,
                        K_index=K_index, parms_iterations=parms_iterations[[ii]], h=h,
                        N_item=N_item[ii], parms_regular_types=parms_regular_types[[ii]],
                        parms_regular_lam=parms_regular_lam[[ii]],
                        parms_regular_alpha=parms_regular_alpha[[ii]],
                        center_group_parms=center_group_parms, msteps=msteps,
                        conv_mstep=conv_mstep, eps=1E-15, L_max=L_max, numdiff=numdiff,
                        mstep_shortcut=mstep_shortcut, weights=weights,
                        item_type=item_type_ii)
            parm_list[[ii]]$b <- res$parms[ b_index ]
            parm_list[[ii]]$a <- res$parms[ a_index ]
            parm_list[[ii]]$psi <- res$parms[ psi_index ]
            parms_values[[ii]] <- res$parms_values
            parms_regularized[[ii]] <- res$parms_regularized
            parms_estimated[[ii]] <- res$parms_estimated

            # if (verbose){
            if (FALSE){
                cat('.')
                utils::flush.console()
            }
        }

        #* parameter centering
        if ( ( ! is.null(center_parms) ) & ( iter <=center_max_iter )  ){
            NL <- length(center_parms)
            for (ll in seq_len(NL) ){
                center_ll <- center_parms[[ll]]
                parms_center <- rep(NA, I)
                for (ii in 1L:I){
                    parms <- c(parm_list[[ii]]$b, parm_list[[ii]]$a )
                    parms_center[ii] <- parms[center_ll[ii]]
                }
                mpc <- stats::weighted.mean(parms_center, w=weights, na.rm=TRUE)
                parms_center <- parms_center - mpc
                for (ii in 1L:I){
                    parms <- c(parm_list[[ii]]$b, parm_list[[ii]]$a )
                    parms[ center_ll[ii] ] <- parms_center[ii]
                    b_index <- parm_index[[ii]]$b_index
                    a_index <- parm_index[[ii]]$a_index
                    parm_list[[ii]]$b <- parms[ b_index ]
                    parm_list[[ii]]$a <- parms[ a_index ]

                    res <- mnlfa_penalty_values_item( parms=parms,
                                parms_iterations=parms_iterations[[ii]],
                                parms_regular_types=parms_regular_types[[ii]],
                                parms_regular_lam=parms_regular_lam[[ii]],
                                N_item=N_item[[ii]],
                                center_group_parms=center_group_parms )
                    parms_values[[ii]] <- res$parms_values
                    parms_regularized[[ii]] <- res$parms_regularized
                    parms_estimated[[ii]] <- res$parms_estimated

                }
            }
        }

        # change in item parameters
        parm_change_item <- max( abs( unlist(parm_list0) - unlist(parm_list) ) )
        if (verbose){
            cat('\nMaximum change in item parameters', '=',
                        round(parm_change_item,6), '\n')
            utils::flush.console()
        }

        #* estimation of trait distribution
        mu <- parm_trait$mu
        sigma <- parm_trait$sigma
        index_mu <- parm_trait$index$mu
        index_sigma <- parm_trait$index$sigma
        fun_trait_regression_unidim <- function(x)
        {
            mu <- x[ index_mu ]
            sigma <- x[ index_sigma ]
            mu_p <- Xdes_mean %*% mu
            sigma_p <- exp( Xdes_sd %*% sigma )
            val <- - mnlfa_rcpp_mstep_trait_unidim_fun( theta=theta, mu_p=mu_p,
                            sigma_p=sigma_p, post=post, weights=weights )
            return(val)
        }
        par <- c(mu, sigma)
        n_parm_trait <- length(par)
        use_rcpp <- TRUE
        #    use_rcpp <- FALSE

        grad_trait_regression_unidim <- function(x)
        {
            grad <- mnlfa_grad_trait_regression_unidim(x=par, post=post, theta=theta,
                        Xdes_mean=Xdes_mean, Xdes_sd=Xdes_sd,
                        index_mu=index_mu, index_sigma=index_sigma,
                        weights=weights, use_rcpp=use_rcpp)
            return(grad)
        }

        objective <- fun_trait_regression_unidim
        gradient <- grad_trait_regression_unidim
        # gradient <- NULL

        # optimizer <- 'nlminb'
        optimizer <- 'optim'

        if (parm_trait$npar>0){

            if (optimizer=='nlminb'){
                res <- stats::nlminb( start=par, objective=objective, gradient=gradient,
                                        control=list(maxit=msteps) )
            }
            if (optimizer=='optim'){
                res <- stats::optim( par=par, fn=objective, gr=gradient,
                                        method='BFGS', control=list(maxit=msteps) )
            }

            parm_trait$mu <- res$par[ index_mu ]
            parm_trait$sigma <- res$par[ index_sigma ]

        }

        parm_change_trait <- max( abs( unlist(parm_trait0) - unlist(parm_trait) ) )
        parm_change <- max(parm_change_item, parm_change_trait)
        if (verbose){
            cat('Maximum change in trait distribution parameters', '=',
                            round(parm_change_trait,6), '\n')
            utils::flush.console()
        }

        #* compute prior distribution
        prior <- mnlfa_compute_prior( parm_trait=parm_trait, Xdes_mean=Xdes_mean,
                        Xdes_sd=Xdes_sd, N=N, theta=theta )

        #* penalty values in optimization
        regular_penalty <- sum( unlist(parms_values) )
        opt_fct <- dev/2 + regular_penalty
        numb_reg_pars <- sum( unlist(parms_regularized) )
        numb_est_pars <- sum( unlist(parms_estimated) ) + n_parm_trait
        BIC <- dev + log(N)*numb_est_pars
        AIC <- dev + 2*numb_est_pars
        if (verbose){
            cat( 'Number of estimated parameters','=', numb_est_pars, '\n')
            cat( 'Number of regularized parameters','=', numb_reg_pars, '\n')
            cat( 'Optimization function value','=', round(opt_fct, 6), '\n')
            cat( 'AIC','=', round(AIC, 2), ' | ')
            cat( 'BIC','=', round(BIC, 2), '\n')
            utils::flush.console()
        }

        iter <- iter + 1
        if ( iter > maxit ){ iterate <- FALSE }
        if ( parm_change < conv ){ iterate <- FALSE }
    }
    #--- end EM algorithm
    #-----------------------------

    converged <- ( iter <=maxit )

    #--- collect all item parameters
    item <- mnlfa_postproc_item( parm_list=parm_list, items=items, item_type=item_type,
                parms_estimated=parms_estimated, parms_regularized=parms_regularized,
                parms_regular_types=parms_regular_types,
                parms_regular_lam=parms_regular_lam,
                parms_regular_alpha=parms_regular_alpha)

    #-- collect parameters of trait distribution
    trait <- mnlfa_postproc_trait(parm_trait=parm_trait)

    #-- information criteria
    W <- sum(weights)
    ic <- list(deviance=dev, np=numb_est_pars, n=N, numb_reg_pars=numb_reg_pars,
                    numb_est_pars=numb_est_pars, W=W)
    ic <- CDM::cdm_calc_information_criteria(ic)

    #--- output
    s2 <- Sys.time()
    time <- list(s1=s1, s2=s2, diff_time=s2-s1)
    res <- list( item=item, trait=trait, parm_list=parm_list, parm_trait=parm_trait,
                    dev=dev, numb_est_pars=numb_est_pars, numb_reg_pars=numb_reg_pars,
                    ic=ic, deviance=dev, prior=prior, post=post,
                    loglike=ll, loglike_case=loglike_case,
                    parms_regular_types=parms_regular_types,
                    parms_regular_lam=parms_regular_lam,
                    parms_iterations=parms_iterations,
                    parm_index=parm_index, regular_penalty=regular_penalty,
                    parms_values=parms_values, parms_regularized=parms_regularized,
                    parms_estimated=parms_estimated, iter=iter-1, converged=converged,
                    call=CALL, time=time )
    class(res) <- 'mnlfa'
    return(res)
}


# cat('mstep trait 2') ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1

