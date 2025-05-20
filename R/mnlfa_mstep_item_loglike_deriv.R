## File Name: mnlfa_mstep_item_loglike_deriv.R
## File Version: 0.286


mnlfa_mstep_item_loglike_deriv <- function(y, y_resp, theta, parms, Xdes_int,
        Xdes_slo, Xdes_res, post, b_index, a_index, psi_index, K_index,
        parms_indices, h, N_item, item_type, eps=1E-15,
        numdiff=TRUE, mstep_shortcut=FALSE, recycle=NULL, weights=NULL)
{

    NP <- length(parms_indices)
    if (NP>1){
        numdiff <- TRUE
    }

    ll0 <- NULL

    if (numdiff){

        args <- list( y=y, y_resp=y_resp, theta=theta, parms=parms, Xdes_int=Xdes_int,
                            Xdes_slo=Xdes_slo, Xdes_res=Xdes_res, post=post,
                            b_index=b_index, a_index=a_index, psi_index=psi_index,
                            K_index=K_index, N_item=N_item, eps=eps, weights=weights,
                            item_type=item_type )
        fct <- 'mnlfa_mstep_item_loglike_eval'

        ll0 <- do.call( what=fct, args=args)

        ll1 <- rep(NA,NP)
        ll2 <- rep(NA,NP)
        parms0 <- parms

        for (pp in 1L:NP){
            i_pp <- parms_indices[pp]
            args$parms <- mnlfa_add_increment(x=parms0, i=i_pp, h=h)
            ll1[pp] <- do.call( what=fct, args=args)
            args$parms <- mnlfa_add_increment(x=parms0, i=i_pp, h=-h)
            ll2[pp] <- do.call( what=fct, args=args)
        }
        #-- compute derivatives
        res <- mnlfa_differences_derivative(ll0=ll0, ll1=ll1, ll2=ll2, h=h)
        res$recycle <- NULL

    }


    #-- analytical derivative

    if (!numdiff){

        res <- mnlfa_mstep_item_loglike_deriv_analytical( parms=parms,
                        a_index=a_index, b_index=b_index, K_index=K_index,
                        psi_index=psi_index, Xdes_int=Xdes_int, Xdes_slo=Xdes_slo,
                        Xdes_res=Xdes_res,
                        theta=theta, parms_indices=parms_indices,
                        N_item=N_item, y=y, y_resp=y_resp, post=post,
                        weights=weights, item_type=item_type)

    }

    D1 <- res$D1
    D2 <- res$D2
    D2_max <- res$D2_max
    recycle <- res$recycle

    #-- updates
    incr <- D1 / D2_max
    #-- output
    res <- list( ll=ll0, D1=D1, D2=D2, D2_max=D2_max, incr=incr, recycle=recycle)
    return(res)
}
