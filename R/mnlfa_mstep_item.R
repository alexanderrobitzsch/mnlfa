## File Name: mnlfa_mstep_item.R
## File Version: 0.609


mnlfa_mstep_item <- function(y, y_resp, theta, parms, Xdes_int,
        Xdes_slo, Xdes_res, post, b_index, a_index, psi_index, K_index, parms_iterations,
        h, N_item, parms_regular_types, parms_regular_lam, parms_regular_alpha,
        center_group_parms, msteps,    conv_mstep, item_type, eps=1E-15, L_max=.25,
        numdiff=TRUE, mstep_shortcut=FALSE,    weights=NULL)
{
    NH <- length(parms_iterations)
    iterate <- TRUE
    iter <- 1

    #--- begin iterations M-steps
    while(iterate){
        parms0 <- parms
        recycle <- NULL
        fct_mstep <- 'mnlfa_mstep_item_loglike_update_parameters'

        for (hh in seq_len(NH) ){
            parms_indices <- parms_iterations[[hh]]
            regular_type <- parms_regular_types[parms_indices][1]
            regular_lam <- parms_regular_lam[parms_indices][1]
            regular_alpha <- parms_regular_alpha[parms_indices][1]
            args <- list( y=y, y_resp=y_resp, theta=theta, parms=parms, Xdes_int=Xdes_int,
                        Xdes_slo=Xdes_slo, Xdes_res=Xdes_res, post=post, b_index=b_index,
                        a_index=a_index, psi_index=psi_index, K_index=K_index,
                        parms_indices=parms_indices, h=h, N_item=N_item,
                        regular_type=regular_type, regular_lam=regular_lam,
                        regular_alpha=regular_alpha,
                        center_group_parms=center_group_parms, eps=eps, L_max=L_max,
                        numdiff=numdiff, mstep_shortcut=mstep_shortcut,
                        weights=weights, item_type=item_type )
            res <- do.call( what=fct_mstep, args=args )
            recycle <- res$recycle
            parms <- res$parms

        }
        change <- max( abs( parms - parms0 ) )
        iter <- iter + 1
        iterate <- ! ( ( iter > msteps ) | ( change < conv_mstep ) )
    }
    #--- end iterations M-steps
    iter_mstep <- iter - 1

    #*** penalty values
    res <- mnlfa_penalty_values_item( parms=parms, parms_iterations=parms_iterations,
                    parms_regular_types=parms_regular_types,
                    parms_regular_lam=parms_regular_lam,
                    parms_regular_alpha=parms_regular_alpha,
                    N_item=N_item, center_group_parms=center_group_parms )
    parms_values <- res$parms_values
    parms_regularized <- res$parms_regularized
    parms_estimated <- res$parms_estimated

    #--- output
    res <- list( parms=parms, iter_mstep=iter_mstep, parms_change=change,
                        parms_values=parms_values, parms_regularized=parms_regularized,
                        parms_estimated=parms_estimated )
    return(res)
}
