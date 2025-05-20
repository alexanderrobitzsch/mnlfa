## File Name: mnlfa_mstep_item_loglike_update_parameters.R
## File Version: 0.337


mnlfa_mstep_item_loglike_update_parameters <- function(y, y_resp, theta, parms,
        Xdes_int, Xdes_slo, Xdes_res, post, b_index, a_index, psi_index, K_index,
        parms_indices, h, N_item, regular_type, regular_lam, regular_alpha,
        center_group_parms,    item_type, eps=1E-15, L_max=.25, numdiff=TRUE,
        mstep_shortcut=FALSE, weights=NULL)
{

    res <- mnlfa_mstep_item_loglike_deriv( y=y, y_resp=y_resp, theta=theta,
                parms=parms, Xdes_int=Xdes_int, Xdes_slo=Xdes_slo, Xdes_res=Xdes_res,
                post=post, b_index=b_index, a_index=a_index, psi_index=psi_index,
                K_index=K_index, parms_indices=parms_indices, h=h, N_item=N_item, eps=eps,
                numdiff=numdiff, mstep_shortcut=mstep_shortcut, recycle=recycle,
                weights=weights, item_type=item_type )
    incr <- res$incr
    L <- max( res$D2_max, L_max)
    ll <- res$ll
    recycle <- res$recycle
    # Newton-Raphson update
    parms[ parms_indices ] <- parms[ parms_indices ] + incr

    # apply threshold operator
    parms <- mnlfa_parameter_regularization( parms=parms, parms_indices=parms_indices,
                    L=L, regular_type=regular_type, regular_lam=regular_lam,
                    regular_tau=NULL, regular_alpha=regular_alpha,
                    center_group_parms=center_group_parms )
    #--- output
    res <- list( parms=parms, ll=ll, recycle=recycle)
    return(res)
}
