## File Name: mnlfa_mstep_item_loglike_eval.R
## File Version: 0.156


mnlfa_mstep_item_loglike_eval <- function(y, y_resp, theta, parms, Xdes_int,
        Xdes_slo, Xdes_res, post, b_index, a_index, psi_index, K_index, N_item,
        item_type, eps=1E-30, weights=NULL )
{
    b_ii <- parms[ b_index ]
    a_ii <- parms[ a_index ]
    if (item_type %in% c('1PL','2PL','NO')){
        b <- Xdes_int %*% b_ii
    }

    if (item_type %in% c('GPCM')){
        K_ii <- parms[ K_index ]
        parms_ii <- list(b=b_ii, K=K_ii)
        b <- mnlfa_compute_moderated_parameter_gpcm_b(Xdes=Xdes_int,
                    parm=parms_ii, value=0)
    }
    a <- Xdes_slo %*% a_ii

    args <- list( a=a, b=b, theta=theta, y=y, y_resp=y_resp )

    if (item_type %in% c('1PL','2PL')){
        fun <- 'mnlfa_rcpp_calc_probs_2pl'
    }
    if (item_type=='NO'){
        psi_ii <- parms[ psi_index ]
        args$psi <- exp( Xdes_res %*% psi_ii )
        fun <- 'mnlfa_rcpp_calc_dnorm'
    }
    if (item_type %in% c('GPCM')){
        fun <- 'mnlfa_rcpp_calc_probs_gpcm'
    }

    like_ii <- do.call(what=fun, args=args)
    log_probs <- log(like_ii+eps)
    ll <- sum( weights*post * log_probs )
    ll <- ll / N_item

    return(ll)
}
