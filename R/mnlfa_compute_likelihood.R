## File Name: mnlfa_compute_likelihood.R
## File Version: 0.098


mnlfa_compute_likelihood <- function(resp, resp_ind, parm_Xdes, parm_list, prior,
        item_type, theta)
{
    N <- nrow(resp)
    TP <- nrow(theta)
    #* compute likelihood
    like <- matrix(1, nrow=N, ncol=TP)
    for (ii in 1L:I){
        item_type_ii <- item_type[ii]
        if (item_type_ii %in% c('1PL','2PL','NO') ){
            b <- mnlfa_compute_moderated_parameter(Xdes=parm_Xdes[[ii]]$Xdes_int,
                    parm=parm_list[[ii]]$b, value=0)
        }
        if (item_type_ii %in% c('GPCM') ){
            b <- mnlfa_compute_moderated_parameter_gpcm_b(Xdes=parm_Xdes[[ii]]$Xdes_int,
                        parm=parm_list[[ii]], value=0)
        }
        a <- mnlfa_compute_moderated_parameter(Xdes=parm_Xdes[[ii]]$Xdes_slo,
                    parm=parm_list[[ii]]$a, value=1)
        if (item_type_ii=='NO'){
            psi <- mnlfa_compute_moderated_parameter(Xdes=parm_Xdes[[ii]]$Xdes_res,
                        parm=parm_list[[ii]]$psi, value=0, trafo='exp')
        }
        y <- resp[,ii]
        y_resp <- resp_ind[,ii]

        args <- list( a=a, b=b, theta=theta, y=y, y_resp=y_resp )

        if (item_type_ii %in% c('1PL','2PL') ){
            fun <- 'mnlfa_rcpp_calc_probs_2pl'
        }

        if (item_type_ii=='NO'){
            args$psi <- psi
            fun <- 'mnlfa_rcpp_calc_dnorm'
        }

        if (item_type_ii=='GPCM'){
            fun <- 'mnlfa_rcpp_calc_probs_gpcm'
        }

        like_ii <- do.call(what=fun, args=args)
        like <- like * like_ii

    }  # end ii
    like_ind <- like
    like <- like * prior

    #-- output
    res <- list(like=like, like_ind=like_ind)
    return(res)
}
