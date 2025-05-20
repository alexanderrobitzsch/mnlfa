## File Name: mnlfa_mstep_item_loglike_deriv_analytical.R
## File Version: 0.178


mnlfa_mstep_item_loglike_deriv_analytical <- function(parms, a_index, b_index, K_index,
        psi_index, Xdes_int, Xdes_slo, Xdes_res, theta, parms_indices, N_item, y, y_resp,
        post, item_type, weights=NULL, use_rcpp=TRUE)
{
    b_ii <- parms[ b_index ]
    a_ii <- parms[ a_index ]
    psi_ii <- parms[ psi_index ]
    K_ii <- parms[ K_index ]

    N <- nrow(Xdes_int)
    thetaM <- sirt::sirt_matrix2(theta, nrow=N)

    if (item_type %in% c('1PL','2PL','NO')){
        b <- Xdes_int %*% b_ii
    }

    if (item_type %in% c('GPCM')){
        K_ii <- parms[ K_index ]
        parms_ii <- list(b=b_ii, K=K_ii)
        b <- mnlfa_compute_moderated_parameter_gpcm_b(Xdes=Xdes_int,
                    parm=parms_ii, value=0)
    }

    a <- as.vector(Xdes_slo %*% a_ii)

    #--- 1PL/2PL
    if (item_type %in% c('1PL','2PL') ){
        U <- stats::plogis(a*thetaM - b )
        if ( parms_indices %in% b_index ){
            ind <- match(parms_indices, b_index)
            fac <- -Xdes_int[, ind]
        }
        if ( parms_indices %in% a_index ){
            ind <- match(parms_indices, a_index)
            fac <- Xdes_slo[, ind]*thetaM
        }
        y_resp_post <- weights*y_resp*post
        D1 <- sum( ( y - U ) * y_resp_post*fac )/N_item
        #- second-order derivative
        D2 <- - sum( U * (1-U) * fac^2 * y_resp_post) / N_item
    }

    #--- normal distribution
    if (item_type %in% c('NO') ){
        U <- b + a*thetaM
        psi <- exp( as.vector(Xdes_res %*% psi_ii) )
        meanstr <- TRUE
        if ( parms_indices %in% b_index ){
            ind <- match(parms_indices, b_index)
            fac <- Xdes_int[, ind]
        }
        if ( parms_indices %in% a_index ){
            ind <- match(parms_indices, a_index)
            fac <- Xdes_slo[, ind]*thetaM
        }
        if ( parms_indices %in% psi_index ){
            ind <- match(parms_indices, psi_index)
            fac <- Xdes_res[, ind]
            meanstr <- FALSE
        }
        if (meanstr){
            H1 <- weights*post*(y-U)/psi^2
        } else {
            Z <- (y-U)^2 / psi^2
            H1 <- weights*post*( Z - 1)
        }
        D1 <- sum( fac*H1*y_resp )/N_item

        if (meanstr){
            H2 <- -weights*post*1/psi^2*fac^2
        } else {
            H2 <- -2*weights*post*fac^2*Z
        }
        D2 <- sum(H2*y_resp)/N_item

    }

    #--- GPCM
    if (item_type %in% c('GPCM') ){

        is_b <- FALSE
        if ( parms_indices %in% b_index ){
            ind <- match(parms_indices, b_index)
            is_b <- TRUE
        }
        if ( parms_indices %in% a_index ){
            ind <- match(parms_indices, a_index)
        }

        y_resp_post <- y_resp * post * weights
        uu <- parms_indices

        args <- list( b=b, a=a, theta=theta, Xdes_int=Xdes_int, Xdes_slo=Xdes_slo, y=y,
                        y_resp=y_resp, ind=ind, is_b=is_b, uu=parms_indices )

        if (use_rcpp){
            fun <- 'mnlfa_rcpp_calc_probs_gpcm_deriv'
        } else {
            fun <- 'mnlfa_mstep_item_loglike_deriv_analytical_gpcm'
        }
        res <- do.call(what=fun, args=args)
        H1 <- res$H1
        H2 <- res$H2

        D1 <- sum( y_resp_post*H1 )/N_item
        D2 <- sum( y_resp_post*H2 )/N_item

    }

    D2_max <- mnlfa_D2max(D2=D2)

    #--- output
    res <- list(D1=D1, D2=D2, D2_max=D2_max, recycle=NULL)
    return(res)
}

