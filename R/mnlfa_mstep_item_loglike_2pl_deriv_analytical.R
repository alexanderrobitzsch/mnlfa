## File Name: mnlfa_mstep_item_loglike_2pl_deriv_analytical.R
## File Version: 0.046


mnlfa_mstep_item_loglike_2pl_deriv_analytical <- function(parms, a_index, b_index,
        Xdes_int, Xdes_slo, theta, parms_indices, N_item, y, y_resp, post,
        weights=NULL)
{
    b_ii <- parms[ b_index ]
    a_ii <- parms[ a_index ]
    N <- nrow(Xdes_int)
    thetaM <- sirt::sirt_matrix2(theta, nrow=N)
    b <- as.vector(Xdes_int %*% b_ii)
    a <- as.vector(Xdes_slo %*% a_ii)
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
    D2_max <- mnlfa_D2max(D2=D2)


    #--- output
    res <- list(D1=D1, D2=D2, D2_max=D2_max, recycle=NULL)
    return(res)
}
