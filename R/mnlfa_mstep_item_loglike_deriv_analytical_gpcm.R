## File Name: mnlfa_mstep_item_loglike_deriv_analytical_gpcm.R
## File Version: 0.05


mnlfa_mstep_item_loglike_deriv_analytical_gpcm <- function(b, a, theta, Xdes_int,
            Xdes_slo, y, y_resp, ind, is_b, uu)
{
    K <- ncol(b)
    N <- nrow(Xdes_int)

    KM <- sirt::sirt_matrix2(1L:K, nrow=N)
    KM2 <- sirt::sirt_matrix2(0:K, nrow=N)
    KM22 <- sirt::sirt_matrix2((0:K)^2, nrow=N)
    ind1 <- cbind( 1L:N, y+1 )
    TP <- nrow(theta)
    H2 <- H1 <- matrix(0, nrow=N, ncol=TP)

    for (tt in 1L:TP){
        M1 <- a*theta[tt,1]*KM - b
        M1 <- cbind(1, exp( M1 ) )
        fj <- M1 / rowSums( M1 )

        if (ind>K | (!is_b) ){
            Efj <- rowSums( KM2*fj )
            Vfj <- rowSums( KM22*fj ) - Efj^2
        }

        if (ind<=K & is_b){
            g1 <- fj[,uu+1]
            g2 <- ifelse(y==uu, g1-1, g1)
            h2 <- -g1*(1-g1)
        }
        if (ind>K & is_b){
            z <- Xdes_int[,ind-K+1]
            g2 <- -z*(y-Efj)
            h2 <- -z^2*Vfj
        }
        if (!is_b){
            z <- Xdes_slo[,ind]
            g2 <- theta[tt,1]*z*(y-Efj)
            h2 <- -theta[tt,1]^2*z^2*Vfj
        }
        H1[,tt] <- g2
        H2[,tt] <- h2
    }

    #-- output
    res <- list(H1=H1, H2=H2)
    return(res)

}
