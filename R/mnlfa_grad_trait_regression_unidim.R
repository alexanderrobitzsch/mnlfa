## File Name: mnlfa_grad_trait_regression_unidim.R
## File Version: 0.063


mnlfa_grad_trait_regression_unidim <- function(x, theta, post, Xdes_mean, Xdes_sd,
        index_mu, index_sigma, weights=NULL, use_rcpp=TRUE)
{
    np <- length(x)
    grad <- rep(NA, np)
    N <- nrow(Xdes_mean)
    thetaM <- sirt::sirt_matrix2(theta, nrow=N)

    mu <- x[ index_mu ]
    sigma <- x[ index_sigma ]
    mu_p <- as.vector( Xdes_mean %*% mu )
    sigma_p <- as.vector( exp( Xdes_sd %*% sigma ) )

    IM <- length(index_mu)
    IS <- length(index_sigma)

    if (!use_rcpp){

        ## derivative with respect to mu is (x-mu)/sigma^2
        H1 <- weights*post * (thetaM - mu_p) / sigma_p^2
        for (hh in 1L:IM){
            grad[hh] <- - sum( H1 * Xdes_mean[,hh] )
        }

        ## derivative with respect to sigma is -1/sigma + (x-mu)^2/sigma^3
        ## multiply it with inner derivative sigma
        H2 <- weights*post * ( - 1 + (thetaM - mu_p)^2 / sigma_p^2    )
        for (hh in 1L:IS){
            grad[IM+hh] <- - sum( H2  * Xdes_sd[,hh] )
        }

    } else {

        grad <- mnlfa_rcpp_mstep_trait_unidim_grad( theta=theta,
                    mu_p=mu_p, sigma_p=sigma_p, post=post, IM=IM, IS=IS,
                    Xdes_mean=Xdes_mean, Xdes_sd=Xdes_sd, weights=weights)
    }

    return(grad)
}
