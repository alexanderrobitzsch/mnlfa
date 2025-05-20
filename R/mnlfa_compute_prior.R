## File Name: mnlfa_compute_prior.R
## File Version: 0.03

mnlfa_compute_prior <- function(parm_trait, Xdes_mean, Xdes_sd, N, theta)
{
    mu_p <- mnlfa_predicted_values(X=Xdes_mean, parm=parm_trait$mu, type='ident')
    sigma_p <- mnlfa_predicted_values(X=Xdes_sd, parm=parm_trait$sigma, type='exp')
    TP <- nrow(theta)
    prior <- matrix(1, nrow=N, ncol=TP)
    for (tt in 1L:TP){
        prior[,tt] <- stats::dnorm(theta[rep(tt,N),1], mean=mu_p, sd=sigma_p )
    }
    prior <- prior/rowSums(prior)
    return(prior)
}
