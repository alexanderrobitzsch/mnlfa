## File Name: mnlfa_sbic_penalty.R
## File Version: 0.02

mnlfa_sbic_penalty <- function(par, eps=0.001)
{
    val <- par^2 / (par^2 + eps)
    return(val)
}
