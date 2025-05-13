## File Name: mnlfa_D2max.R
## File Version: 0.01

mnlfa_D2max <- function(D2, eps=1e-3)
{
    res <- max(abs(D2))*(1+eps)
    return(res)
}
