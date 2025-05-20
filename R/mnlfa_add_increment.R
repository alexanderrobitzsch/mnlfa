## File Name: mnlfa_add_increment.R
## File Version: 0.01

mnlfa_add_increment <- function(x, i, h)
{
    x[i] <- x[i] + h
    return(x)
}
