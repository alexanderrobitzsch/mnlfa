## File Name: mnlfa_assign_values_vector.R
## File Version: 0.02

mnlfa_assign_values_vector <- function(x, i, val)
{
    z <- x[i]
    if (length(x)<i){
        z <- val
    }
    return(z)
}
