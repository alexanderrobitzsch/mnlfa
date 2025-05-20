## File Name: mnlfa_predicted_values.R
## File Version: 0.043

mnlfa_predicted_values <- function(X, parm, type="ident")
{
    if (ncol(X)>0){
        y <- X %*% parm
    } else {
        y <- 0
    }
    if (type=='exp'){
        y <- exp(y)
    }
    return(as.vector(y))
}
