## File Name: mnlfa_compute_moderated_parameter.R
## File Version: 0.055

mnlfa_compute_moderated_parameter <- function(Xdes, parm, value=0, trafo="none")
{
    if (ncol(Xdes)>0){
        res <- Xdes %*% parm
    } else {
        res <- matrix(value, nrow=nrow(Xdes),ncol=1)
    }
    if (trafo=='exp'){
        res <- exp(res)
    }
    return(res)
}
