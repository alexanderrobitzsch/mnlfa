## File Name: mnlfa_proc_item_parameters_parameter_settings.R
## File Version: 0.046

mnlfa_proc_item_parameters_parameter_settings <- function(v1, regular_type, regular_lam,
        regular_alpha, item_type)
{

    np_b <- length(v1$b)
    np_a <- length(v1$a)
    np_psi <- length(v1$psi)
    NP <- np_b + np_a + np_psi
    reg1 <- rep('none', NP)
    reg3 <- reg2 <- rep(0, NP)
    names(reg1) <- c( names(v1$b), names(v1$a), names(v1$psi) )
    names(reg3) <- names(reg2) <- names(reg1)

    NH <- 2
    if (item_type=='NO'){
        NH <- 3
    }

    for (hh in 1L:NH){
        if (hh==1){
            ind <- 1 + seq_len( np_b - 1 )
            if (item_type=='GPCM'){
                K <- as.numeric(v1['K'])
                if (np_b > K){
                    ind <- K + seq_len( np_b - K)
                } else {
                    ind <- NULL
                }
            }
        }
        if (hh==2){
            ind <- np_b + 1 + seq_len( np_a - 1 )
        }
        if (hh==3){
            ind <- np_b + np_a + 1 + seq_len( np_psi - 1 )
        }
        reg1[ind] <- mnlfa_assign_values_vector(x=regular_type, i=hh, val='none')
        reg2[ind] <- mnlfa_assign_values_vector(x=regular_lam, i=hh, val=0)
        reg3[ind] <- mnlfa_assign_values_vector(x=regular_alpha, i=hh, val=0)
    }

    #- parameter indices
    b_index <- seq_len(np_b)
    # if (offset_int){ b_index <- NULL; }
    a_index <- np_b + seq_len(np_a)
    psi_index <- np_b + np_a + seq_len(np_psi)
    K_index <- np_b + np_a + np_psi + 1

    h1 <- list(b_index=b_index, a_index=a_index, psi_index=psi_index,
                    K_index=K_index)

    #-- output
    res <- list(h1=h1, reg1=reg1, reg2=reg2, reg3=reg3, np_a=np_a, np_b=np_b,
                    np_psi=np_psi)
    return(res)
}
