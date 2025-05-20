## File Name: mnlfa_postproc_item.R
## File Version: 0.089


mnlfa_postproc_item <- function(parm_list, items, item_type, parms_estimated,
    parms_regularized, parms_regular_types, parms_regular_lam, parms_regular_alpha)
{
    item <- NULL
    I <- length(items)

    if (sum(unlist(parms_regular_alpha))>0){
        is_alpha <- TRUE
    } else {
        is_alpha <- FALSE
    }

    for (ii in 1L:I){
        parm_ii <- parm_list[[ii]]
        parm_ii <- c( parm_ii$b, parm_ii$a, parm_ii$psi )
        item_ii <- paste(items[ii])
        np_ii <- length(parm_ii)
        sel <- 1L:np_ii
        dfr <- data.frame( itemid=ii, intid=sel, item=rep(item_ii,np_ii))
        dfr$item_type <- item_type[ii]
        dfr$parm <- names(parm_ii)
        dfr$est <- as.vector(parm_ii)
        dfr$estim <- parms_estimated[[ii]][sel]
        dfr$regul <- parms_regularized[[ii]][sel]
        dfr$reg_type <- parms_regular_types[[ii]][sel]
        dfr$reg_lam <- parms_regular_lam[[ii]][sel]
        if (is_alpha){
            dfr$reg_alpha <- parms_regular_alpha[[ii]][sel]
        }
        rownames(dfr) <- NULL
        item <- rbind( item, dfr )
    }
    return(item)
}
