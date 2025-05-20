## File Name: mnlfa_create_inits_item_parameters.R
## File Version: 0.102

mnlfa_create_inits_item_parameters <- function(Xdes, item_ii, des1, formula, ref_value=0,
        par="b", parlab="int", K=NULL)
{
    nc <- ncol(Xdes)
    offset <- FALSE
    if (nrow(Xdes)==1){
        nc <- 0
    }
    if (nc > 0){
        g1 <- rep(0,ncol(Xdes) )
        g1[1] <- ref_value
        names(g1) <- paste0(item_ii, '_', par, '_', colnames(Xdes) )

        item_type <- attr(des1, 'item_type')
        if ( ( item_type %in% c('GPCM') ) & ( par=='b')  ){
            P2 <- ncol(Xdes)+K-1
            g1 <- rep(0,P2)
            g1[1L:K] <- 1.5*seq(-1,1,length=K)
            names(g1)[1L:K] <- paste0(item_ii, '_', par, '_Cat', 1L:K )
            if (P2>K){
                names(g1)[(K+1):P2] <- paste0(item_ii, '_', par, '_', colnames(Xdes)[-1])
            }
        }

    } else {
        g1 <- c(ref_value)
        # off <- stats::model.offset(x=formula)
        names(g1) <- paste0(item_ii, '_', par, '_offset' )
        Xdes_label <- paste0('Xdes_', parlab)
        des1[[ Xdes_label ]] <- matrix(0, nrow=nrow(Xdes), ncol=1)
        colnames(des1[[ Xdes_label ]] ) <- paste0(item_ii, '_', par, '_offset' )
        offset <- TRUE
    }
    #-- result list
    res <- list(g1=g1, offset=offset, des1=des1)
    return(res)
}
