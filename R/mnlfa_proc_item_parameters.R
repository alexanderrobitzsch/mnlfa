## File Name: mnlfa_proc_item_parameters.R
## File Version: 0.479


mnlfa_proc_item_parameters <- function(dat, formula_int, formula_slo, formula_res,
    item_type, parms_regular_types, parms_regular_lam, parms_regular_alpha,
    regular_type, regular_lam, regular_alpha, parms_iterations,
    parm_list_init=NULL, K=NULL)
{
    I <- length(item_type)
    items <- names(item_type)

    #* lists of formulas
    formula_int <- mnlfa_proc_item_parameters_formula_list(
                            formula_parm=formula_int, items=items)
    formula_slo <- mnlfa_proc_item_parameters_formula_list(
                            formula_parm=formula_slo, items=items)
    formula_res <- mnlfa_proc_item_parameters_formula_list(
                            formula_parm=formula_res, items=items)

    #** parameter lists
    parm_list <- list()
    parm_Xdes <- list()
    parm_index <- list()
    create_type <- FALSE
    if (is.null(parms_regular_types)){
        parms_regular_types    <- list()
        create_type <- TRUE
    }
    create_lam <- FALSE
    if (is.null(parms_regular_lam)){
        parms_regular_lam <- list()
        create_lam <- TRUE
    }
    create_alpha <- FALSE
    if (is.null(parms_regular_alpha)){
        parms_regular_alpha <- list()
        create_alpha <- TRUE
    }
    create_iter <- FALSE
    if (is.null(parms_iterations)){
        parms_iterations <- list()
        create_iter <- TRUE
    }

    for (ii in 1L:I){

        item_type_ii <- item_type[ii]
        item_ii <- names(item_type)[ii]

        Xdes_int <- stats::model.matrix( object=formula_int[[ii]], data=dat)
        Xdes_slo <- stats::model.matrix( object=formula_slo[[ii]], data=dat)

        if (item_type_ii %in% c('NO') ){
            Xdes_res <- stats::model.matrix( object=formula_res[[ii]], data=dat)
        } else {
            Xdes_res <- NULL
        }
        des1 <- list(Xdes_int=Xdes_int, Xdes_slo=Xdes_slo, Xdes_res=Xdes_res)
        attr(des1, 'item') <- item_ii
        attr(des1, 'item_type') <- item_type_ii


        #* 2PL or 1PL
        if (item_type_ii %in% c('1PL','2PL','NO','GPCM') ){
            # design matrices

            v1 <- list()
            #- item intercepts
            res <- mnlfa_create_inits_item_parameters(Xdes=Xdes_int, item_ii=item_ii,
                            formula=formula_int[[ii]], des1=des1, ref_value=0,
                            par='b', parlab='int', K=K[ii])
            v1[['b']] <- res$g1
            des1 <- res$des1
            offset_int <- res$offset

            #- item slopes
            res <- mnlfa_create_inits_item_parameters(Xdes=Xdes_slo, item_ii=item_ii,
                            formula=formula_slo[[ii]], des1=des1, ref_value=1,
                            par='a', parlab='slo')
            v1[['a']] <- res$g1
            des1 <- res$des1
            offset_slo <- res$offset

            #- item residual SDs
            if (item_type_ii=='NO'){
                res <- mnlfa_create_inits_item_parameters(Xdes=Xdes_res, item_ii=item_ii,
                                formula=formula_res[[ii]], des1=des1, ref_value=0,
                                par='psi', parlab='res')
                v1[['psi']] <- res$g1
                des1 <- res$des1
                offset_res <- res$offset
            } else {
                v1$psi <- NULL
            }

            if (item_type_ii %in% c('GPCM')){
                u1 <- K[ii]
                names(u1) <- paste0(item_ii, '_K')
                v1[['K']] <- u1
            }

            #- regularization parameters settings
            res <- mnlfa_proc_item_parameters_parameter_settings(v1=v1,
                            regular_type=regular_type, regular_lam=regular_lam,
                            regular_alpha=regular_alpha, item_type=item_type_ii)
            h1 <- res$h1
            reg1 <- res$reg1
            reg2 <- res$reg2
            reg3 <- res$reg3
            np_a <- res$np_a
            np_b <- res$np_b
            np_psi <- res$np_psi

        }  # if item_type_ii=...

        parm_list[[ii]] <- v1
        parm_Xdes[[ii]] <- des1
        parm_index[[ii]] <- h1
        if (create_type){
            parms_regular_types[[ii]] <- reg1
        }
        if (create_lam){
            parms_regular_lam[[ii]] <- reg2
        }
        if (create_alpha){
            parms_regular_alpha[[ii]] <- reg3
        }

        if (create_iter){
            if (item_type[[ii]] %in% c('2PL','NO','GPCM')){
                parms_iterations[[ii]] <- as.list(1L:(np_b+np_a+np_psi))
                if (offset_int & ( ! offset_slo) ){
                    parms_iterations[[ii]] <- as.list(2L:(np_b+np_a+np_psi))
                }
                if ( (!offset_int) & ( offset_slo)){
                    parms_iterations[[ii]] <- as.list(1L:(np_b+np_psi))
                }
                if (offset_int & offset_slo){
                    parms_iterations[[ii]] <- as.list(NULL)
                }
            }
            if (item_type[[ii]]=='1PL'){
                parms_iterations[[ii]] <- as.list(1L:(np_b))
                if (offset_int){
                    parms_iterations[[ii]] <- as.list(NULL)
                }
            }
        } #- end create iter

    }  ## end item ii

    if (!is.null(parm_list_init)){
        parm_list <- parm_list_init
    }

    names(parm_list) <- items
    names(parm_Xdes) <- items
    names(parms_regular_types) <- items
    names(parms_regular_lam) <- items
    names(parms_iterations) <- items
    names(parm_index) <- items

    #--- output
    res <- list( parm_list=parm_list, parm_Xdes=parm_Xdes,
                    parms_regular_types=parms_regular_types,
                    parms_regular_lam=parms_regular_lam,
                    parms_regular_alpha=parms_regular_alpha,
                    parms_iterations=parms_iterations,
                    parm_index=parm_index)
    return(res)
}


