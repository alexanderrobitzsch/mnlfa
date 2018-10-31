## File Name: summary.mnlfa.R
## File Version: 0.08


summary.mnlfa <- function( object, file=NULL, ... )
{
    CDM::osink( file=file, suffix=paste0( "__SUMMARY.Rout") )

    cat("-----------------------------------------------------------------------------\n")

    #-- print package
    # mnlfa_print_summary_package(pack="mnlfa")
    cat("\n")

    #-- print call
    CDM::cdm_print_summary_call(object=object)

    #-- print computation time
    CDM::cdm_print_summary_computation_time(object=object)

    cat("Moderated factor analysis\n")

    cat("\n-----------------------------------------------------------------------------\n")
    cat( "Number of iterations=", object$iter, "\n" )
    if ( ! object$converged ){ cat("Maximum number of iterations was reached.\n") }

    cat( "\nDeviance","=", round( object$deviance, 2 ), " | " )
    cat( "Log Likelihood","=", round( -object$deviance/2, 2 ), "\n" )
    cat( "Penalty","=", round( object$regular_penalty, 2 ), "\n" )

    cat( "Number of persons","=", object$ic$n, "\n" )

    cat( "Number of estimated parameters","=", object$ic$np, "\n" )
    cat( "Number of regularized parameters","=", object$ic$numb_reg_pars, "\n\n" )

    #-- information criteria
    CDM::cdm_print_summary_information_criteria(object=object)

    cat("-----------------------------------------------------------------------------\n")
    cat("Item Parameters \n")
    obji <- object$item
    CDM::cdm_print_summary_data_frame(obji, digits=3)

    cat("-----------------------------------------------------------------------------\n")
    cat("Trait Distribution Parameters \n\n")
    obji <- object$parm_trait
    cat("Mean parameters\n")
    print( round(obji$mu, digits=4))
    cat("\nLog standard deviation parameters\n")
    print( round(obji$sigma, digits=4))
        
    CDM::csink( file=file )
}
#*******************************************************
