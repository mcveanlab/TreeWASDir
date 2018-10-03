#' Likelihood surfaces for 400 SNPs in the HES analysis of the UK Biobank
#' 
load.HES.lk.data <- function( ) {

    file = system.file( 'HES_LLK.rdata', package = 'TreeWASDir' )
    load( file )
    out <- list( d = d, res = res )
    return(out)
}
    

#' List containing parameters
#' 
load.pars <- function( ) {

    file = system.file( 'HES_pars.rdata', package = 'TreeWASDir' )
    load( file )
    out <- pars
    return(out)
}
    



