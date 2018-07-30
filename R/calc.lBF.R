#' calculate BF
#'
#' @param pars A list object as returned by the function set.parameters
#' @param data.sub A 3-dimension array (N. SNPs x N. phenotypes x 3) with state likelihoods.
#' @param w0 An integer with column index for the null state.
#' @param log10 Return Bayes Factor in log 10 base.
#'
#' @return lBF Tree Bayes Factor
#'
calc.lBF <- function(
    pars = NULL,
    data.sub = NULL,
    w0 = 2,
    log10 = TRUE
) {

    if( is.null(pars) | is.null(data.sub) ) {
        stop("Missing input data.\n")
    }
    
    if (ncol(data.sub)==2) w0 <- 1;

    llk.full <- calc.llk.tree(pars, data.sub);

    q00 <- pars$p.stay + pars$p.switch * (1-pars$pi1);
    lq00 <- log(q00);
    n.trans <- nrow(data.sub)-1;
    l.p.null <- log(1-pars$pi1)+n.trans*lq00;
    l.lk.null <- sum(data.sub[,w0]);
    tmp <- c(llk.full, l.p.null+l.lk.null);
    mx <- max(tmp);
    tmp <- exp(tmp-mx);
    lBF <- mx+log(tmp[1]-tmp[2])-l.lk.null-log(1-exp(l.p.null));
    if (log10) lBF <- lBF/log(10);
    
    return(lBF);
}

