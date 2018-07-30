#' Function to get marginal posterior on -/0/+ profile
#'
#' @param pars A list object as returned by the function set.parameters
#' @param data.sub A 3-dimension array (N. SNPs x N. phenotypes x 3) with state likelihoods.
#' @return pp. A 2-dimension array (N. Phenotypes x 3) with posterior probabilties.
#' 

marginal.posterior.profile <- function(
    pars = NULL,
    data.sub = NULL
) {

    ## Get forward and G matrices
    tmp <- calc.llk.tree(pars, data.sub, returnForwardMatrix=T, returnGMatrix=T);
    f <- tmp$f;
    g <- tmp$g;
    
    ## Build parents list and reverse F traversal order
    parents <- rep(0, pars$n.phenos);
    for (i in pars$t.path) parents[pars$ontology[[i]]] <- i;
    ord <- c(rev(pars$t.path), pars$terminals);

    ## Construct backward matrix
    b <- array(0, dim(f));
    b[ord[1],] <- log(pars$stat.dist);
    for (i in ord[-1]) {
        r.i <- b[parents[i],]+f[parents[i],]-g[i,];
        mx <- max(r.i);
        tmp <- mx+log(sum(exp(r.i-mx)));
        tmp <- tmp+pars$lp.switch+log(pars$stat.dist);
        tmp2 <- -pars$theta.tree+b[parents[i],]+f[parents[i],]-g[i,];
        tmp3 <- cbind(tmp2, tmp);
        mx <- apply(tmp3, 1, max);
        b[i,] <- log(rowSums(exp(tmp3-mx)))+mx;
    }
    
    ## Posteriors
    tmp <- b+f;
    mx <- apply(tmp, 1, max);
    pp <- exp(tmp-mx);
    pp <- pp/rowSums(pp);
    return(pp);
}
