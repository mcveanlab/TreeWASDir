#' Function to calculate integrated likelihood for set of variants up tree
#'
#' @param pars A list object as returned by the function set.parameters
#' @param data.sub A 3-dimension array (N. SNPs x N. phenotypes x 3) with state likelihoods.
#' @param returnForwardMatrix Boolean
#' @param returnGMatrix Boolean
#'
#' @return Either a numeric object with the Tree likelihood or a list with FWD and G matrices

calc.llk.tree <- function(
    pars,
    data.sub,
    returnForwardMatrix = FALSE,
    returnGMatrix = FALSE
) {

    ## Get integrated likelihood at nodes - will be overwritten for internal nodes
    mx <- apply(data.sub, 1, max);
    d <- exp(data.sub-mx);
    llk.integrated <- log(d %*% pars$stat.dist)+mx;
    if (returnGMatrix) g <- array(0, dim(d));
    
    for (i in pars$t.path) {
        ## Emissions at node
        emiss.node<-data.sub[i,];			
        data.sub[i,]<-0;
        for (j in pars$ontology[[i]]) {
            tmp1 <- cbind(
                data.sub[j,]-pars$theta.tree,
                llk.integrated[j]+pars$lp.switch
            );
            mx1 <- apply(tmp1, 1, max);
            tmp2 <- mx1+log(rowSums(exp(tmp1-mx1)));
            data.sub[i,] <- data.sub[i,]+tmp2;
            if (returnGMatrix) g[j,] <- tmp2;
        }
        data.sub[i,] <- data.sub[i,]+emiss.node;
        mx <- max(data.sub[i,]);
        llk.integrated[i] <- mx+log(sum(exp(data.sub[i,]-mx)*pars$stat.dist));
    }
    if (returnForwardMatrix) {
        if (!returnGMatrix) {
            return(data.sub);
        } else {
            return(list(f=data.sub, g=g));
        }
    } else {
        return(llk.integrated[i]);
    }
}
