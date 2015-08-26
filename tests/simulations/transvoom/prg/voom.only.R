voom.only = function (y, design = NULL, col.loc=NULL, span = 0.5, iter=30L, ...) 
{
    if(!is.matrix(y)) y=as.matrix(y)
	if(is.null(col.loc)) col.loc = apply(y, 2L, median)
	col.loc0 = col.loc; 
	col.loc=scale(col.loc, center=TRUE, scale=FALSE)
	
    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }

	y=y-rep(col.loc, each=NROW(y))
	
    fit <- lmFit(y, design, ...)
    if (is.null(fit$Amean)) 
        fit$Amean <- rowMeans(y, na.rm = TRUE)
    sx <- fit$Amean # + mean(log2(lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma)

    l <- lowess(sx, sy, f = span, iter=iter)
    f <- approxfun(l, rule = 2)
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, 
            j, drop = FALSE])
    }
    else {
        fitted.values <- fit$coef %*% t(fit$design)
    }
	
	adj.fit = fitted.values + rep(col.loc, each=NROW(y))
    w <- 1/f(adj.fit)^4
	structure(w, dim=dim(y), design=design, col.loc=col.loc)
}
