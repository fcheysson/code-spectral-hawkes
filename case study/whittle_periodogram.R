whittle_periodogram = function (periodogram, kern, binsize = NULL, trunc = 5L, init, 
          ...) 
{
    if (is.character(kern)) {
        kern = match.arg(tolower(kern), c("exponential", 
                                          "symmetricexponential", "gaussian", "powerlaw", 
                                          "pareto3", "pareto2", "pareto1"))
        switch(kern, exponential = {
            kern = new(Exponential)
        }, symmetricexponential = {
            kern = new(SymmetricExponential)
        }, gaussian = {
            kern = new(Gaussian)
        }, powerlaw = {
            kern = new(PowerLaw)
        }, pareto3 = {
            kern = new(Pareto3)
        }, pareto2 = {
            kern = new(Pareto2)
        }, pareto1 = {
            kern = new(Pareto1)
        })
    }
    else if (!any(sapply(paste0("Rcpp_", c("Exponential", 
                                           "SymmetricExponential", "Gaussian", "PowerLaw", 
                                           "Pareto3", "Pareto2", "Pareto1")), 
                         function(class_) {
                             is(kern, class_)
                         }))) 
        stop("'kern' must be a valid kernel.")
    if (!is.null(binsize)) 
        kern$binsize = binsize
    n <- length(periodogram)
    I = periodogram
    nlopt_fn <- function(param) {
        kern$param <- param
        return(kern$whittle(I, trunc))
    }
    x0 = init
    optargs = list(hessian = TRUE, lower = rep(0.01, length(kern$param)), 
                   upper = c(Inf, 0.99, rep(Inf, length(kern$param) - 2)), 
                   method = "L-BFGS-B")
    optargs = modifyList(optargs, list(...))
    opt <- do.call(optim, c(list(par = x0, fn = nlopt_fn), optargs))
    kern$param = opt$par
    output = kern
    return(output)
}