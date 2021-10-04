# !diagnostics off

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(ggplot2)
library(hawkesbow)
library(latex2exp)
library(foreach)
library(doParallel)

numCores = 3


# Exponential kernel ------------------------------------------------------


n.sim = 1e3

burn.in = 100
Ts = c(50, 100, 150, 250, 400, 650, 1000, 2000, 4000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
rate = 1

cl = makePSOCKcluster(numCores) # use multicore, set to the number of our cores
clusterExport(cl, c("burn.in", "Ts", "binsizes", "lambda", "mu", "rate"))
registerDoParallel(cl)

estimates = foreach (sim = icount(n.sim), .combine = rbind, .multicombine = TRUE, .inorder = FALSE,
                     .packages = c("hawkesbow"), .verbose = FALSE) %dopar% {
                         
                         # # Counter display
                         # if (sim %% 10 == 0) cat(sim, ". ", sep="")
                         
                         estimates = data.frame(T = numeric(length(Ts) * (1+length(binsizes))),
                                                simulation = integer(length(Ts) * (1+length(binsizes))),
                                                binsize = character(length(Ts) * (1+length(binsizes))),
                                                baseline = numeric(length(Ts) * (1+length(binsizes))),
                                                reprod = numeric(length(Ts) * (1+length(binsizes))),
                                                rate = numeric(length(Ts) * (1+length(binsizes))),
                                                message = character(length(Ts) * (1+length(binsizes))))
                         
                         model = new(Exponential)
                         it = 1
                         
                         for (Tend in Ts) {
                             
                             # Simulate Hawkes process
                             x = hawkes(burn.in + Tend, fun = lambda, repr = mu, family = "exp", rate = rate)
                             x$p = x$p[x$p > burn.in] - burn.in
                             x$end = Tend
                             
                             # MLE
                             opt = mle(x$p, "Exponential", Tend, init = c(lambda, mu, rate), lb = rep(0.05, 3), ub = c(50, .95, 50))
                             
                             estimates$T[it] = Tend
                             estimates$simulation[it] = sim
                             estimates$binsize[it] = "MLE"
                             estimates$baseline[it] = opt$par[1]
                             estimates$reprod[it] = opt$par[2]
                             estimates$rate[it] = opt$par[3]
                             estimates$message[it] = opt$opt$message
                             
                             it = it + 1
                             
                             for (binsize in binsizes) {
                                 
                                 # Whittle
                                 y = discrete(x, binsize = binsize)
                                 opt = whittle(y, model, binsize, lower = rep(0.05, 3), upper = c(50, .95, 50), control = list(ndeps = rep(1e-5, 3)))
                                 
                                 estimates$T[it] = Tend
                                 estimates$simulation[it] = sim
                                 estimates$binsize[it] = as.character(binsize)
                                 estimates$baseline[it] = opt$par[1]
                                 estimates$reprod[it] = opt$par[2]
                                 estimates$rate[it] = opt$par[3]
                                 estimates$message[it] = opt$opt$message
                                 
                                 it = it + 1
                                 
                             }
                             
                         }
                         
                         return(estimates)
                         
                     } %>% as_tibble() %>% mutate(message = as.factor(message))

stopCluster(cl)

save(estimates, file = "exponential_convergence2.RData")

estimates.long = estimates %>% mutate(T = as.factor(T)) %>%
    gather("parameter", "estimate", baseline:rate) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "rate")))

trueval = tibble(T = rep(Ts, each = (1+length(binsizes)) * 3),
                 binsize = rep(c("MLE", as.character(binsizes)), length(Ts) * 3),
                 parameter = rep(rep(c("baseline", "reprod", "rate"), each = 1 + length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, rate), each = 1 + length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "rate")))

ggplot(estimates.long %>% 
           filter(T %in% paste0("T = ", c(100, 1000))), 
       aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=filter(trueval, T %in% paste0("T = ", c(100, 1000))), size=4, colour = "firebrick") +
    scale_x_discrete(labels=c("baseline"=expression(hat(eta)), "reprod"=expression(hat(mu)), "rate"=expression(hat(beta)))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,4)) +
    xlab("Parameters") +
    ylab("Estimates") +
    # labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = 1 e^{-1t}$ on $(0, T)$ | true values are larger points")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("exponential_estimates.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("exponential_estimates.png", width=13.6, height=7.6, units="in")

estimates.variance = estimates %>%
    gather("parameter", "estimate", baseline:rate) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("bin size = ", binsizes))),
                         levels = c("Regular MLE", paste0("bin size = ", binsizes))),
        parameter = factor(parameter, levels = c("baseline", "reprod", "rate")),
        trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, rate))
    ) %>%
    group_by(T, binsize, parameter) %>%
    summarise(bias.squared = mean(estimate - trueval)^2, 
              variance = var(estimate), 
              MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>%
    filter(T == 100) %>%
    group_by(parameter) %>%
    summarise(min = min(MSE), max = max(MSE))

ggplot(estimates.variance %>% filter(T != 50),
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point(size = 2) +
    geom_abline(intercept = -.1, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .0025, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    annotate(geom = "label",
             x = 90,
             y = sqrt(estimates.T100$min * estimates.T100$max),
             label = c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta)),
             hjust = 1, size = 7
    ) +
    annotate(geom = "linerange", x = sqrt(90*100), ymin = estimates.T100$min * .925, ymax = estimates.T100$max * 1.075,
             size = 2, color = c("grey80", "grey50", "grey50")) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Bin size",
                       labels = c("MLE", "0.25", "0.5", "1", "2"),
                       values = c("firebrick", viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta))) +
    xlab("T") +
    ylab("Mean Square Error") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("exponential_convergence.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("exponential_convergence.png", width=13.6, height=7.6, units="in")

estimates.standard = estimates.long %>%
    inner_join(trueval %>% rename(trueval = estimate)) %>%
    inner_join(estimates.long %>% group_by(T, binsize, parameter) %>% summarise(stdev = sd(estimate))) %>%
    mutate(stdEstimate = (estimate - trueval) / stdev) %>%
    select(- estimate, - trueval, - stdev)

ggplot(estimates.standard %>% filter(T == "T = 1000"), aes(sample = stdEstimate)) +
    qqplotr::stat_qq_band() +
    qqplotr::stat_qq_point() +
    qqplotr::stat_qq_line(col = "red") +
    facet_grid(parameter ~ binsize) +
    theme_bw()

estimates.T100 = estimates.variance %>%
    filter(T == 100) %>%
    group_by(parameter) %>%
    summarise(min = min(bias.squared), max = max(bias.squared))

estimates.at100 = estimates.variance %>% 
    filter(T == 100)

ggplot(estimates.variance %>% filter(T != 50),
       aes(x=T, y=bias.squared, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point() +
    geom_abline(intercept = -.1, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .0025, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    annotate(geom = "segment",
             x = as.numeric(as.character(
                 plyr::revalue(estimates.at100$parameter,
                               c("baseline"=82.65, "reprod"=90, "rate"=98))
             )),
             xend = 100,
             y = estimates.at100$bias.squared,
             yend = estimates.at100$bias.squared,
    ) +
    annotate(geom = "label",
             x = as.numeric(as.character(
                 plyr::revalue(estimates.at100$parameter,
                               c("baseline"=82.65, "reprod"=90, "rate"=98))
             )),
             y = estimates.at100$bias.squared,
             label = plyr::revalue(
                 estimates.at100$parameter,
                 c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta))
             ),
             hjust = 1, parse = TRUE
    ) +
    # annotate(geom = "label",
    #          x = 90,
    #          y = sqrt(estimates.T100$min * estimates.T100$max),
    #          label = c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta)),
    #          hjust = 1, size = 7
    # ) +
    # annotate(geom = "linerange", x = sqrt(90*100), ymin = estimates.T100$min * .925, ymax = estimates.T100$max * 1.075,
    #          size = 2, color = c("grey80", "grey50", "grey50")) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Bin size",
                       labels = c("MLE", "0.25", "0.5", "1", "2"),
                       values = c("firebrick", viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta))) +
    xlab("T") +
    ylab("Bias Squared") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))


# FOR THE ESTIMATION OF POWER LAW KERNELS WITH a FIXED --------------------


whittle = function(counts, kern, binsize = NULL, trunc = 5L, init = NULL, ...) {
    
    # Check that the argument 'kern' is either a string that matches of the kernels implemented
    if (is.character(kern)) {
        kern = match.arg(tolower(kern),
                         c("exponential", "symmetricexponential", "gaussian",
                           "powerlaw", "pareto3", "pareto2", "pareto1"))
        switch(kern,
               exponential = {kern = new(Exponential)},
               symmetricexponential = {kern = new(SymmetricExponential)},
               gaussian = {kern = new(Gaussian)},
               powerlaw = {kern = new(PowerLaw)},
               pareto3 = {kern = new(Pareto3)},
               pareto2 = {kern = new(Pareto2)},
               pareto1 = {kern = new(Pareto1)})
    } else if ( # or that it refers to a valid hawkes kernel
        !any(sapply(
            paste0("Rcpp_", c("Exponential", "SymmetricExponential", "Gaussian", "PowerLaw", "Pareto3", "Pareto2", "Pareto1")),
            function(class_) {is(kern, class_)}
        ))
    ) stop("'kern' must be a valid kernel.")
    
    if (!is.null(binsize)) kern$binsize = binsize
    
    # Periodogram
    n <- length(counts)
    dft <- fft(counts - mean(counts))
    I <- Mod(dft)^2 / n
    I <- I[-1]      # remove value corresponding to density in 0 (not well defined for centered processes)
    
    # Whittle pseudo likelihood function (for optim)
    nlopt_fn <- function(param) {
        kern$param[1:3] <- param
        return( kern$whittle(I, trunc) )
    }
    
    # Sensible initialisation
    if (is.null(init)) {
        ymean = mean(counts)
        # For PowerLaw
        if (is(kern, "Rcpp_PowerLaw")) {
            wmin = Inf
            mu = .25
            shape = 1
            for (mu_ in c(.25, .5, .75)) {
                for (shape_ in 1:10/2) {
                    kern$param[1] = ymean * (1 - mu_)
                    kern$param[2] = mu_
                    kern$param[3] = shape_
                    whit = kern$whittle(I, trunc = trunc)
                    if (whit < wmin) {
                        mu = mu_
                        shape = shape_
                        wmin = whit
                    }
                }
            }
            x0 = c(ymean * (1 - mu), mu, shape)
        } else if (is(kern, "Rcpp_Exponential")) {
            wmin = Inf
            mu = .25
            rate = 1
            for (mu_ in c(.25, .5, .75)) {
                for (rate_ in 1:5) {
                    kern$param[1] = ymean * (1 - mu_)
                    kern$param[2] = mu_
                    kern$param[3] = rate_
                    whit = kern$whittle(I, trunc = trunc)
                    if (whit < wmin) {
                        mu = mu_
                        rate = rate_
                        wmin = whit
                    }
                }
            }
            x0 = c(ymean * (1 - mu), mu, rate)
        } else {
            x0 = c(ymean * 2,
                   .5,
                   runif(length(kern$param)-2, 0, 10))
        }
    } else x0 = init
    
    optargs = list(hessian = TRUE,
                   lower = rep(.01, length(kern$param) - 1),
                   upper = c(Inf, .99, rep(Inf, length(kern$param) - 3)),
                   method = "L-BFGS-B")
    
    optargs = modifyList(optargs, list(...))
    opt <- do.call(optim, c(list(par = x0, fn = nlopt_fn), optargs))
    
    # Create output object
    kern$param[1:3] = opt$par
    
    output = list(
        par = opt$par,
        kernel = kern,
        counts = counts,
        binsize = kern$binsize,
        opt = opt
    )
    
    class(output) = "HawkesModel"
    
    return( output )
}

mle = function(events, kern, end, init = NULL, opts = NULL, ...) {
    
    # Check that the argument 'kern' is either a string that matches of the kernels implemented
    if (is.character(kern)) {
        kern = match.arg(tolower(kern),c("exponential", "powerlaw"))
        switch(kern,
               exponential = {model = new(Exponential)},
               powerlaw = {model = new(PowerLaw)})
    } else if ( # or that it refers to a valid hawkes model
        !any(sapply(
            paste0("Rcpp_", c("Exponential", "PowerLaw")),
            function(class_) {is(model, class_)}
        ))
    ) stop("'kern' must be a valid kernel.")
    
    # Likelihood function (for nloptr)
    nlopt_fn <- function(param) {
        model$param[1:3] <- param
        loglikngrad = lapply(X = model$loglikngrad(events, end), FUN = "-")
        loglikngrad$gradient = loglikngrad$gradient[1:3]
        return( loglikngrad )
    }
    
    if (is.null(opts))
        opts <- list("algorithm" = "NLOPT_LD_LBFGS")
    else {
        if (is.null(opts[["algorithm"]]))
            opts <- c(opts, "algorithm" = "NLOPT_LD_LBFGS")
        if (is.null(opts[["xtol_rel"]]))
            opts <- c(opts, "xtol_rel" = 1e-04)
    }
    
    if (is.null(init))
        x0 = c(runif(1, 0, 2),
               runif(1, 0, 1),
               runif(length(model$param)-3, 0, 10))
    else x0 = init
    
    optargs = list(lb = rep(.0001, length(model$param) - 1),
                   ub = c(Inf, .9999, rep(Inf, length(model$param) - 3)),
                   opts = opts)
    optargs = modifyList(optargs, list(...))
    opt <- do.call(nloptr::nloptr, c(list(x0 = x0, eval_f = nlopt_fn), optargs))
    
    # Create output object
    model$param[1:3] = opt$solution
    
    output = list(
        par = opt$solution,
        model = model,
        events = events,
        end = end,
        opt = opt
    )
    
    class(output) = "HawkesModel"
    
    return( output )
}


# POWER LAW KERNEL, GAMMA = 2.5 -------------------------------------------


n.sim = 1e3

burn.in = 100
Ts = c(100, 150, 250, 400, 650, 1000, 2000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
shape = 2.5
scale = 1.5

cl = makePSOCKcluster(numCores, outfile = "log.txt") # use multicore, set to the number of our cores
clusterExport(cl, c("burn.in", "Ts", "binsizes", "lambda", "mu", "shape", "scale"))
registerDoParallel(cl)

estimates = foreach (sim = icount(n.sim), .combine = rbind, .multicombine = TRUE, .inorder = FALSE,
                     .packages = c("hawkesbow"), .errorhandling = "remove", .verbose = FALSE) %dopar% {
                         
                         # Counter display
                         if (sim %% 5 == 0) cat(sim, ". ", sep="")
                         
                         estimates = data.frame(T = numeric(length(Ts) * (1+length(binsizes))),
                                                simulation = integer(length(Ts) * (1+length(binsizes))),
                                                binsize = character(length(Ts) * (1+length(binsizes))),
                                                baseline = numeric(length(Ts) * (1+length(binsizes))),
                                                reprod = numeric(length(Ts) * (1+length(binsizes))),
                                                shape = numeric(length(Ts) * (1+length(binsizes))),
                                                scale = numeric(length(Ts) * (1+length(binsizes))),
                                                message = character(length(Ts) * (1+length(binsizes))))
                         
                         model = new(PowerLaw)
                         model$param[4] = scale
                         it = 1
                         
                         for (Tend in Ts) {
                             
                             # Simulate Hawkes process
                             x = hawkes(burn.in + Tend, fun = lambda, repr = mu, family = "powerlaw", shape = shape, scale = scale)
                             x$p = x$p[x$p > burn.in] - burn.in
                             x$end = Tend
                             
                             # MLE
                             opt = mle(x$p, model, Tend, lb = rep(0.05, 3), ub = c(50, 0.95, 50), init = c(lambda, mu, shape))
                             
                             estimates$T[it] = Tend
                             estimates$simulation[it] = sim
                             estimates$binsize[it] = "MLE"
                             estimates$baseline[it] = opt$par[1]
                             estimates$reprod[it] = opt$par[2]
                             estimates$shape[it] = opt$par[3]
                             estimates$scale[it] = opt$par[4]
                             estimates$message[it] = opt$opt$message
                             
                             it = it + 1
                             
                             for (binsize in binsizes) {
                                 
                                 # Whittle
                                 y = discrete(x, binsize = binsize)
                                 opt = whittle(y, model, binsize, lower = rep(0.05, 3), upper = c(50, 0.95, 50), 
                                               control = list(fnscale = 1e3, parscale = c(1,1,5), ndeps = rep(1e-4, 3)))
                                 
                                 estimates$T[it] = Tend
                                 estimates$simulation[it] = sim
                                 estimates$binsize[it] = as.character(binsize)
                                 estimates$baseline[it] = opt$par[1]
                                 estimates$reprod[it] = opt$par[2]
                                 estimates$shape[it] = opt$par[3]
                                 estimates$scale[it] = opt$par[4]
                                 estimates$message[it] = opt$opt$message
                                 
                                 it = it + 1
                                 
                             }
                             
                         }
                         
                         return(estimates)
                         
                     } %>% as_tibble() %>% mutate(message = as.factor(message))

stopCluster(cl)

save(estimates, file = "powerlaw25_no_a_convergence3.RData")

estimates.long = estimates %>%
    select(-scale) %>% 
    gather("parameter", "estimate", baseline:shape) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts)), levels = paste0("T = ", Ts)),
           parameter = factor(parameter, levels = c("baseline", "reprod", "shape")))

trueval = tibble(T = rep(Ts, each = (1+length(binsizes)) * 3),
                 binsize = rep(c("MLE", as.character(binsizes)), length(Ts) * 3),
                 parameter = rep(rep(c("baseline", "reprod", "shape"), each = 1 + length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, shape), each = 1 + length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "shape")))

ggplot(estimates.long %>% 
           filter(T %in% paste0("T = ", c(100, 1000))), 
       aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=filter(trueval, T %in% paste0("T = ", c(100, 1000))), size=4, colour = "firebrick") +
    scale_x_discrete(labels=c("baseline"=expression(hat(eta)), "reprod"=expression(hat(mu)), "shape"=expression(hat(gamma)))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,4)) +
    xlab("Parameters") +
    ylab("Estimates") +
    # labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = \\gamma a^{\\gamma} (a+t)^{-\\gamma-1}$ with $\\gamma = 2.5$ and $a = 1.5$ on $(0, T)$ | true values are larger points")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("powerlaw25_no_a_estimates.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("powerlaw25_no_a_estimates.png", width=13.6, height=7.6, units="in")

estimates.variance = estimates %>%
    select(-scale) %>% 
    gather("parameter", "estimate", baseline:shape) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("bin size = ", binsizes))),
                         levels = c("Regular MLE", paste0("bin size = ", binsizes))),
        parameter = factor(parameter, levels = c("baseline", "reprod", "shape")),
        trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, shape))
    ) %>%
    group_by(T, binsize, parameter) %>%
    summarise(bias.squared = mean(estimate - trueval)^2, 
              variance = var(estimate), 
              MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>%
    filter(T == 100) %>%
    group_by(parameter) %>%
    summarise(min = min(MSE), max = max(MSE))

estimates.at100 = estimates.variance %>% 
    filter(T == 100)

## Modify estimates.T100 and estimates.at100 to get both of them
estimates.T100$max[2] = 0.0441
estimates.T100$min[1] = 0.132
estimates.at100 = estimates.at100 %>% filter(MSE < 0.1, MSE > 0.07)

ggplot(estimates.variance %>% filter(T != 50),
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point(size = 2) +
    geom_abline(intercept = -.1, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .0025, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    annotate(geom = "segment",
             x = as.numeric(as.character(
                 plyr::revalue(estimates.at100$parameter,
                               c("baseline"=82.65, "reprod"=90, "shape"=98))
             )),
             xend = 100,
             y = estimates.at100$MSE,
             yend = estimates.at100$MSE,
             ) +
    annotate(geom = "label",
             x = as.numeric(as.character(
                 plyr::revalue(estimates.at100$parameter,
                               c("baseline"=82.65, "reprod"=90, "shape"=98))
             )),
             y = estimates.at100$MSE,
             label = plyr::revalue(
                 estimates.at100$parameter,
                 c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma))
             ),
             hjust = 1, parse = TRUE
    ) +
    annotate(geom = "label",
             x = 90,
             y = sqrt(estimates.T100$min * estimates.T100$max),
             label = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma)),
             hjust = 1, size = 7
    ) +
    annotate(geom = "linerange", x = sqrt(90*100), ymin = estimates.T100$min * .925, ymax = estimates.T100$max * 1.075,
             size = 2, color = c("grey80", "grey50", "grey50")) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Bin size",
                       labels = c("MLE", "0.25", "0.5", "1", "2"),
                       values = c("firebrick", viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma))) +
    # guides(linetype = guide_legend(order = 1), shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
    xlab("T") +
    ylab("Mean Square Error") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("powerlaw25_no_a_convergence.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("powerlaw25_no_a_convergence.png", width=13.6, height=7.6, units="in")


# POWER LAW KERNEL, GAMMA = 1.5 -------------------------------------------


n.sim = 1e3

burn.in = 100
Ts = c(100, 150, 250, 400, 650, 1000, 2000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
shape = 1.5
scale = 1.5

cl = makePSOCKcluster(numCores, outfile = "log.txt") # use multicore, set to the number of our cores
clusterExport(cl, c("burn.in", "Ts", "binsizes", "lambda", "mu", "shape", "scale"))
registerDoParallel(cl)

estimates = foreach (sim = icount(n.sim), .combine = rbind, .multicombine = TRUE, .inorder = FALSE,
                     .packages = c("hawkesbow"), .errorhandling = "remove", .verbose = FALSE) %dopar% {
                         
                         # Counter display
                         if (sim %% 5 == 0) cat(sim, ". ", sep="")
                         
                         estimates = data.frame(T = numeric(length(Ts) * (1+length(binsizes))),
                                                simulation = integer(length(Ts) * (1+length(binsizes))),
                                                binsize = character(length(Ts) * (1+length(binsizes))),
                                                baseline = numeric(length(Ts) * (1+length(binsizes))),
                                                reprod = numeric(length(Ts) * (1+length(binsizes))),
                                                shape = numeric(length(Ts) * (1+length(binsizes))),
                                                scale = numeric(length(Ts) * (1+length(binsizes))),
                                                message = character(length(Ts) * (1+length(binsizes))))
                         
                         model = new(PowerLaw)
                         model$param[4] = scale
                         it = 1
                         
                         for (Tend in Ts) {
                             
                             # Simulate Hawkes process
                             x = hawkes(burn.in + Tend, fun = lambda, repr = mu, family = "powerlaw", shape = shape, scale = scale)
                             x$p = x$p[x$p > burn.in] - burn.in
                             x$end = Tend
                             
                             # MLE
                             opt = mle(x$p, model, Tend, lb = rep(0.05, 3), ub = c(50, 0.95, 50), init = c(lambda, mu, shape))
                             
                             estimates$T[it] = Tend
                             estimates$simulation[it] = sim
                             estimates$binsize[it] = "MLE"
                             estimates$baseline[it] = opt$par[1]
                             estimates$reprod[it] = opt$par[2]
                             estimates$shape[it] = opt$par[3]
                             estimates$scale[it] = opt$par[4]
                             estimates$message[it] = opt$opt$message
                             
                             it = it + 1
                             
                             for (binsize in binsizes) {
                                 
                                 # Whittle
                                 y = discrete(x, binsize = binsize)
                                 opt = whittle(y, model, binsize, lower = rep(0.05, 3), upper = c(50, 0.95, 50), 
                                               control = list(fnscale = 1e3, parscale = c(1,1,5), ndeps = rep(1e-4, 3)))
                                 
                                 estimates$T[it] = Tend
                                 estimates$simulation[it] = sim
                                 estimates$binsize[it] = as.character(binsize)
                                 estimates$baseline[it] = opt$par[1]
                                 estimates$reprod[it] = opt$par[2]
                                 estimates$shape[it] = opt$par[3]
                                 estimates$scale[it] = opt$par[4]
                                 estimates$message[it] = opt$opt$message
                                 
                                 it = it + 1
                                 
                             }
                             
                         }
                         
                         return(estimates)
                         
                     } %>% as_tibble() %>% mutate(message = as.factor(message))

stopCluster(cl)

save(estimates, file = "powerlaw15_no_a_convergence.RData")

estimates.long = estimates %>%
    select(-scale) %>% 
    gather("parameter", "estimate", baseline:shape) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts)), levels = paste0("T = ", Ts)),
           parameter = factor(parameter, levels = c("baseline", "reprod", "shape")))

trueval = tibble(T = rep(Ts, each = (1+length(binsizes)) * 3),
                 binsize = rep(c("MLE", as.character(binsizes)), length(Ts) * 3),
                 parameter = rep(rep(c("baseline", "reprod", "shape"), each = 1 + length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, shape), each = 1 + length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "shape")))

ggplot(estimates.long %>% 
           filter(T %in% paste0("T = ", c(100, 1000))), 
       aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=filter(trueval, T %in% paste0("T = ", c(100, 1000))), size=4, colour = "firebrick") +
    scale_x_discrete(labels=c("baseline"=expression(hat(eta)), "reprod"=expression(hat(mu)), "shape"=expression(hat(gamma)))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,4)) +
    xlab("Parameters") +
    ylab("Estimates") +
    # labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = \\gamma a^{\\gamma} (a+t)^{-\\gamma-1}$ with $\\gamma = 1.5$ and $a = 1.5$ on $(0, T)$ | true values are larger points")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("powerlaw15_no_a_estimates.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("powerlaw15_no_a_estimates.png", width=13.6, height=7.6, units="in")

estimates.variance = estimates %>%
    select(-scale) %>% 
    gather("parameter", "estimate", baseline:shape) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("bin size = ", binsizes))),
                         levels = c("Regular MLE", paste0("bin size = ", binsizes))),
        parameter = factor(parameter, levels = c("baseline", "reprod", "shape")),
        trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, shape))
    ) %>%
    group_by(T, binsize, parameter) %>%
    summarise(bias.squared = mean(estimate - trueval)^2, 
              variance = var(estimate), 
              MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>%
    filter(T == 100) %>%
    group_by(parameter) %>%
    summarise(min = min(MSE), max = max(MSE))

estimates.at100 = estimates.variance %>% 
    filter(T == 100)

ggplot(estimates.variance %>% filter(T != 50),
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point(size = 2) +
    geom_abline(intercept = .2, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .005, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    # annotate(geom = "segment",
    #          x = as.numeric(as.character(
    #              plyr::revalue(estimates.at100$parameter,
    #                            c("baseline"=82.65, "reprod"=90, "shape"=98))
    #          )),
    #          xend = 100,
    #          y = estimates.at100$MSE,
    #          yend = estimates.at100$MSE,
    #          color = as.character(plyr::mapvalues(
    #              estimates.at100$parameter,
    #              from = c("baseline", "reprod", "shape"),
    #              to = viridis::magma(3, end = 0.75)
    #          ))) +
    # annotate(geom = "label",
    #          x = as.numeric(as.character(
    #              plyr::revalue(estimates.at100$parameter,
    #                            c("baseline"=82.65, "reprod"=90, "shape"=98))
    #          )),
    #          y = estimates.at100$MSE,
    #          label = plyr::revalue(
    #              estimates.at100$parameter,
    #              c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma))
    #          ),
    #          color = as.character(plyr::mapvalues(
    #              estimates.at100$parameter,
    #              from = c("baseline", "reprod", "shape"),
    #              to = viridis::magma(3, end = 0.75)
    #          )),
    #          hjust = 1, parse = TRUE
    # ) +
    annotate(geom = "label",
             x = 90,
             y = sqrt(estimates.T100$min * estimates.T100$max),
             label = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma)),
             hjust = 1, size = 7
    ) +
    annotate(geom = "linerange", x = sqrt(90*100), ymin = estimates.T100$min * .925, ymax = estimates.T100$max * 1.075,
             size = 2, color = c("grey80", "grey50", "grey50")) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Bin size",
                       labels = c("MLE", "0.25", "0.5", "1", "2"),
                       values = c("firebrick", viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma))) +
    # guides(linetype = guide_legend(order = 1), shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
    xlab("T") +
    ylab("Mean Square Error") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("powerlaw15_no_a_convergence.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("powerlaw15_no_a_convergence.png", width=13.6, height=7.6, units="in")


# POWER LAW KERNEL, GAMMA = 0.5 -------------------------------------------


n.sim = 1e3

burn.in = 100
Ts = c(100, 150, 250, 400, 650, 1000, 2000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
shape = 0.5
scale = 1.5

cl = makePSOCKcluster(numCores, outfile = "log.txt") # use multicore, set to the number of our cores
clusterExport(cl, c("burn.in", "Ts", "binsizes", "lambda", "mu", "shape", "scale"))
registerDoParallel(cl)

estimates = foreach (sim = icount(n.sim), .combine = rbind, .multicombine = TRUE, .inorder = FALSE,
                     .packages = c("hawkesbow"), .errorhandling = "remove", .verbose = FALSE) %dopar% {
                         
                         # Counter display
                         if (sim %% 5 == 0) cat(sim, ". ", sep="")
                         
                         estimates = data.frame(T = numeric(length(Ts) * (1+length(binsizes))),
                                                simulation = integer(length(Ts) * (1+length(binsizes))),
                                                binsize = character(length(Ts) * (1+length(binsizes))),
                                                baseline = numeric(length(Ts) * (1+length(binsizes))),
                                                reprod = numeric(length(Ts) * (1+length(binsizes))),
                                                shape = numeric(length(Ts) * (1+length(binsizes))),
                                                scale = numeric(length(Ts) * (1+length(binsizes))),
                                                message = character(length(Ts) * (1+length(binsizes))))
                         
                         model = new(PowerLaw)
                         model$param[4] = scale
                         it = 1
                         
                         for (Tend in Ts) {
                             
                             # Simulate Hawkes process
                             x = hawkes(burn.in + Tend, fun = lambda, repr = mu, family = "powerlaw", shape = shape, scale = scale)
                             x$p = x$p[x$p > burn.in] - burn.in
                             x$end = Tend
                             
                             # MLE
                             opt = mle(x$p, model, Tend, lb = rep(0.05, 3), ub = c(50, 0.95, 50), init = c(lambda, mu, shape))
                             
                             estimates$T[it] = Tend
                             estimates$simulation[it] = sim
                             estimates$binsize[it] = "MLE"
                             estimates$baseline[it] = opt$par[1]
                             estimates$reprod[it] = opt$par[2]
                             estimates$shape[it] = opt$par[3]
                             estimates$scale[it] = opt$par[4]
                             estimates$message[it] = opt$opt$message
                             
                             it = it + 1
                             
                             for (binsize in binsizes) {
                                 
                                 # Whittle
                                 y = discrete(x, binsize = binsize)
                                 opt = whittle(y, model, binsize, lower = rep(0.05, 3), upper = c(50, 0.95, 50), 
                                               control = list(fnscale = 1e3, parscale = c(1,1,5), ndeps = rep(1e-4, 3)))
                                 
                                 estimates$T[it] = Tend
                                 estimates$simulation[it] = sim
                                 estimates$binsize[it] = as.character(binsize)
                                 estimates$baseline[it] = opt$par[1]
                                 estimates$reprod[it] = opt$par[2]
                                 estimates$shape[it] = opt$par[3]
                                 estimates$scale[it] = opt$par[4]
                                 estimates$message[it] = opt$opt$message
                                 
                                 it = it + 1
                                 
                             }
                             
                         }
                         
                         return(estimates)
                         
                     } %>% as_tibble() %>% mutate(message = as.factor(message))

stopCluster(cl)

save(estimates, file = "powerlaw05_no_a_convergence.RData")

estimates.long = estimates %>%
    select(-scale) %>% 
    gather("parameter", "estimate", baseline:shape) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts)), levels = paste0("T = ", Ts)),
           parameter = factor(parameter, levels = c("baseline", "reprod", "shape")))

trueval = tibble(T = rep(Ts, each = (1+length(binsizes)) * 3),
                 binsize = rep(c("MLE", as.character(binsizes)), length(Ts) * 3),
                 parameter = rep(rep(c("baseline", "reprod", "shape"), each = 1 + length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, shape), each = 1 + length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("bin size = ", binsizes))),
                            levels = c("Regular MLE", paste0("bin size = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "shape")))

ggplot(estimates.long %>% 
           filter(T %in% paste0("T = ", c(100, 1000))), 
       aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=filter(trueval, T %in% paste0("T = ", c(100, 1000))), size=4, colour = "firebrick") +
    scale_x_discrete(labels=c("baseline"=expression(hat(eta)), "reprod"=expression(hat(mu)), "shape"=expression(hat(gamma)))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,4)) +
    xlab("Parameters") +
    ylab("Estimates") +
    # labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = \\gamma a^{\\gamma} (a+t)^{-\\gamma-1}$ with $\\gamma = 0.5$ and $a = 1.5$ on $(0, T)$ | true values are larger points")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("powerlaw05_no_a_estimates.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("powerlaw05_no_a_estimates.png", width=13.6, height=7.6, units="in")

estimates.variance = estimates %>%
    select(-scale) %>% 
    gather("parameter", "estimate", baseline:shape) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("bin size = ", binsizes))),
                         levels = c("Regular MLE", paste0("bin size = ", binsizes))),
        parameter = factor(parameter, levels = c("baseline", "reprod", "shape")),
        trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, shape))
    ) %>%
    group_by(T, binsize, parameter) %>%
    summarise(bias.squared = mean(estimate - trueval)^2, 
              variance = var(estimate), 
              MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>%
    filter(T == 100) %>%
    group_by(parameter) %>%
    summarise(min = min(MSE), max = max(MSE))

estimates.at100 = estimates.variance %>% 
    filter(T == 100)

ggplot(estimates.variance %>% filter(T != 50),
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point(size = 2) +
    geom_abline(intercept = .75, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .0175, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    # annotate(geom = "segment",
    #          x = as.numeric(as.character(
    #              plyr::revalue(estimates.at100$parameter,
    #                            c("baseline"=82.65, "reprod"=90, "shape"=98))
    #          )),
    #          xend = 100,
    #          y = estimates.at100$MSE,
    #          yend = estimates.at100$MSE,
    #          color = as.character(plyr::mapvalues(
    #              estimates.at100$parameter,
    #              from = c("baseline", "reprod", "shape"),
    #              to = viridis::magma(3, end = 0.75)
    #          ))) +
    # annotate(geom = "label",
    #          x = as.numeric(as.character(
    #              plyr::revalue(estimates.at100$parameter,
    #                            c("baseline"=82.65, "reprod"=90, "shape"=98))
    #          )),
    #          y = estimates.at100$MSE,
    #          label = plyr::revalue(
    #              estimates.at100$parameter,
    #              c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma))
    #          ),
    #          color = as.character(plyr::mapvalues(
    #              estimates.at100$parameter,
    #              from = c("baseline", "reprod", "shape"),
    #              to = viridis::magma(3, end = 0.75)
    #          )),
    #          hjust = 1, parse = TRUE
    # ) +
    annotate(geom = "label",
             x = 90,
             y = sqrt(estimates.T100$min * estimates.T100$max),
             label = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma)),
             hjust = 1, size = 7
    ) +
    annotate(geom = "linerange", x = sqrt(90*100), ymin = estimates.T100$min * .925, ymax = estimates.T100$max * 1.075,
             size = 2, color = c("grey80", "grey50", "grey50")) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Bin size",
                       labels = c("MLE", "0.25", "0.5", "1", "2"),
                       values = c("firebrick", viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "shape"=expression(gamma))) +
    # guides(linetype = guide_legend(order = 1), shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
    xlab("T") +
    ylab("Mean Square Error") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("powerlaw05_no_a_convergence.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("powerlaw05_no_a_convergence.png", width=13.6, height=7.6, units="in")
