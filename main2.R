# !diagnostics off

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(ggplot2)
library(hawkesbow)
library(latex2exp)
library(PtProcess)

##########################
## EXPONENTIAL KERNEL

n.sim = 1e3

Ts = c(100, 1000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
rate = 1
init = c(.5, .8, 2)

model = new(Exponential)
estimates = tibble(T = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   simulation = integer(n.sim * length(Ts) * (1+length(binsizes))),
                   binsize = character(n.sim * length(Ts) * (1+length(binsizes))),
                   baseline = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   reprod = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   rate = numeric(n.sim * length(Ts) * (1+length(binsizes))))

it = 1

for (T in Ts) {

    # Counter display
    display = paste0("Starting simulations for T = ", T)
    cat("\n", rep("=", nchar(display)), "\n", display, "\n", sep = "")

    for (sim in 1:n.sim) {

        # Counter display
        if (sim %% 10 == 0) cat(sim, ". ", sep="")

        # Simulate Hawkes process
        x = hawkes(T, fun = lambda, repr = mu, family = "exp", rate = rate)

        # MLE
        data = new(ContinuousData, x$p, x$T)
        model$param = init
        opt = mle(model, data)

        estimates$T[it] = T
        estimates$simulation[it] = sim
        estimates$binsize[it] = "MLE"
        estimates$baseline[it] = opt$solution[1]
        estimates$reprod[it] = opt$solution[2]
        estimates$rate[it] = opt$solution[3]

        it = it + 1

        for (binsize in binsizes) {

            # Whittle
            y = discrete(x, binsize = binsize)
            data = new(DiscreteData, y, binsize)
            model$param = init
            opt = whittle(model, data)

            estimates$T[it] = T
            estimates$simulation[it] = sim
            estimates$binsize[it] = as.character(binsize)
            estimates$baseline[it] = opt$par[1]
            estimates$reprod[it] = opt$par[2]
            estimates$rate[it] = opt$par[3]

            it = it + 1

        }

    }

}

estimates.long = estimates %>% gather("parameter", "estimate", baseline:rate) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("binsize = ", binsizes))),
                            levels = c("Regular MLE", paste0("binsize = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "rate")))

trueval = tibble(T = rep(Ts, each = (1+length(binsizes)) * length(init)),
                 binsize = rep(c("MLE", as.character(binsizes)), length(Ts) * length(init)),
                 parameter = rep(rep(c("baseline", "reprod", "rate"), each = 1 + length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, rate), each = 1 + length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("binsize = ", binsizes))),
                            levels = c("Regular MLE", paste0("binsize = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "rate")))

ggplot(estimates.long, aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=trueval, col="red", linetype=3) +
    scale_x_discrete(labels=c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,4)) +
    xlab("Parameters") +
    ylab("Estimates") +
    labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = 1 e^{-1t}$ on $(0, T)$ | true values in red")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("exponential_estimates.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("exponential_estimates.png", width=13.6, height=7.6, units="in")


############################
## CONVERGENCE SPEED

n.sim = 1e3

Ts = c(50, 100, 150, 250, 400, 650, 1000, 2000, 4000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
rate = 1
init = c(.5, .8, 2)

model = new(Exponential)
estimates = tibble(T = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   simulation = integer(n.sim * length(Ts) * (1+length(binsizes))),
                   binsize = character(n.sim * length(Ts) * (1+length(binsizes))),
                   baseline = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   reprod = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   rate = numeric(n.sim * length(Ts) * (1+length(binsizes))))

it = 1

for (T in Ts) {
    
    # Counter display
    display = paste0("Starting simulations for T = ", T)
    cat("\n", rep("=", nchar(display)), "\n", display, "\n", sep = "")
    
    for (sim in 1:n.sim) {
        
        # Counter display
        if (sim %% 10 == 0) cat(sim, ". ", sep="")
        
        # Simulate Hawkes process
        x = hawkes(T, fun = lambda, repr = mu, family = "exp", rate = rate)
        
        # MLE
        data = new(ContinuousData, x$p, x$T)
        model$param = init
        opt = mle(model, data)
        
        estimates$T[it] = T
        estimates$simulation[it] = sim
        estimates$binsize[it] = "MLE"
        estimates$baseline[it] = opt$solution[1]
        estimates$reprod[it] = opt$solution[2]
        estimates$rate[it] = opt$solution[3]
        
        it = it + 1
        
        for (binsize in binsizes) {
            
            # Whittle
            y = discrete(x, binsize = binsize)
            data = new(DiscreteData, y, binsize)
            model$param = init
            opt = whittle(model, data)
            
            estimates$T[it] = T
            estimates$simulation[it] = sim
            estimates$binsize[it] = as.character(binsize)
            estimates$baseline[it] = opt$par[1]
            estimates$reprod[it] = opt$par[2]
            estimates$rate[it] = opt$par[3]
            
            it = it + 1
            
        }
        
    }
    
}

estimates.long = estimates %>% mutate(T = as.factor(T)) %>% 
    gather("parameter", "estimate", baseline:rate) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("binsize = ", binsizes))),
                            levels = c("Regular MLE", paste0("binsize = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "rate")))

trueval = tibble(T = rep(Ts, each = (1+length(binsizes)) * length(init)),
                 binsize = rep(c("MLE", as.character(binsizes)), length(Ts) * length(init)),
                 parameter = rep(rep(c("baseline", "reprod", "rate"), each = 1 + length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, rate), each = 1 + length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                            to = c("Regular MLE", paste0("binsize = ", binsizes))),
                            levels = c("Regular MLE", paste0("binsize = ", binsizes))),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "rate")))

ggplot(estimates.long, aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=trueval, col="red") +
    scale_x_discrete(labels=c("baseline"=expression(eta), "reprod"=expression(mu), "rate"=expression(beta))) +
    facet_grid(binsize ~ T) +
    coord_cartesian(ylim=c(0,3)) +
    xlab("Parameters") +
    ylab("Estimates") +
    labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = 1 e^{-1t}$ on $(0, T)$ | true values in red")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

estimates.variance = estimates %>% 
    gather("parameter", "estimate", baseline:rate) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("binsize = ", binsizes))),
                         levels = c("Regular MLE", paste0("binsize = ", binsizes))),
       parameter = factor(parameter, levels = c("baseline", "reprod", "rate")),
       trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, rate))
    ) %>% 
    group_by(T, binsize, parameter) %>% 
    summarise(variance = var(estimate), MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>% 
    filter(T == 100) %>% 
    group_by(parameter) %>% 
    summarise(min = min(MSE), max = max(MSE))

ggplot(estimates.variance %>% filter(T != 50), 
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point() +
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
    scale_color_manual(name = "Binsize", 
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

ggsave("exponential_convergence3.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("exponential_convergence3.png", width=13.6, height=7.6, units="in")

estimates.standard = estimates.long %>% 
    inner_join(trueval %>% rename(trueval = estimate)) %>% 
    inner_join(estimates.long %>% group_by(T, binsize, parameter) %>% summarise(stdev = sd(estimate))) %>% 
    mutate(stdEstimate = (estimate - trueval) / stdev) %>% 
    select(- estimate, - trueval, - stdev)

ggplot(estimates.standard %>% filter(T == "T = 1000"), aes(sample = stdEstimate)) +
    geom_qq() +
    geom_qq_line(col = "red") +
    facet_grid(parameter ~ binsize) +
    theme_bw()

############################
## PARETO3 KERNEL

n.sim = 1e3

Ts = c(100, 150, 250, 400, 650, 1000, 2000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
rate = 1

lambda = lambda
mu = mu
a = 2 / (3 * rate)
init = c(.5, .8, 1)

model = new(Pareto3)
estimates = tibble(T = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   simulation = integer(n.sim * length(Ts) * (1+length(binsizes))),
                   binsize = character(n.sim * length(Ts) * (1+length(binsizes))),
                   baseline = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   reprod = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   position = numeric(n.sim * length(Ts) * (1+length(binsizes))))

it = 1

for (T in Ts) {

    # Counter display
    display = paste0("Starting simulations for T = ", T)
    cat("\n", rep("=", nchar(display)), "\n", display, "\n", sep = "")

    for (sim in 1:n.sim) {

        # Counter display
        if (sim %% 10 == 0) cat(sim, ". ", sep="")

        # Simulate Hawkes process
        x = hawkes(T, fun = lambda, repr = mu, family = "pareto", lambda = 3, a = a)

        # # MLE
        # data = new(ContinuousData, x$p, x$T)
        # model$param = init
        # opt = mle(model, data)
        #
        estimates$T[it] = T
        estimates$simulation[it] = sim
        estimates$binsize[it] = "MLE"
        # estimates$baseline[it] = opt$solution[1]
        # estimates$reprod[it] = opt$solution[2]
        # estimates$position[it] = opt$solution[3]

        it = it + 1

        tryCatch(expr = {
            for (binsize in binsizes) {

                # Whittle
                y = discrete(x, binsize = binsize)
                data = new(DiscreteData, y, binsize)
                model$param = init
                opt = whittle(model, data, lower = c(.0001, .0001, .01))

                estimates$T[it] = T
                estimates$simulation[it] = sim
                estimates$binsize[it] = as.character(binsize)
                estimates$baseline[it] = opt$par[1]
                estimates$reprod[it] = opt$par[2]
                estimates$position[it] = opt$par[3]

                it = it + 1
            }},
            error = function(cond) {
                break
        })

    }

}

estimates.long = estimates %>% 
    filter(binsize != "MLE") %>% 
    gather("parameter", "estimate", baseline:position) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = binsizes,
                                            to = paste0("binsize = ", binsizes)),
                            levels = paste0("binsize = ", binsizes)),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts)),
                      levels = paste0("T = ", Ts)),
           parameter = factor(parameter, levels = c("baseline", "reprod", "position")))

trueval = tibble(T = rep(Ts, each = length(binsizes) * length(init)),
                 binsize = rep(as.character(binsizes), length(Ts) * length(init)),
                 parameter = rep(rep(c("baseline", "reprod", "position"), each = length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, a), each = length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = binsizes,
                                            to = paste0("binsize = ", binsizes)),
                            levels = paste0("binsize = ", binsizes)),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "position")))

ggplot(estimates.long %>% filter(T %in% paste0("T = ", c(100, 1000))), 
       aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=trueval %>% filter(T %in% paste0("T = ", c(100, 1000))), 
               col="red", linetype=3) +
    scale_x_discrete(labels=c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,2)) +
    xlab("Parameters") +
    ylab("Estimates") +
    labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = 3 (2/3)^3 t^{-4}$ on $(0, T)$ | true values in red")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("pareto3_estimates3.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("pareto3_estimates3.png", width=13.6, height=7.6, units="in")

estimates.variance = estimates %>% 
    gather("parameter", "estimate", baseline:position) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("binsize = ", binsizes))),
                         levels = c("Regular MLE", paste0("binsize = ", binsizes))),
        parameter = factor(parameter, levels = c("baseline", "reprod", "position")),
        trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, a))
    ) %>% 
    group_by(T, binsize, parameter) %>% 
    summarise(variance = var(estimate), MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>% 
    filter(T == 100, binsize != "Regular MLE")

ggplot(estimates.variance %>% filter(binsize != "Regular MLE"), 
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point() +
    geom_abline(intercept = -.1, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .0025, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    annotate(geom = "segment", 
             x = as.numeric(as.character(
                 plyr::revalue(estimates.T100$parameter, 
                               c("baseline"=82.65, "reprod"=90, "position"=98))
                 )),
             xend = 100,
             y = estimates.T100$MSE,
             yend = estimates.T100$MSE,
             color = as.character(plyr::mapvalues(
                 estimates.T100$parameter,
                 from = c("baseline", "reprod", "position"),
                 to = viridis::magma(3, end = 0.75)
             ))) +
    annotate(geom = "label", 
             x = as.numeric(as.character(
                 plyr::revalue(estimates.T100$parameter,
                               c("baseline"=82.65, "reprod"=90, "position"=98))
                 )),
             y = estimates.T100$MSE,
             label = plyr::revalue(
                 estimates.T100$parameter,
                 c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))
             ),
             color = as.character(plyr::mapvalues(
                 estimates.T100$parameter,
                 from = c("baseline", "reprod", "position"),
                 to = viridis::magma(3, end = 0.75)
             )),
             hjust = 1, parse = TRUE
    ) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Binsize", 
                       labels = c("0.25", "0.5", "1", "2"),
                       values = c(viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))) +
    xlab("T") +
    ylab("Mean Square Error") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("pareto3_convergence3.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("pareto3_convergence3.png", width=13.6, height=7.6, units="in")

############################
## PARETO2 KERNEL

n.sim = 1e3

Ts = c(100, 150, 250, 400, 650, 1000, 2000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
rate = 1

lambda = lambda
mu = mu
a = 1 / (2 * rate)
init = c(.5, .8, 1)

model = new(Pareto2)
estimates = tibble(T = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   simulation = integer(n.sim * length(Ts) * (1+length(binsizes))),
                   binsize = character(n.sim * length(Ts) * (1+length(binsizes))),
                   baseline = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   reprod = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   position = numeric(n.sim * length(Ts) * (1+length(binsizes))))

it = 1

for (T in Ts) {
    
    # Counter display
    display = paste0("Starting simulations for T = ", T)
    cat("\n", rep("=", nchar(display)), "\n", display, "\n", sep = "")
    
    for (sim in 1:n.sim) {
        
        # Counter display
        if (sim %% 10 == 0) cat(sim, ". ", sep="")
        
        # Simulate Hawkes process
        x = hawkes(T, fun = lambda, repr = mu, family = "pareto", lambda = 2, a = a)
        
        # # MLE
        # data = new(ContinuousData, x$p, x$T)
        # model$param = init
        # opt = mle(model, data)
        #
        estimates$T[it] = T
        estimates$simulation[it] = sim
        estimates$binsize[it] = "MLE"
        # estimates$baseline[it] = opt$solution[1]
        # estimates$reprod[it] = opt$solution[2]
        # estimates$position[it] = opt$solution[3]
        
        it = it + 1
        
        tryCatch(expr = {
            for (binsize in binsizes) {
                
                # Whittle
                y = discrete(x, binsize = binsize)
                data = new(DiscreteData, y, binsize)
                model$param = init
                opt = whittle(model, data, lower = c(.0001, .0001, .01))
                
                estimates$T[it] = T
                estimates$simulation[it] = sim
                estimates$binsize[it] = as.character(binsize)
                estimates$baseline[it] = opt$par[1]
                estimates$reprod[it] = opt$par[2]
                estimates$position[it] = opt$par[3]
                
                it = it + 1
            }},
            error = function(cond) {
                break
            })
        
    }
    
}

estimates.long = estimates %>% 
    filter(binsize != "MLE") %>% 
    gather("parameter", "estimate", baseline:position) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = binsizes,
                                            to = paste0("binsize = ", binsizes)),
                            levels = paste0("binsize = ", binsizes)),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts)),
                      levels = paste0("T = ", Ts)),
           parameter = factor(parameter, levels = c("baseline", "reprod", "position")))

trueval = tibble(T = rep(Ts, each = length(binsizes) * length(init)),
                 binsize = rep(as.character(binsizes), length(Ts) * length(init)),
                 parameter = rep(rep(c("baseline", "reprod", "position"), each = length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, a), each = length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = binsizes,
                                            to = paste0("binsize = ", binsizes)),
                            levels = paste0("binsize = ", binsizes)),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "position")))

ggplot(estimates.long %>% filter(T %in% paste0("T = ", c(100, 1000))), 
       aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=trueval %>% filter(T %in% paste0("T = ", c(100, 1000))), 
               col="red", linetype=3) +
    scale_x_discrete(labels=c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,2)) +
    xlab("Parameters") +
    ylab("Estimates") +
    labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = 2 (1/2)^2 t^{-3}$ on $(0, T)$ | true values in red")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("pareto2_estimates.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("pareto2_estimates.png", width=13.6, height=7.6, units="in")

estimates.variance = estimates %>% 
    gather("parameter", "estimate", baseline:position) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("binsize = ", binsizes))),
                         levels = c("Regular MLE", paste0("binsize = ", binsizes))),
        parameter = factor(parameter, levels = c("baseline", "reprod", "position")),
        trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, a))
    ) %>% 
    group_by(T, binsize, parameter) %>% 
    summarise(variance = var(estimate), MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>% 
    filter(T == 100, binsize != "Regular MLE")

ggplot(estimates.variance %>% filter(binsize != "Regular MLE"), 
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point() +
    geom_abline(intercept = -.1, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .0025, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    annotate(geom = "segment", 
             x = as.numeric(as.character(
                 plyr::revalue(estimates.T100$parameter, 
                               c("baseline"=82.65, "reprod"=90, "position"=98))
             )),
             xend = 100,
             y = estimates.T100$MSE,
             yend = estimates.T100$MSE,
             color = as.character(plyr::mapvalues(
                 estimates.T100$parameter,
                 from = c("baseline", "reprod", "position"),
                 to = viridis::magma(3, end = 0.75)
             ))) +
    annotate(geom = "label", 
             x = as.numeric(as.character(
                 plyr::revalue(estimates.T100$parameter,
                               c("baseline"=82.65, "reprod"=90, "position"=98))
             )),
             y = estimates.T100$MSE,
             label = plyr::revalue(
                 estimates.T100$parameter,
                 c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))
             ),
             color = as.character(plyr::mapvalues(
                 estimates.T100$parameter,
                 from = c("baseline", "reprod", "position"),
                 to = viridis::magma(3, end = 0.75)
             )),
             hjust = 1, parse = TRUE
    ) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Binsize", 
                       labels = c("0.25", "0.5", "1", "2"),
                       values = c(viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))) +
    xlab("T") +
    ylab("Mean Square Error") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("pareto2_convergence.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("pareto2_convergence.png", width=13.6, height=7.6, units="in")

############################
## PARETO1 KERNEL

n.sim = 1e3

Ts = c(100, 150, 250, 400, 650, 1000, 2000)
binsizes = c(.25, .5, 1, 2)

lambda = 1
mu = .5
rate = 1

lambda = lambda
mu = mu
a = 1 / (3 * rate)
init = c(.5, .8, 1)

model = new(Pareto1)
estimates = tibble(T = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   simulation = integer(n.sim * length(Ts) * (1+length(binsizes))),
                   binsize = character(n.sim * length(Ts) * (1+length(binsizes))),
                   baseline = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   reprod = numeric(n.sim * length(Ts) * (1+length(binsizes))),
                   position = numeric(n.sim * length(Ts) * (1+length(binsizes))))

it = 1

for (T in Ts) {
    
    # Counter display
    display = paste0("Starting simulations for T = ", T)
    cat("\n", rep("=", nchar(display)), "\n", display, "\n", sep = "")
    
    for (sim in 1:n.sim) {
        
        # Counter display
        if (sim %% 10 == 0) cat(sim, ". ", sep="")
        
        # Simulate Hawkes process
        x = hawkes(T, fun = lambda, repr = mu, family = "pareto", lambda = 1, a = a)
        
        # # MLE
        # data = new(ContinuousData, x$p, x$T)
        # model$param = init
        # opt = mle(model, data)
        #
        estimates$T[it] = T
        estimates$simulation[it] = sim
        estimates$binsize[it] = "MLE"
        # estimates$baseline[it] = opt$solution[1]
        # estimates$reprod[it] = opt$solution[2]
        # estimates$position[it] = opt$solution[3]
        
        it = it + 1
        
        tryCatch(expr = {
            for (binsize in binsizes) {
                
                # Whittle
                y = discrete(x, binsize = binsize)
                data = new(DiscreteData, y, binsize)
                model$param = init
                opt = whittle(model, data, lower = c(.01, .01, .01), upper = c(Inf, 0.99, Inf))
                
                estimates$T[it] = T
                estimates$simulation[it] = sim
                estimates$binsize[it] = as.character(binsize)
                estimates$baseline[it] = opt$par[1]
                estimates$reprod[it] = opt$par[2]
                estimates$position[it] = opt$par[3]
                
                it = it + 1
            }},
            error = function(cond) {
                break
            })
        
    }
    
}

estimates.long = estimates %>% 
    filter(binsize != "MLE") %>% 
    gather("parameter", "estimate", baseline:position) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = binsizes,
                                            to = paste0("binsize = ", binsizes)),
                            levels = paste0("binsize = ", binsizes)),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts)),
                      levels = paste0("T = ", Ts)),
           parameter = factor(parameter, levels = c("baseline", "reprod", "position")))

trueval = tibble(T = rep(Ts, each = length(binsizes) * length(init)),
                 binsize = rep(as.character(binsizes), length(Ts) * length(init)),
                 parameter = rep(rep(c("baseline", "reprod", "position"), each = length(binsizes)), length(Ts)),
                 estimate = rep(rep(c(lambda, mu, a), each = length(binsizes)), length(Ts))) %>%
    mutate(binsize = factor(plyr::mapvalues(binsize, from = binsizes,
                                            to = paste0("binsize = ", binsizes)),
                            levels = paste0("binsize = ", binsizes)),
           T = factor(plyr::mapvalues(T, from = Ts, to = paste0("T = ", Ts))),
           parameter = factor(parameter, levels = c("baseline", "reprod", "position")))

ggplot(estimates.long %>% filter(T %in% paste0("T = ", c(100, 1000))), 
       aes(x=parameter, y=estimate)) +
    geom_boxplot(fill="grey92") +
    geom_point(data=trueval %>% filter(T %in% paste0("T = ", c(100, 1000))), 
               col="red", linetype=3) +
    scale_x_discrete(labels=c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))) +
    facet_grid(T ~ binsize) +
    coord_cartesian(ylim=c(0,2)) +
    xlab("Parameters") +
    ylab("Estimates") +
    labs(title = TeX("$\\eta = 1,\\,\\mu = 0.5,\\, h^* (t) = 1 (1/3)^1 t^{-2}$ on $(0, T)$ | true values in red")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("pareto1_estimates.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("pareto1_estimates.png", width=13.6, height=7.6, units="in")

estimates.variance = estimates %>% 
    gather("parameter", "estimate", baseline:position) %>%
    mutate(
        binsize = factor(plyr::mapvalues(binsize, from = c("MLE", binsizes),
                                         to = c("Regular MLE", paste0("binsize = ", binsizes))),
                         levels = c("Regular MLE", paste0("binsize = ", binsizes))),
        parameter = factor(parameter, levels = c("baseline", "reprod", "position")),
        trueval = ifelse(parameter == "baseline", lambda, ifelse(parameter == "reprod", mu, a))
    ) %>% 
    group_by(T, binsize, parameter) %>% 
    summarise(variance = var(estimate), MSE = mean((estimate - trueval)^2))

estimates.T100 = estimates.variance %>% 
    filter(T == 100, binsize != "Regular MLE")

ggplot(estimates.variance %>% filter(binsize != "Regular MLE"), 
       aes(x=T, y=MSE, shape=parameter, linetype=parameter, color=binsize)) +
    geom_line() +
    geom_point() +
    geom_abline(intercept = -.1, slope = -1, col = "grey40", lty = "longdash") +
    annotate(geom = "text", x = 250, y = .0025, angle = -11.25, label =" Slope = -1", col = "grey40", size = 5) +
    annotate(geom = "segment", 
             x = as.numeric(as.character(
                 plyr::revalue(estimates.T100$parameter, 
                               c("baseline"=82.65, "reprod"=90, "position"=98))
             )),
             xend = 100,
             y = estimates.T100$MSE,
             yend = estimates.T100$MSE,
             color = as.character(plyr::mapvalues(
                 estimates.T100$parameter,
                 from = c("baseline", "reprod", "position"),
                 to = viridis::magma(3, end = 0.75)
             ))) +
    annotate(geom = "label", 
             x = as.numeric(as.character(
                 plyr::revalue(estimates.T100$parameter,
                               c("baseline"=82.65, "reprod"=90, "position"=98))
             )),
             y = estimates.T100$MSE,
             label = plyr::revalue(
                 estimates.T100$parameter,
                 c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))
             ),
             color = as.character(plyr::mapvalues(
                 estimates.T100$parameter,
                 from = c("baseline", "reprod", "position"),
                 to = viridis::magma(3, end = 0.75)
             )),
             hjust = 1, parse = TRUE
    ) +
    scale_x_log10() + # scale_x_log10(breaks = Ts, labels = Ts) +
    annotation_logticks(sides = "b", scaled = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_fixed(ratio = .25) +
    scale_color_manual(name = "Binsize", 
                       labels = c("0.25", "0.5", "1", "2"),
                       values = c(viridis::viridis(4, end = .75)),
                       guide = guide_legend(reverse = TRUE)) +
    scale_linetype_manual(name = "Parameter",
                          labels = c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a)),
                          values = c("solid", "dashed", "dotted")) +
    scale_shape_discrete(name = "Parameter",
                         labels = c("baseline"=expression(eta), "reprod"=expression(mu), "position"=expression(a))) +
    xlab("T") +
    ylab("Mean Square Error") +
    # labs(title = TeX("TO BE FILLED")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size=20))

ggsave("pareto1_convergence.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("pareto1_convergence.png", width=13.6, height=7.6, units="in")

# CHECK SPECTRAL DENSITY FUNCTION FOR PARETO3
T = 1e3
binsize = 1

lambda = 1
mu = .5
rate = 1
a = 5

n = T / binsize
I = rep(0, n)
for (k in 1:1e3) {
    x = hawkes(T, fun = lambda, repr = mu, family = "pareto", lambda = 2, a = a)
    y = discrete(x, binsize=binsize)
    I = I + Mod(fft(y - mean(y)))^2 / length(y)
}
data = new(DiscreteData, y, binsize)
model = new(Pareto2)
model$attach(data)
model$param=c(lambda, mu, a)
z<-0:(n-1)/n*2*pi
plot(z, I/1e3, type="l")
lines(z, model$f1(z, 5), col=2)
