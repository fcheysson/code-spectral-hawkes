# !diagnostics off

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(ggplot2)
library(hawkesbow)
library(latex2exp)
library(lubridate)

################ Example links for Japan weekly datasets #####################
# https://www.niid.go.jp/niid/images/idwr/data-e/idwr-e2020/202008/zensu08.csv
# https://www.niid.go.jp/niid/images/idwr/data-e/idwr-e2015/201501/zensu01.csv
# https://www.niid.go.jp/niid/images/idwr/data-e/idwr-e2014/1452/zensu52.csv
# https://www.niid.go.jp/niid/images/idwr/data-e/idwr-e2012/1238/zensu38.csv
# https://www.niid.go.jp/niid/images/idwr/data-e/1237/zensu37.csv
# https://www.niid.go.jp/niid/images/idwr/data-e/1231/zensu31.csv
##############################################################################

format2Digits = function(x) formatC(x, width = 2, format = "d", flag = "0")

# number of weeks: 52*7 (2013-2019) + 1 (2015 has 53 weeks) + 8 (2020) + 22 (2012)
n = 52*7 + 1 + 8 + 22

year = 2020
week = 08 %>% format2Digits()
prefecture = "Tokyo"
measles = tibble(Year = integer(n), Week = integer(n), Count = integer(n))
it = 1
while (TRUE) { # It will break when an URL is no longer readable
    
    # Display notification
    if ((year != 2015 && week == "52") || (year == 2015 && week == "53"))
        cat("\n", year, ": ", sep="")
    cat(week, ". ", sep="")
    
    # URL link
    if (year >= 2015)
        url = paste0("https://www.niid.go.jp/niid/images/idwr/data-e/idwr-e", year, "/", year, week, "/zensu", week, ".csv")
    else if (year >=2013 || week >= 38)
        url = paste0("https://www.niid.go.jp/niid/images/idwr/data-e/idwr-e", year, "/", year-2000, week, "/zensu", week, ".csv")
    else 
        url = paste0("https://www.niid.go.jp/niid/images/idwr/data-e/", year-2000, week, "/zensu", week, ".csv")
    
    # Read what's needed
    connection = url(url)
    df = read.csv(connection, header = FALSE, stringsAsFactors = FALSE, sep = ",")
    prefectRow = grep(prefecture, df[,1])
    measlesCol = grep("Measles", df[4,])
    
    # Save what's needed
    measles$Year[it] = year
    measles$Week[it] = week %>% as.integer()
    measles$Count[it] = df[prefectRow, measlesCol] %>% as.integer()
    
    # Reduce week by one
    if (week != "01") {
        week = (as.numeric(week) - 1) %>% format2Digits()
    } else {
        year = year - 1
        week = ifelse(year == "2015", 53, 52) %>% format2Digits()
    }
    
    it = it + 1
    
}

measles = measles %>% drop_na() %>% filter(Year != 0)
measles = measles %>% 
    mutate(Date = ISOweek::ISOweek2date(paste0(Year, "-W", Week %>% format2Digits(), "-4"))) %>% 
    select(Date, everything()) %>% 
    arrange(Date)

ggplot(measles, aes(x=Date, y=Count)) +
    geom_segment(aes(xend=Date, yend=0)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) + 
    theme_bw() +
    theme(text = element_text(size=20))

ggsave("japan_measles.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("japan_measles.png", width=13.6, height=7.6, units="in")

model = new(Gaussian)
model$binsize = 7
ymean = mean(measles$Count)
reprs = c(.25, .5, .75)
means = c(-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0)
vars = c(1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
opt = NULL
value = Inf
for (mu in reprs) {
    for (nu in means) {
        for (sigma2 in vars) {
            init = c(ymean * (1 - mu), mu, nu, sigma2)
            temp = whittle(measles$Count, model, init = init, lower = c(.01, .01, -100, .01), upper = c(100, .99, 100, 100))
            if (temp$opt$value < value) {
                opt = temp
                value = temp$opt$value
            }
        }
    }
}

opt$opt$par

## Goodness-of-fit testing
library(foreach)
library(doParallel)
library(hawkesbow)
library(SuperGauss)
source("whittle_periodogram.R")

y = measles$Count
n = length(y)
Nsup = ceiling((n-1)/2)
Ninf = - floor((n-1)/2)
omega = 2 * pi * Ninf:Nsup / n

kern = function(x) ifelse(abs(x) <= pi, 1.5 * (1.0 - (x / pi)^2), 0.0)
kernh = function(x, h) kern(x / h)
h = 0.10

k = kernh(2* pi * 0:(n-1) / n, h)
q = function(J) 1 / (n*h) * SuperGauss::toep.mult(k, J)

I = Mod(fft(y - mean(y)))^2 / n
I = c(I[(Nsup+2):n], I[1:(Nsup+1)])
f = opt$kernel$f1(omega, trunc = 10L)
S = 2*pi*sqrt(h) * sum( q(I/f - 1)^2 )

init = opt$par

n.sim = 1e3
cl = makePSOCKcluster(3) # use multicore, set to the number of our cores
clusterExport(cl, c("n", "Nsup", "omega", "h", "k", "f", "q", "init"))
registerDoParallel(cl)
Ss = foreach (i = icount(n.sim), .combine = c, .multicombine = TRUE, .inorder = FALSE,
              .packages = c("SuperGauss", "hawkesbow")) %dopar% {
                  Iboots = f * rexp(n)
                  Iboots[n - Nsup] = f[n - Nsup] * rchisq(1, df = 1)
                  model = whittle_periodogram(c(Iboots[(n-Nsup+1):n], Iboots[1:(n-Nsup-1)]), 
                                              kern = "gauss", binsize = 7, init = init)
                  fboots = model$f1(omega, trunc = 10L)
                  return( 2*pi*sqrt(h) * sum( q(Iboots/fboots - 1)^2 ) )
              }
stopCluster(cl)

muh = 1/sqrt(h) * 12*pi/5
tau2 = 2672*pi^2/385
# plot(Ss)
plot(z<-0:100, dnorm(z, mean = muh, sd = sqrt(tau2)), col=2, type="l")
lines(density(Ss))
abline(v=S, col="grey")
cat("mean(Ss) = ", mean(Ss), " and     sd(Ss) = ", sd(Ss),
    "\n     muh = ", muh,      " and sqrt(tau2) = ", sqrt(tau2), 
    sep = "")
c(muh - 1.96*sqrt(tau2), muh + 1.96*sqrt(tau2))

omega = 2 * pi * 1:Nsup / n
f = opt$kernel$f1(omega, trunc = 10L)
k = kernh(2* pi * 0:(Nsup-1) / n, h)
q = function(J) 1 / (n*h) * SuperGauss::toep.mult(k, J)
I = (Mod(fft(y - mean(y)))^2 / n)
I = I[2:(Nsup+1)]
Q = 2*pi*n*h * q(I/f - 1)^2 / (12*pi/5)
plot(omega, I/f, type="b")
lines(omega, q(I/f), type="l", col=2)
lines(omega, q(I/f) / q(rep(1, Nsup)), type="l", col=2)    # Renormalised
plot(omega, Q, type="l")
abline(h=qchisq(.95, 1))

# Plots
omega = 2 * pi * 1:Nsup / n
f = opt$kernel$f1(omega, trunc = 10L)
I = (Mod(fft(y - mean(y)))^2 / n)
I = I[2:(Nsup+1)]
qtib = tibble(freq = omega, value = I/f, type = "I/f")
h = 0.05
k = kernh(2* pi * 0:(Nsup-1) / n, h)
q = function(J) 1 / (n*h) * SuperGauss::toep.mult(k, J)
qtib = bind_rows(qtib, 
                 tibble(freq = omega, value = q(I/f) / q(rep(1, Nsup)), type = "h = 0.05"))
h = 0.10
k = kernh(2* pi * 0:(Nsup-1) / n, h)
q = function(J) 1 / (n*h) * SuperGauss::toep.mult(k, J)
qtib = bind_rows(qtib, 
                 tibble(freq = omega, value = q(I/f) / q(rep(1, Nsup)), type = "h = 0.10"))

ggplot(mapping = aes(x = freq, y = value)) +
    geom_line(data=filter(qtib, type == "I/f"), col="grey70") +
    geom_point(data=filter(qtib, type == "I/f"), col="grey70") +
    geom_line(aes(group = type, linetype = type), data = filter(qtib, type != "I/f")) +
    scale_x_continuous(breaks=seq(0, pi, by=pi/4), labels=sapply(c("$0$", "$\\pi/4$", "$\\pi/2$", "$3\\pi/4$", "$\\pi$"), TeX), limits=pi*c(0,1)) +
    scale_linetype_manual(name = "Bandwidth",
                          values = c("dashed", "solid")) + 
    xlab("Frequency") +
    ylab("") +
    guides(linetype = FALSE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text = element_text(size=20))

ggsave("goodnessOfFit1.eps", device="ps", width=13.6, height=7.6, units="in")

# Plot 2
omega = 2 * pi * 1:Nsup / n
f = opt$kernel$f1(omega, trunc = 10L)
h = 0.10
k = kernh(2* pi * 0:(Nsup-1) / n, h)
q = function(J) 1 / (n*h) * SuperGauss::toep.mult(k, J)
I = (Mod(fft(y - mean(y)))^2 / n)
I = I[2:(Nsup+1)]
Q = 2*pi*n*h * q(I/f - 1)^2 / (12*pi/5)
Qtib = tibble(freq = omega, value = Q)

ggplot(data = Qtib, mapping = aes(x = freq, y = value)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks=seq(0, pi, by=pi/4), labels=sapply(c("$0$", "$\\pi/4$", "$\\pi/2$", "$3\\pi/4$", "$\\pi$"), TeX), limits=pi*c(0,1)) +
    xlab("Frequency") +
    ylab("") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text = element_text(size=20))

ggsave("goodnessOfFit2.eps", device="ps", width=13.6, height=7.6, units="in")

## WEEKPLOT FOR MEASLES

ggplot(mutate(measles, Week = Week) %>% filter(Week <= 52), aes(x = Week, group = Week, y = Count)) +
    geom_boxplot() +
    stat_smooth(aes(group = 1), method = "loess", col = "blue", span = 0.1, se = FALSE) +
    scale_x_continuous(breaks = seq(0, 52, by = 13))

ggsave("tokyoMeasles_weekplot.png")
