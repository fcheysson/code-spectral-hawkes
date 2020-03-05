# !diagnostics off

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(ggplot2)
library(hawkesbow)
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
    # scale_y_continuous(breaks = 2*0:5, labels = 2*0:5, ) + 
    theme_bw() +
    theme(text = element_text(size=20))

ggsave("japan_measles.eps", device="ps", width=13.6, height=7.6, units="in")
ggsave("japan_measles.png", width=13.6, height=7.6, units="in")

model = new(Gaussian)
data = new(DiscreteData, measles$Count, 7)
means = c(-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0)
vars = c(1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
opt = NULL
value = Inf
for (nu in means) {
    for (sigma2 in vars) {
        model$param = c(1.0, .5, nu, sigma2)
        temp = whittle(model, data, lower = c(.01, .01, -100, .01), upper = c(100, .99, 100, 100))
        if (temp$value < value) {
            opt = temp
            value = temp$value
        }
    }
}

opt$par
