##################################################
# Chapter 1 examples
##################################################
library(tidyverse)
source('grab_meps.R')

##################################################
# Visualize histogram of total expenditures
##################################################
totexp_panel <- function(year, varname){
    dset <- grab_meps(year)
    dset <- dset %>% mutate(SCALEDTOTEXP17=TOTEXP17/10000,
                            LOGTOTEXP17=log(TOTEXP17+1))
    dset <- dset %>% select(SCALEDTOTEXP17, LOGTOTEXP17)
    if(varname == 'SCALEDTOTEXP17'){
        xmax <- 10
        xstep <- 1
        ymax <- 35
        ystep <- 5
        xlabel <- '\nTotal expenditure ($10,000)'
        bpoints <- c(seq(0, xmax, length=15), max(dset[['SCALEDTOTEXP17']]))
    } else {
        xmax <- 16
        xstep <- 2
        ymax <- 6
        ystep <- 1
        xlabel <- '\nLog(Total expenditure+$1)'
        bpoints <- 20
    }
    hdat <- hist(dset[[varname]], breaks=bpoints, plot=FALSE)
    hset <- with(hdat, data.frame(mids, counts))
    hset <- hset %>% mutate(scaledcounts=counts/1000)
    gg <- ggplot(hset)
    gg <- gg+geom_bar(aes(x=mids, y=scaledcounts), stat='identity')
    gg <- gg+scale_x_continuous(name=xlabel,
                                breaks=seq(0, xmax, by=xstep),
                                limits=c(0, xmax))
    gg <- gg+scale_y_continuous(name='Participants (1,000)\n',
                                breaks=seq(0, ymax, by=ystep),
                                limits=c(0, ymax))
    print(gg)
}
totexp_panel(2017, 'SCALEDTOTEXP17')
totexp_panel(2017, 'LOGTOTEXP17')

