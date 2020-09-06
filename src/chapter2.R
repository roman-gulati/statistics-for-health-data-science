##################################################
# Chapter 2 examples
##################################################
library(tidyverse)
library(scales)
library(viridis)

source('grab_meps.R')

set.seed(12345)

##################################################
# Visualize statistical distributions
##################################################
distribution_plot <- function(dset, dfun, darg, name){
    theme_update(legend.position=c(0.8, 0.8))
    gg <- ggplot(dset, aes(x=x))
    for(index in 1:length(darg)){
        iarg <- darg[[index]]
        params <- names(iarg)
        params <- sapply(params,
                         switch,
                         'size'='n',
                         'prob'='p',
                         'lambda'='lambda',
                         'shape'='alpha',
                         'rate'='beta',
                         'mean'='mu',
                         'sd'='sigma',
                         'meanlog'='mu',
                         'sdlog'='sigma')
        label <- paste(paste(params, unlist(iarg), sep='=='), collapse='*\',\'~~')
        gg <- gg+stat_function(data=dset %>% mutate(label=label),
                               aes(colour=label),
                               size=0.75,
                               fun=dfun,
                               n=nrow(dset),
                               args=iarg,
                               geom='path')
        if(grepl('poisson|binomial', name))
            gg <- gg+stat_function(data=dset %>% mutate(label=label),
                                   aes(colour=label),
                                   fun=dfun,
                                   n=nrow(dset),
                                   args=iarg,
                                   geom='point')
    }
    gg <- gg+scale_colour_viridis(name='Parameters',
                                  discrete=TRUE,
                                  labels=parse_format())
    gg <- gg+labs(x='\nx', y='f(x)\n')
    print(gg)
}

##################################################
# Binomial distribution
##################################################
binomial_plot <- function(){
    dset <- tibble(x=seq(0, 20))
    darg <- list(list(size=10, prob=0.5),
                 list(size=10, prob=0.1),
                 list(size=20, prob=0.5),
                 list(size=20, prob=0.1))
    distribution_plot(dset, dbinom, darg, 'binomial')
}
binomial_plot()

##################################################
# Poisson distribution
##################################################
poisson_plot <- function(){
    dset <- tibble(x=seq(0, 20))
    darg <- list(list(lambda=0.5),
                 list(lambda=1),
                 list(lambda=5),
                 list(lambda=9))
    distribution_plot(dset, dpois, darg, name='poisson')
}
poisson_plot()

##################################################
# Gamma distribution
##################################################
gamma_plot <- function(){
    dset <- tibble(x=seq(0, 20, by=0.1))
    darg <- list(list(shape=1, rate=1),
                 list(shape=1, rate=0.5),
                 list(shape=3, rate=1),
                 list(shape=3, rate=0.5),
                 list(shape=9, rate=1),
                 list(shape=9, rate=0.5))
    distribution_plot(dset, dgamma, darg, name='gamma')
}
gamma_plot()

##################################################
# Normal distribution
##################################################
normal_plot <- function(){
    dset <- tibble(x=seq(-5, 5, by=0.01))
    darg <- list(list(mean=0, sd=0.25),
                 list(mean=0, sd=0.5),
                 list(mean=0, sd=1),
                 list(mean=0, sd=2))
    distribution_plot(dset, dnorm, darg, name='normal')
}
normal_plot()

##################################################
# Log-normal distribution
##################################################
lognormal_plot <- function(){
    dset <- tibble(x=seq(0, 5, by=0.01))
    darg <- list(list(meanlog=0, sdlog=0.25),
                 list(meanlog=0, sdlog=0.5),
                 list(meanlog=0, sdlog=1),
                 list(meanlog=0, sdlog=2))
    distribution_plot(dset, dlnorm, darg, name='lognormal')
}
lognormal_plot()

##################################################
# Load MEPS 2017 data
##################################################
dset <- grab_meps(2017)

##################################################
# Transform and select variables
##################################################
dset <- dset %>% mutate(LOGTOTEXP17=log(TOTEXP17),
                        LOGOPTEXP17=log(OPTEXP17),
                        LOGIPTEXP17=log(IPTEXP17),
                        SCALEDIPTEXP17=IPTEXP17/10000)
dset <- dset %>% select(TOTEXP17,
                        OPTEXP17,
                        IPTEXP17,
                        LOGTOTEXP17,
                        LOGOPTEXP17,
                        LOGIPTEXP17,
                        SCALEDIPTEXP17,
                        RTHLTH31,
                        OPTOTV17,
                        STRKDX)

##################################################
# Quantify expenses by perceived health status
##################################################
expenses_by_health_status_table <- function(dset){
    rset <- dset %>% group_by(RTHLTH31)
    rset <- rset %>% summarize(Count=n(), Costs=mean(TOTEXP17))
    rset <- rset %>% filter(RTHLTH31 > 0)
    rset <- rset %>% transmute('Perceived health status'=factor(RTHLTH31,
                                                                levels=seq(5),
                                                                labels=c('Excellent',
                                                                         'Very good',
                                                                         'Good',
                                                                         'Fair',
                                                                         'Poor')),
                               'Percent of sample'=percent_format(accuracy=1)(Count/sum(Count)),
                               'Total expenditures'=dollar_format(accuracy=1)(Costs))
    return(rset)
}
expenses_by_health_status_table(dset)

##################################################
# Quantify total, inpatient, and outpatient expenses
# on natural and log-transformed scales
##################################################
natural_vs_logged_table <- function(dset){
    total.natural <- dset %>% filter(TOTEXP17 > 0) %>% summarize(Costs=mean(TOTEXP17))
    outpt.natural <- dset %>% filter(OPTEXP17 > 0) %>% summarize(Costs=mean(OPTEXP17))
    inpt.natural <- dset %>% filter(IPTEXP17 > 0) %>% summarize(Costs=mean(IPTEXP17))
    total.logged <- dset %>% filter(TOTEXP17 > 0) %>% summarize(Costs=mean(LOGTOTEXP17))
    outpt.logged <- dset %>% filter(OPTEXP17 > 0) %>% summarize(Costs=mean(LOGOPTEXP17))
    inpt.logged <- dset %>% filter(IPTEXP17 > 0) %>% summarize(Costs=mean(LOGIPTEXP17))
    lset <- tibble('Expense type'=c('Total', 'Outpatient', 'Inpatient'),
                   'Total expenses'=dollar_format(accuracy=1)(unlist(c(total.natural, outpt.natural, inpt.natural))),
                   'Logged expenses'=unlist(c(total.logged, outpt.logged, inpt.logged)))
    lset <- lset %>% mutate('Retransformed expenses'=dollar_format(accuracy=1)(exp(`Logged expenses`)))
    return(lset)
}
natural_vs_logged_table(dset)

##################################################
# Visualize inpatient expenses
##################################################
inpatient_expenses_plot <- function(dset){
    iset <- dset %>% filter(SCALEDIPTEXP17 > 0)
    hdat <- hist(iset[['SCALEDIPTEXP17']],
                 breaks=c(seq(0, 10, length=50), max(iset[['SCALEDIPTEXP17']])),
                 plot=FALSE)
    hset <- with(hdat, data.frame(mids, counts))
    inpt.summ <- iset %>% summarize(Mean=mean(IPTEXP17),
                                    p25=quantile(IPTEXP17, probs=0.25),
                                    p50=quantile(IPTEXP17, probs=0.50),
                                    p75=quantile(IPTEXP17, probs=0.75))
    gg <- ggplot(hset)
    gg <- gg+geom_bar(aes(x=mids, y=counts), stat='identity')
    gg <- gg+geom_vline(aes(xintercept=inpt.summ$Mean/10000),
                        colour='orange',
                        size=1)
    gg <- gg+geom_vline(aes(xintercept=inpt.summ$p50/10000),
                        colour='skyblue',
                        size=1)
    gg <- gg+scale_x_continuous(name='\nTotal inpatient expenditures ($10,000)',
                                limits=c(0, 10.05),
                                breaks=seq(0, 10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Persons\n',
                                limits=c(0, 300),
                                breaks=seq(0, 300, by=50),
                                expand=c(0, 0))
    print(gg)
    cat('Mean expenses:\n')
    print(inpt.summ['Mean'])
    cat('Median and 1st and 3rd quartile expenses:\n')
    print(inpt.summ[c('p25', 'p50', 'p75')])
}
inpatient_expenses_plot(dset)

##################################################
# Visualize outpatient visits
##################################################
observed_outpatient_visits_plot <- function(dset){
    dset <- dset %>% filter(STRKDX == 1)
    hdat <- hist(dset[['OPTOTV17']],
                 breaks=seq(0, 100),
                 plot=FALSE)
    hset <- with(hdat, data.frame(mids, density))
    gg <- ggplot(hset)
    gg <- gg+ggtitle('Observed')
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+scale_x_continuous(name='\nTotal outpatient visits',
                                breaks=seq(0, 10),
                                limits=c(0, 10))
    gg <- gg+scale_y_continuous(name='Frequency\n',
                                limits=c(0, 1),
                                labels=percent_format(accuracy=1))
    print(gg)
}
observed_outpatient_visits_plot(dset)

simulated_outpatient_visits_plot <- function(dset, nreps=1000){
    dset <- dset %>% filter(STRKDX == 1)
    outpt.mean <- dset %>% summarize(Mean=mean(OPTOTV17))
    sset <- tibble(visits=rpois(n=nreps, lambda=unlist(outpt.mean)))
    hdat <- with(sset, hist(visits, breaks=seq(0, 100), plot=FALSE))
    hset <- with(hdat, data.frame(mids, density))
    gg <- ggplot(hset)
    gg <- gg+ggtitle('Simulated')
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+scale_x_continuous(name='\nTotal outpatient visits',
                                breaks=seq(0, 10),
                                limits=c(0, 10))
    gg <- gg+scale_y_continuous(name='Frequency\n',
                                limits=c(0, 1),
                                labels=percent_format(accuracy=1))
    outpt.var <- dset %>% summarize(Var=var(OPTOTV17))
    cat('Observed mean:\n')
    print(outpt.mean)
    cat('Observed variance:\n')
    print(outpt.var)
    cat('Simulated variance:\n')
    print(var(sset$visits))
    cat('Observed variance/Observed mean\n')
    print(outpt.var/outpt.mean)
    print(gg)
}
simulated_outpatient_visits_plot(dset)

