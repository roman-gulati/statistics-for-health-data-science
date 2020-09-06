##################################################
# Create figures and tables for Chapter 2
##################################################
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(viridis)
library(xtable)

set.seed(12345)

##################################################
# Load MEPS 2017 data
##################################################
source('shared.R')
dset <- grab_meps(2017)

##################################################
# Visualize statistical distributions
##################################################
distribution_plot <- function(dset, dfun, darg, filename, ext='pdf', saveit=FALSE){
    gg_theme(legend.position=c(0.8, 0.8))
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
        if(grepl('poisson|binomial', filename))
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
    if(saveit){
        filename <- paste('02', filename, sep='-')
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
}

##################################################
# Binomial distribution
##################################################
binomial_plot <- function(saveit=FALSE){
    dset <- tibble(x=seq(0, 20))
    darg <- list(list(size=10, prob=0.5),
                 list(size=10, prob=0.1),
                 list(size=20, prob=0.5),
                 list(size=20, prob=0.1))
    distribution_plot(dset,
                      dbinom,
                      darg,
                      filename='binomial',
                      saveit=saveit)
}
#binomial_plot(saveit=TRUE)

##################################################
# Multinomial distribution (ternary plot)
##################################################
multinomial_plot <- function(ext='pdf', saveit=FALSE){
    library(ggtern)
    n <- 100
    p <- c(0.35, 0.5, 0.15)
    x1 <- seq(0, n)
    x2 <- seq(0, n)
    x3 <- seq(0, n)
    eset <- tidyr::expand_grid(x1, x2, x3)
    eset <- eset %>% filter(x1+x2+x3 == n)
    eset <- eset %>% mutate(d=apply(eset, 1, dmultinom, prob=p))
    gg_theme(legend.position='right')
    gg <- ggtern(data=eset, aes(x=x1, y=x2, z=x3))
    gg <- gg+geom_point(aes(colour=d), size=1)
    gg <- gg+scale_colour_viridis(name='',
                                  begin=0.9,
                                  end=0,
                                  limits=c(0, 0.01),
                                  expand=c(0, 0),
                                  breaks=seq(0, 0.01, by=0.002),
                                  guide=guide_colorbar(title='Density',
                                                       ticks=FALSE,
                                                       draw.ulim=FALSE,
                                                       draw.llim=FALSE))
    gg <- gg+labs(x='Type 1', y='Type 2', z='Type 3')
    print(gg)
    if(saveit){
        filename <- '02-multinomial'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
}
#multinomial_plot(saveit=TRUE)

##################################################
# Poisson distribution
##################################################
poisson_plot <- function(saveit=FALSE){
    dset <- tibble(x=seq(0, 20))
    darg <- list(list(lambda=0.5),
                 list(lambda=1),
                 list(lambda=5),
                 list(lambda=9))
    distribution_plot(dset,
                      dpois,
                      darg,
                      filename='poisson',
                      saveit=saveit)
}
#poisson_plot(saveit=TRUE)

##################################################
# Gamma distribution
##################################################
gamma_plot <- function(saveit=FALSE){
    dset <- tibble(x=seq(0, 20, by=0.1))
    darg <- list(list(shape=1, rate=1),
                 list(shape=1, rate=0.5),
                 list(shape=3, rate=1),
                 list(shape=3, rate=0.5),
                 list(shape=9, rate=1),
                 list(shape=9, rate=0.5))
    distribution_plot(dset,
                      dgamma,
                      darg,
                      filename='gamma',
                      saveit=saveit)
}
#gamma_plot(saveit=TRUE)

##################################################
# Normal distribution
##################################################
normal_plot <- function(saveit=FALSE){
    dset <- tibble(x=seq(-5, 5, by=0.01))
    darg <- list(list(mean=0, sd=0.25),
                 list(mean=0, sd=0.5),
                 list(mean=0, sd=1),
                 list(mean=0, sd=2))
    distribution_plot(dset,
                      dnorm,
                      darg,
                      filename='normal',
                      saveit=saveit)
}
#normal_plot(saveit=TRUE)

##################################################
# Log-normal distribution
##################################################
lognormal_plot <- function(saveit=FALSE){
    dset <- tibble(x=seq(0, 5, by=0.01))
    darg <- list(list(meanlog=0, sdlog=0.25),
                 list(meanlog=0, sdlog=0.5),
                 list(meanlog=0, sdlog=1),
                 list(meanlog=0, sdlog=2))
    distribution_plot(dset,
                      dlnorm,
                      darg,
                      filename='lognormal',
                      saveit=saveit)
}
#lognormal_plot(saveit=TRUE)

##################################################
# Transform variables and specify design for MEPS 2017
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
expenses_by_health_status_table <- function(dset, saveit=FALSE){
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
    if(saveit){
        filename <- '02-expenses_by_health_status'
        filename <- paste(filename, 'tex', sep='.')
        Caption <- 'Total medical expenditures per person in MEPS 2017 by perceived health status.'
        print(xtable(rset,
                     align='cccc',
                     digits=0,
                     label='condmean',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              sanitize.text.function=sanitize,
              include.rownames=FALSE,
              hline.after=0)
    }
}
#expenses_by_health_status_table(dset, saveit=TRUE)

##################################################
# Quantify total, inpatient, and outpatient expenses
# on natural and log-transformed scales
##################################################
natural_vs_logged_table <- function(dset, saveit=FALSE){
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
    if(saveit){
        filename <- '02-natural_vs_logged_expenses'
        filename <- paste(filename, 'tex', sep='.')
        Caption <- 'Total, outpatient and inpatient medical expenses per person in MEPS 2017 calculated directly, after log-transformation, and after retransformation.'
        print(xtable(lset,
                     align='ccccc',
                     digits=2,
                     label='tablelog',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              sanitize.text.function=sanitize,
              include.rownames=FALSE,
              hline.after=0)
    }
}
#natural_vs_logged_table(dset, saveit=TRUE)

##################################################
# Visualize inpatient expenses
##################################################
inpatient_expenses_plot <- function(dset, ext='pdf', saveit=FALSE){
    iset <- dset %>% filter(SCALEDIPTEXP17 > 0)
    hdat <- hist(iset[['SCALEDIPTEXP17']],
                 breaks=c(seq(0, 10, length=50), max(iset[['SCALEDIPTEXP17']])),
                 plot=FALSE)
    hset <- with(hdat, data.frame(mids, counts))
    inpt.summ <- iset %>% summarize(Mean=mean(IPTEXP17),
                                    p25=quantile(IPTEXP17, probs=0.25),
                                    p50=quantile(IPTEXP17, probs=0.50),
                                    p75=quantile(IPTEXP17, probs=0.75))
    gg_theme()
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
    if(saveit){
        filename <- '02-inpatient_expenses'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
    cat('Mean expenses:\n')
    print(inpt.summ['Mean'])
    cat('Median and 1st and 3rd quartile expenses:\n')
    print(inpt.summ[c('p25', 'p50', 'p75')])
}
#inpatient_expenses_plot(dset, saveit=TRUE)

##################################################
# Visualize outpatient visits
##################################################
observed_outpatient_visits_panel <- function(dset){
    hdat <- hist(dset[['OPTOTV17']],
                 breaks=seq(0, 100),
                 plot=FALSE)
    hset <- with(hdat, data.frame(mids, density))
    gg_theme()
    gg <- ggplot(hset)
    gg <- gg+ggtitle('Observed')
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+scale_x_continuous(name='\nTotal outpatient visits',
                                breaks=seq(0, 10),
                                limits=c(0, 10))
    gg <- gg+scale_y_continuous(name='Frequency\n',
                                limits=c(0, 1),
                                labels=percent_format(accuracy=1))
    return(gg)
}

simulated_outpatient_visits_panel <- function(dset, nreps=1000){
    outpt.mean <- dset %>% summarize(Mean=mean(OPTOTV17))
    sset <- tibble(visits=rpois(n=nreps, lambda=unlist(outpt.mean)))
    hdat <- with(sset, hist(visits, breaks=seq(0, 100), plot=FALSE))
    hset <- with(hdat, data.frame(mids, density))
    gg_theme()
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
    return(gg)
}

outpatient_visits_plot <- function(dset, ext='pdf', saveit=FALSE){
    oset <- dset %>% filter(STRKDX == 1)
    opanel <- observed_outpatient_visits_panel(oset)
    spanel <- simulated_outpatient_visits_panel(oset)
    Layout <- grid.layout(ncol=2, nrow=1, heights=unit(5, 'null'))
    if(saveit){
        filename <- '02-outpatient_visits'
        filename <- paste(filename, ext, sep='.')
        get(ext)(here('figures', filename),
                 height=5,
                 width=10)
        on.exit(graphics.off())
    }
    grid.newpage()
    pushViewport(viewport(layout=Layout))
    print(opanel, vp=viewport(layout.pos.col=1, layout.pos.row=1))
    print(spanel, vp=viewport(layout.pos.col=2, layout.pos.row=1))
}
#outpatient_visits_plot(dset, saveit=TRUE)

