##################################################
# Chapter 6 examples
##################################################
library(tidyverse)
library(scales)
library(grid)
library(viridis)
library(lmtest)
library(quantreg)

source('grab_meps.R')

set.seed(12345)

##################################################
# Load MEPS 2017 data
##################################################
grab_and_curate_meps <- function(year, minage=18){
    dset <- grab_meps(year)
    dset <- dset %>% select(Age=AGE17X,
                            Sex=SEX,
                            Race=RACETHX,
                            Diabetes=DIABDX,
                            Inpatient=IPTEXP17,
                            Outpatient=OPTEXP17,
                            Total=TOTEXP17)
    dset <- dset %>% mutate(ScaledAge=Age-minage,
                            LogInpatient=log(Inpatient),
                            LogOutpatient=log(Outpatient),
                            LogTotal=log(Total),
                            AnyInpatient=Inpatient > 0,
                            AnyOutpatient=Outpatient > 0,
                            AnyTotal=Total > 0,
                            ScaledInpatient=Inpatient/10000,
                            ScaledOutpatient=Outpatient/10000,
                            ScaledTotal=Total/10000)
    dset <- dset %>% mutate_all(.funs=function(x) replace(x, x < 0, NA))
    dset <- dset %>% mutate(Sex=factor(Sex,
                                       levels=c(1, 2),
                                       labels=c('Male', 'Female')),
                            Race=factor(Race,
                                        levels=c(2, 1, 3, 4, 5),
                                        labels=c('Non-Hispanic White',
                                                 'Hispanic',
                                                 'Non-Hispanic Black',
                                                 'Non-Hispanic Asian',
                                                 'Non-Hispanic Mixed/Other')),
                            Race=fct_collapse(Race, 'Non-Hispanic Mixed/Other'=c('Non-Hispanic Asian',
                                                                                 'Non-Hispanic Mixed/Other')),
                            Diabetes=factor(Diabetes,
                                            levels=c(2, 1),
                                            labels=c('No', 'Yes')))
    cat('Remove',
        dset %>% filter(Age < minage) %>% nrow(),
        'participants age <', minage, '\n')
    dset <- dset %>% filter(between(Age, minage, Inf))
    return(dset)
}
dset <- grab_and_curate_meps(2017)

##################################################
# Visualize inpatient, outpatient, and total costs
##################################################
cost_histogram <- function(dset, varname, xmax=5, xstep=1, ymax=25, ystep=5){
    xlabel <- switch(varname,
                     'Inpatient'='\nInpatient costs ($10,000)',
                     'Outpatient'='\nOutpatient costs ($10,000)',
                     'Total'='\nTotal costs ($10,000)')
    scaledvarname <- paste0('Scaled', varname)
    vardata <- dset[[scaledvarname]]
    posdata <- vardata[vardata > 0]
    bpoints <- c(seq(min(posdata), xmax, length=10), max(posdata))
    hdat <- hist(posdata, breaks=bpoints, plot=FALSE)
    hset <- with(hdat, data.frame(mids, counts))
    hset <- hset %>% mutate(scaledcounts=counts/1e3)
    gg <- ggplot(hset)
    gg <- gg+geom_bar(aes(x=mids, y=scaledcounts),
                      width=0.45,
                      stat='identity')
    gg <- gg+geom_bar(data=tibble(mids=0,
                                  scaledcounts=sum(vardata == 0)/1e3),
                      aes(x=mids, y=scaledcounts),
                      width=0.1,
                      fill='orange',
                      stat='identity')
    gg <- gg+scale_x_continuous(name=xlabel,
                                breaks=seq(0, xmax, by=xstep),
                                limits=c(-0.1, xmax))
    gg <- gg+scale_y_continuous(name='Participants (1,000)\n',
                                breaks=seq(0, ymax, by=ystep),
                                limits=c(0, ymax))
    print(gg)
}
cost_histogram(dset, 'Inpatient')
cost_histogram(dset, 'Outpatient')
cost_histogram(dset, 'Total')

##################################################
# Summarize inpatient, outpatient, and total costs
##################################################
summary_helper2 <- function(dset, varname, anyvarname=NULL){
    if(is.null(anyvarname))
        sset <- dset %>% summarize('No. of participants'=comma_format(accuracy=1)(nrow(dset)),
                                   'Mean costs'=dollar_format(accuracy=1)(mean(!!sym(varname))),
                                   'Median costs'=dollar_format(accuracy=1)(median(!!sym(varname))))
    else
        sset <- dset %>% summarize('No. of participants'=comma_format(accuracy=1)(nrow(dset)),
                                   'Mean costs'=dollar_format(accuracy=1)(mean(!!sym(varname))),
                                   'Median costs'=dollar_format(accuracy=1)(median(!!sym(varname))),
                                   'Percent with zero costs'=percent_format(accuracy=1)(1-mean(!!sym(anyvarname))))
    sset <- data.frame(Group=ifelse(is.null(anyvarname), 'Non-zero costs', 'Overall'),
                       Statistic=names(sset),
                       t(sset),
                       row.names=NULL,
                       check.names=FALSE,
                       stringsAsFactors=FALSE)
    sset <- sset %>% rename(!!sym(varname):='t(sset)')
    return(sset)
}

summary_helper <- function(dset, varname){
    anyvarname <- paste0('Any', varname)
    sset <- dset %>% summary_helper2(varname, anyvarname)
    gset <- dset %>% filter(!!sym(anyvarname))
    gset <- gset %>% summary_helper2(varname)
    sset <- bind_rows(sset, gset)
    return(sset)
}

summary_statistics <- function(dset){
    iset <- summary_helper(dset, 'Inpatient')
    oset <- summary_helper(dset, 'Outpatient')
    tset <- summary_helper(dset, 'Total')
    sset <- full_join(iset, oset, by=c('Group', 'Statistic'))
    sset <- full_join(sset, tset, by=c('Group', 'Statistic'))
    return(sset)
}
summary_statistics(dset)

##################################################
# Illustrate two (log)normal distributions
##################################################
normal_density <- function(darg, qcutoff=0.2, rsize=10000, log=FALSE){
    if(log){
        dfun <- dlnorm
        qfun <- qlnorm
        rfun <- rlnorm
        xmin <- 0
        xmax <- 15
        ymax <- 1.5
        ystep <- 0.3
        ptitle <- 'After transformation'
    } else {
        dfun <- dnorm
        qfun <- qnorm
        rfun <- rnorm
        xmin <- -10
        xmax <- 10
        ymax <- 0.4
        ystep <- 0.1
        ptitle <- 'Before transformation'
    }
    tails <- lapply(darg, function(x) qfun(qcutoff, x$mean, x$sd, lower.tail=FALSE))
    means <- lapply(darg, function(x) mean(rfun(rsize, x$mean, x$sd)))
    dset <- tibble(x=seq(xmin, xmax, by=0.01))
    theme_update(legend.position=c(0.8, 0.8))
    gg <- ggplot(dset, aes(x=x))
    gg <- gg+ggtitle(ptitle)
    for(i in 1:length(darg)){
        iarg <- darg[[i]]
        ilab <- paste0('N(0,', iarg$sd^2, ')')
        iset <- dset %>% mutate(Label=ilab)
        iset <- iset %>% mutate(Mean=means[[i]])
        gg <- gg+geom_path(data=iset,
                           aes(colour=Label),
                           size=0.75,
                           stat='function',
                           fun=dfun,
                           n=nrow(dset),
                           args=iarg)
        gg <- gg+geom_area(data=iset,
                           stat='function',
                           xlim=c(tails[[i]], xmax),
                           aes(fill=Label),
                           alpha=0.3,
                           fun=dfun,
                           args=iarg,
                           show.legend=FALSE)
        gg <- gg+geom_vline(data=iset,
                            aes(xintercept=Mean,
                                colour=Label),
                            alpha=0.5,
                            size=1,
                            show.legend=FALSE)
    }
    gg <- gg+scale_colour_viridis(name='',
                                  discrete=TRUE,
                                  begin=0.2,
                                  end=0.7)
    gg <- gg+scale_fill_viridis(name='',
                                discrete=TRUE,
                                begin=0.2,
                                end=0.7)
    gg <- gg+scale_y_continuous(name='f(x)\n',
                                limits=c(0, ymax),
                                breaks=seq(0, ymax, by=ystep))
    gg <- gg+labs(x='\nx')
    print(gg)
}
darg <- list(list(mean=0, sd=1), list(mean=0, sd=2))
normal_density(darg, log=FALSE)
normal_density(darg, log=TRUE)

##################################################
# Illustrate linear regression fit to log total costs
##################################################
format_model <- function(fit, dset){
    cset <- coef(summary(fit))
    cset <- data.frame(Variable=rownames(cset), cset, check.names=FALSE)
    pname <- grep('^Pr[(][>][|]', names(cset), value=TRUE)
    cset <- cset %>% transmute(Variable=sub('(Intercept)',
                                             'Intercept',
                                             Variable,
                                             fixed=TRUE),
                               Variable=gsub(paste0('(',
                                                     paste(names(dset), collapse='|'),
                                                     ')'),
                                              '\\1 = ',
                                              Variable),
                               Variable=sub('ScaledAge', 'Age$-$18', Variable),
                               Variable=sub(' = $', '', Variable),
                               Coefficient=Estimate,
                               SE=`Std. Error`,
                               'P-value'=paste0('$',
                                                !!sym(pname),
                                                '$'))
    return(cset)
}

lognormal_regression_table <- function(dset, adjust=FALSE){
    dset <- dset %>% filter(AnyTotal)
    dset <- dset %>% filter(!is.na(Diabetes))
    if(adjust)
        fit <- lm(LogTotal~ScaledAge+Sex+Race+Diabetes, data=dset)
    else
        fit <- lm(LogTotal~Diabetes, data=dset)
    fset <- format_model(fit, dset)
    dset <- dset %>% mutate(Residuals=residuals(fit))
    gset <- dset %>% group_by(Diabetes)
    gset <- gset %>% summarize(Variance=var(Residuals))
    vset <- bind_rows(tibble(Diabetes='Overall', Variance=var(dset$Residuals)), gset)
    print(vset)
    variance_diff <- gset %>% with(diff(Variance))
    diabetes_effect <- exp(coef(fit)['DiabetesYes']+variance_diff/2)
    cat('Estimated effect of prior diabetes diagnosis:', diabetes_effect, '\n')
    print(bptest(fit))
}
lognormal_regression_table(dset, adjust=FALSE)
lognormal_regression_table(dset, adjust=TRUE)

##################################################
# Visualize gamma distributions
##################################################
gamma_densities <- function(dset, dfun, darg){
    dset <- tibble(x=seq(0, 20, by=0.1))
    darg <- list(list(shape=1, rate=1),
                 list(shape=1, rate=0.5),
                 list(shape=3, rate=1),
                 list(shape=3, rate=0.5),
                 list(shape=9, rate=1),
                 list(shape=9, rate=0.5))
    theme_update(legend.position=c(0.8, 0.8),
                 legend.text=element_text(hjust=0))
    gg <- ggplot(dset, aes(x=x))
    for(index in 1:length(darg)){
        iarg <- darg[[index]]
        params <- names(iarg)
        params <- sapply(params,
                         switch,
                         'shape'='alpha',
                         'rate'='beta')
        label <- paste(paste(params, unlist(iarg), sep='=='), collapse='*\',\'~~')
        gg <- gg+stat_function(data=dset %>% mutate(label=label),
                               aes(colour=label),
                               size=0.75,
                               fun=dgamma,
                               n=nrow(dset),
                               args=iarg,
                               geom='path')
    }
    gg <- gg+scale_colour_viridis(name='Parameters',
                                  discrete=TRUE,
                                  labels=parse_format())
    gg <- gg+labs(x='\nx', y='f(x)\n')
    print(gg)
}
gamma_densities()

##################################################
# Compare lognormal and gamma regression models
##################################################
lognormal_gamma_regression_table <- function(dset){
    dset <- dset %>% filter(AnyTotal)
    dset <- dset %>% filter(!is.na(Diabetes))
    gfit <- glm(Total~ScaledAge+Sex+Race+Diabetes,
                family=Gamma(link='log'),
                data=dset)
    nfit <- lm(LogTotal~ScaledAge+Sex+Race+Diabetes,
               data=dset)
    gset <- format_model(gfit, dset)
    nset <- format_model(nfit, dset)
    mset <- bind_rows(data.frame(Model='Gamma',
                                 gset,
                                 stringsAsFactors=FALSE,
                                 check.names=FALSE),
                      data.frame(Model='Lognormal',
                                 nset,
                                 stringsAsFactors=FALSE,
                                 check.names=FALSE))
    return(mset)
}
lognormal_gamma_regression_table(dset)

##################################################
# Visualize mixture of normal densities
##################################################
mixture_normal_density <- function(probability, size=1e4){
    size1 <- size*probability
    size2 <- size*(1-probability)
    sample1 <- rnorm(size1, mean=0, sd=1)
    sample2 <- rnorm(size2, mean=5, sd=3)
    label <- paste0('p == ', probability)
    dset <- tibble(x=c(sample1, sample2),
                   label=label)
    theme_update(legend.position=c(0.8, 0.8))
    gg <- ggplot(dset)
    gg <- gg+geom_density(aes(x=x, fill=label),
                          colour=NA,
                          adjust=2,
                          size=0.75)
    gg <- gg+geom_vline(aes(xintercept=mean(x)), size=0.75)
    gg <- gg+scale_x_continuous(name='\nx',
                                breaks=seq(-10, 20, by=5),
                                limits=c(-11, 21))
    gg <- gg+scale_y_continuous(name='f(x)\n',
                                breaks=seq(0, 0.3, by=0.1),
                                limits=c(0, 0.3))
    gg <- gg+scale_fill_viridis(name='',
                                discrete=TRUE,
                                begin=0.45,
                                labels=parse_format())
    gg <- gg+labs(y='f(x)\n')
    print(gg)
}
mixture_normal_density(probability=0.2)
mixture_normal_density(probability=0.7)

##################################################
# Two-part model of total medical expenditures
##################################################
two_part_regression_table <- function(dset){
    dset <- dset %>% filter(!is.na(Diabetes))
    bfit <- glm(AnyTotal~ScaledAge+Sex+Race+Diabetes,
                family=binomial(link='logit'),
                data=dset)
    gfit <- glm(Total~ScaledAge+Sex+Race+Diabetes,
                family=Gamma(link='log'),
                data=dset,
                subset=AnyTotal)
    bset <- format_model(bfit, dset)
    gset <- format_model(gfit, dset)
    mset <- bind_rows(data.frame(Model='Logistic',
                                 bset,
                                 stringsAsFactors=FALSE,
                                 check.names=FALSE),
                      data.frame(Model='Gamma',
                                 gset,
                                 stringsAsFactors=FALSE,
                                 check.names=FALSE))
    bpred0 <- predict(bfit, newdata=dset %>% mutate(Diabetes='No'), type='response')
    bpred1 <- predict(bfit, newdata=dset %>% mutate(Diabetes='Yes'), type='response')
    gpred0 <- predict(gfit, newdata=dset %>% mutate(Diabetes='No'), type='response')
    gpred1 <- predict(gfit, newdata=dset %>% mutate(Diabetes='Yes'), type='response')
    mae <- mean(bpred1*gpred1-bpred0*gpred0)
    mme <- mean(bpred1*gpred1/(bpred0*gpred0))
    cat('Marginal additive effect:', round(mae), '\n')
    cat('Marginal multiplicative effect:', round(mme, 3), '\n')
    lfit <- lm(Total~ScaledAge+Sex+Race+Diabetes, data=dset)
    cat('Linear regression estimate:', round(coef(lfit)['DiabetesYes']), '\n')
}
two_part_regression_table(dset)

###########################################################
# Median regression for costs and log(costs)
# Note: solution is nonunique because of binary covariate
###########################################################
format_quantile_model <- function(fit, dset, group){
    cset <- coef(summary(fit))
    cset <- data.frame(Variable=rownames(cset), cset, check.names=FALSE, row.names=NULL)
    cset <- cset %>% transmute(Group=group,
                               Variable=sub('(Intercept)', 'Intercept', Variable, fixed=TRUE),
                               Variable=gsub(paste0('(', paste(names(dset), collapse='|'), ')'), '\\1 = ', Variable),
                               Variable=sub('ScaledAge', 'Age$-$18', Variable),
                               Variable=sub(' = $', '', Variable),
                               Coefficient=Value,
                               SE=`Std. Error`,
                               'P-value'=paste0('$', `Pr(>|t|)`, '$'))
    return(cset)
}

median_regression_analysis <- function(dset){
    aset <- dset %>% filter(!is.na(Diabetes))
    pset <- aset %>% filter(AnyTotal)
    afit <- rq(Total~ScaledAge+Sex+Race+Diabetes, data=aset, tau=0.5)
    amod <- format_quantile_model(afit, aset, 'Overall')
    pfit <- rq(Total~ScaledAge+Sex+Race+Diabetes, data=pset, tau=0.5)
    pmod <- format_quantile_model(pfit, aset, 'Non-zero costs')
    mset <- bind_rows(amod, pmod)
    return(mset)
}
median_regression_analysis(dset)

