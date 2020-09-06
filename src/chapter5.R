##################################################
# Chapter 5 examples
##################################################
library(tidyverse)
library(grid)
library(scales)
library(viridis)
library(foreign)
library(MASS)
library(pscl)
library(tableone)

##################################################
# Basic Health Plan data
##################################################
grab_and_curate_bph <- function(filename){
    dset <- read.dta(filename, convert.factors=FALSE)
    dset <- dset %>% filter(between(age, 18, Inf))
    dset <- dset %>% mutate(sex=factor(sex, levels=c(0, 1), labels=c('Women', 'Men')),
                            raceind=factor(raceind, levels=c(0, 1), labels=c('White', 'Other')),
                            nchroniccat=factor(nchronic),
                            nchroniccat=fct_collapse(nchroniccat, '4'=as.character(seq(4, max(nchronic, na.rm=TRUE)))),
                            provider=factor(provider, labels=c('Spokane HMO',
                                                               'Clallam',
                                                               'Pierce',
                                                               'Spokane IPA')))
    return(dset)
}
dset <- grab_and_curate_bph('basic_health_plan.dta')

##################################################
# Prostate cancer mortality data
##################################################
grab_and_curate_seer <- function(filename){
    dset <- read.dta(filename, convert.factors=FALSE)
    dset <- dset %>% rename(Race='race',
                            Year='yeardeath',
                            Rate='deathrate',
                            Deaths='deathtotal',
                            Population='population')
    dset <- dset %>% mutate(Race=factor(Race,
                                        levels=c('white', 'black'),
                                        labels=c('White', 'Black')))
    dset <- dset %>% dplyr::select(Race, Year, Deaths, Population, Rate)
    return(dset)
}
sset <- grab_and_curate_seer('seer_prostate_cancer_mortality.dta')

##################################################
# Visualize Poisson and negative binomial distributions
# with the same mean
##################################################
process_params <- function(iarg){
    params <- names(iarg)
    params <- sapply(params,
                     switch,
                     'mu'='mu',
                     'size'='alpha',
                     'lambda'='mu')
    return(paste(paste(params, unlist(iarg), sep='=='), collapse='*\',\'~~'))
}

scale_linear <- function(aarg, aargs){
    amin <- min(aargs)
    amax <- max(aargs)
    return((aarg-amin)/(amax-amin))
}

poisson_negative_binomial_plot <- function(dset=tibble(x=seq(0, 20)),
                                           dfun1=dpois,
                                           dfun2=dnbinom,
                                           darg=list(list(lambda=1),
                                                     list(lambda=5),
                                                     list(lambda=9)),
                                           aargs=exp(seq(0, 3, by=0.5))){
    theme_update(legend.position=c(0.8, 0.8))
    gg <- ggplot(dset, aes(x=x))
    for(index in 1:length(darg)){
        iarg <- darg[[index]]
        label <- process_params(iarg)
        gg <- gg+stat_function(data=dset %>% mutate(label=label),
                               aes(colour=label),
                               linetype='solid',
                               size=0.75,
                               fun=dfun1,
                               n=nrow(dset),
                               args=iarg,
                               geom='path')
        gg <- gg+stat_function(data=dset %>% mutate(label=label),
                               aes(colour=label),
                               fun=dfun1,
                               n=nrow(dset),
                               args=iarg,
                               geom='point')
        if(!is.null(dfun2)){
            iarg$mu <- iarg$lambda
            for(aarg in aargs){
                iarg$size <- aarg
                iarg$lambda <- NULL
                gg <- gg+stat_function(data=dset %>% mutate(label=label),
                                       aes(colour=label),
                                       alpha=scale_linear(aarg, aargs),
                                       size=0.75,
                                       fun=dfun2,
                                       n=nrow(dset),
                                       args=iarg,
                                       geom='path')
                gg <- gg+stat_function(data=dset %>% mutate(label=label),
                                       aes(colour=label),
                                       alpha=scale_linear(aarg, aargs),
                                       fun=dfun2,
                                       n=nrow(dset),
                                       args=iarg,
                                       geom='point')
            }
        }
    }
    gg <- gg+scale_colour_viridis(name='Parameters',
                                  discrete=TRUE,
                                  labels=parse_format())
    gg <- gg+labs(x='\nx', y='f(x)\n')
    print(gg)
}
poisson_negative_binomial_plot(dfun2=NULL)
poisson_negative_binomial_plot(dfun2=dnbinom)

##################################################
# Descriptive tables
##################################################
bhp_table <- function(dset){
    dset <- dset %>% mutate(provider=factor(provider, levels=c('Spokane HMO',
                                                               'Spokane IPA',
                                                               'Pierce',
                                                               'Clallam')),
                            raceind=fct_explicit_na(raceind, na_level='Unknown'),
                            nchroniccat=fct_explicit_na(nchroniccat, na_level='Unknown'))
    dset <- dset %>% rename(Age='age',
                            Enrollment='bhpmonth',
                            Sex='sex',
                            Race='raceind',
                            Conditions='nchroniccat',
                            Provider='provider',
                            Visits='outvis')
    table1 <- CreateTableOne(vars=c('Visits',
                                    'Enrollment',
                                    'Age',
                                    'Sex',
                                    'Race',
                                    'Conditions'),
                             strata='Provider',
                             data=dset,
                             factorVars=c('Sex',
                                          'Race',
                                          'Conditions'))
    tset <- print(table1, printToggle=FALSE)
    tset <- subset(tset, select=-c(p, test))
    tset <- data.frame(Characteristic=rownames(tset),
                       tset,
                       row.names=NULL,
                       check.names=FALSE)
    tset <- tset %>% mutate(Characteristic=sub('%',
                                               '\\%',
                                               Characteristic,
                                               fixed=TRUE),
                            Characteristic=sub('  ',
                                               '\\quad',
                                               Characteristic,
                                               fixed=TRUE))
    return(tset)
}
bhp_table(dset)

seer_table_helper <- function(dset, race, ivars=c('Year', 'Race')){
    dset <- dset %>% filter(Race == race)
    dset <- dset %>% mutate(Deaths=label_comma(accuracy=1)(Deaths),
                            Population=label_comma(accuracy=1)(Population))
    mvars <- setdiff(names(dset), ivars)
    names(dset)[match(mvars, names(dset))] <- paste0(mvars, ' (', unique(dset$Race), ')')
    dset <- dset %>% dplyr::select(-Race)
    return(dset)
}

seer_table <- function(dset){
    dset <- dset %>% mutate(Race=substr(Race, 0, 1))
    wset <- dset %>% seer_table_helper('W')
    bset <- dset %>% seer_table_helper('B')
    xset <- full_join(wset, bset, by='Year')
    return(xset)
}
seer_table(sset)

##################################################
# Visualize outpatient visits
##################################################
bhp_plot <- function(dset){
    fset <- dset %>% dplyr::select(outvis, provider)
    fset <- fset %>% group_by(provider)
    fset <- fset %>% summarize(Mean=mean(outvis),
                               SD=sd(outvis))
    gg <- ggplot(fset, aes(x=Mean, y=SD))
    gg <- gg+geom_point(aes(colour=provider), size=4)
    gg <- gg+geom_text(aes(colour=provider, label=provider),
                       vjust=2,
                       size=6)
    gg <- gg+scale_x_continuous(name='\nMean number of outpatient visits',
                                limits=c(0, 7),
                                breaks=seq(0, 7),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Standard deviation of number of outpatient visits\n',
                                limits=c(0, 12),
                                breaks=seq(0, 12, by=2),
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(discrete=TRUE, begin=0.6, end=0.1)
    print(gg)
}
bhp_plot(dset)

##################################################
# Compare histogram of outpatient visits to Poisson
# distribution with the same mean
##################################################
bhp_histogram <- function(dset){
    bpoints <- c(seq(-0.5, 40.5), 500)
    hdat <- hist(dset[['outvis']], breaks=bpoints, plot=FALSE)
    hset <- with(hdat, data.frame(mids, density))
    hset <- hset %>% mutate(poisson=dpois(mids, lambda=mean(dset$outvis)))
    gg <- ggplot(hset)
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+geom_path(aes(x=mids, y=poisson),
                       colour='tomato',
                       size=0.75)
    gg <- gg+geom_point(aes(x=mids, y=poisson),
                        colour='tomato')
    gg <- gg+scale_x_continuous(name='\nNumber of outpatient visits',
                                limits=c(-1, 40),
                                breaks=seq(0, 40, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Frequency\n',
                                limits=c(0, 0.25),
                                breaks=seq(0, 0.25, by=0.05),
                                expand=c(0, 0))
    print(gg)
}
bhp_histogram(dset)

##################################################
# Fit Poisson and negative binomial models
##################################################
format_model_helper <- function(cset, dset){
    cset <- data.frame(Predictor=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% transmute(Predictor=sub('(Intercept)', 'Intercept', Predictor, fixed=TRUE),
                               Predictor=gsub(paste0('(', paste(names(dset), collapse='|'), ')'), '\\1 = ', Predictor),
                               Predictor=sub('age =', 'age', Predictor),
                               Coefficient=sub('-', '$-$', sprintf('%4.2f', Estimate), fixed=TRUE),
                               'Relative risk'=sprintf('%4.2f', exp(Estimate)),
                               '95% CI'=paste0('(',
                                               sprintf('%4.2f', exp(Estimate-1.96*`Std. Error`)),
                                               ', ',
                                               sprintf('%4.2f', exp(Estimate+1.96*`Std. Error`)),
                                               ')'),
                               'P-value'=sub('(.*)', '$\\1 $', `Pr(>|z|)`))
    cset <- cset %>% filter(grepl('^provider', Predictor))
    cset <- cset %>% mutate(Predictor=sub('^provider = ', '', Predictor))
    cset <- cset %>% rename(Provider='Predictor')
    return(cset)
}

format_model <- function(fit, dset){
    cset <- coef(summary(fit))
    if(inherits(fit, 'zeroinfl')){
        zset <- format_model_helper(cset$zero, dset)
        cset <- format_model_helper(cset$count, dset)
        cset <- bind_rows(data.frame(Component='Zero model',
                                     zset,
                                     stringsAsFactors=FALSE,
                                     check.names=FALSE),
                          data.frame(Component='Count model',
                                     cset,
                                     stringsAsFactors=FALSE,
                                     check.names=FALSE))
        cset <- cset %>% dplyr::select(-Coefficient)
        cset <- cset %>% mutate(Component=factor(Component, levels=c('Zero model',
                                                                     'Count model')))
        cset <- cset %>% arrange(Component)
    } else {
        cset <- format_model_helper(cset, dset)
    }
    cset <- cset %>% mutate(Provider=factor(Provider, levels=c('Spokane IPA',
                                                               'Pierce',
                                                               'Clallam')))
    cset <- cset %>% arrange(Provider)
    return(cset)
}

bhp_model_table <- function(dset,
                            distribution='Poisson',
                            ext='pdf',
                            exposure=FALSE){
    response <- 'outvis'
    predictor <- c('raceind', 'sex', 'age', 'provider', 'nchroniccat')
    if(exposure & !grepl('zero', distribution)){
        mset <- dset %>% dplyr::select(all_of(c(response, predictor, 'bhpmonth')))
        predictor <- c(predictor, 'offset(log(bhpmonth))')
    }
    predictor <- paste(predictor, collapse='+')
    equation <- paste(response, predictor, sep='~')
    if(distribution == 'Poisson')
        fit <- glm(as.formula(equation),
                   family=poisson(link='log'),
                   data=dset)
    if(distribution == 'negative binomial')
        fit <- glm.nb(as.formula(equation),
                      data=dset)
    if(distribution == 'zero-inflated Poisson')
        fit <- zeroinfl(as.formula(equation),
                        offset=log(bhpmonth),
                        dist='poisson',
                        data=dset)
    if(distribution == 'zero-inflated negative binomial')
        fit <- zeroinfl(as.formula(equation),
                        offset=log(bhpmonth),
                        dist='negbin',
                        data=dset)
    fset <- format_model(fit, dset)
    if(exposure & distribution == 'Poisson'){
        mset <- mset %>% drop_na()
        mset <- mset %>% mutate(predicted=predict(fit, type='response'),
                                group=cut(predicted, c(0:14, Inf)))
        mset <- mset %>% group_by(group)
        sset <- mset %>% summarize(mean=mean(predicted),
                                   variance=var(outvis))
        gg <- ggplot(sset)
        gg <- gg+geom_point(aes(x=mean, y=variance),
                            colour='purple',
                            size=2)
        gg <- gg+geom_abline(slope=1, colour='seagreen')
        gg <- gg+geom_text(data=sset %>% filter(group == '(11,12]'),
                           aes(x=12, y=mean),
                           label='x = y',
                           colour='seagreen',
                           angle=5,
                           vjust=-2,
                           size=6)
        gg <- gg+scale_x_continuous(name='\nMean number of outpatient visits',
                                    limits=c(0, 15),
                                    breaks=seq(0, 14, by=2),
                                    expand=c(0, 0))
        gg <- gg+scale_y_continuous(name='Variance of number of outpatient visits\n',
                                    limits=c(0, 100),
                                    breaks=seq(0, 100, by=20),
                                    expand=c(0, 0))
        print(gg)
    }
    tibble(Model=distribution,
           Exposure=ifelse(exposure, 'Yes', 'No'),
           AIC=AIC(fit),
           BIC=BIC(fit))
}
pne <- bhp_model_table(dset, distribution='Poisson', exposure=FALSE)
pwe <- bhp_model_table(dset, distribution='Poisson', exposure=TRUE)
nbwe <- bhp_model_table(dset, distribution='negative binomial', exposure=TRUE)
zpwe <- bhp_model_table(dset, distribution='zero-inflated Poisson', exposure=TRUE)
znbwe <- bhp_model_table(dset, distribution='zero-inflated negative binomial', exposure=TRUE)
abset <- bind_rows(pne, pwe, nbwe, zpwe, znbwe)

##################################################
# Prostate cancer mortality
##################################################
seer_model_helper <- function(dset, exposure=FALSE, interact=FALSE, fitonly=FALSE){
    response <- 'Deaths'
    predictor <- c('Race', 'Year')
    if(exposure)
        predictor <- c(predictor, 'offset(log(Population))')
    if(interact)
        predictor <- c(predictor, 'Race:Year')
    predictor <- paste(predictor, collapse='+')
    equation <- paste(response, predictor, sep='~')
    fit <- glm(as.formula(equation),
               family=poisson(link='log'),
               data=dset)
    if(fitonly)
        return(fit)
    cset <- coef(summary(fit))
    cset <- data.frame(Predictor=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% transmute(Predictor=sub('(Intercept)', 'Intercept', Predictor, fixed=TRUE),
                               Predictor=gsub(paste0('(', paste(names(dset), collapse='|'), ')'), '\\1 = ', Predictor),
                               Predictor=sub('Year =', 'Year', Predictor),
                               Predictor=sub(':', '$\\times $', Predictor, fixed=TRUE),
                               'Relative risk'=sprintf('%4.2f', exp(Estimate)),
                               '95% CI'=paste0('(',
                                               sprintf('%4.2f', exp(Estimate-1.96*`Std. Error`)),
                                               ', ',
                                               sprintf('%4.2f', exp(Estimate+1.96*`Std. Error`)),
                                               ')'),
                               'P-value'=sub('(.*)', '$\\1 $', `Pr(>|z|)`))
    cset <- cset %>% filter(Predictor != 'Intercept')
    cset <- data.frame(Model=case_when(!exposure ~ 1,
                                       exposure & !interact ~ 2,
                                       exposure & interact ~ 3),
                       Exposure=ifelse(exposure, 'Yes', 'No'),
                       Interaction=ifelse(interact, 'Yes', 'No'),
                       cset,
                       stringsAsFactors=FALSE,
                       check.names=FALSE)
    return(cset)
}

seer_model_table <- function(dset){
    mset <- dset %>% mutate(Year=Year-min(Year))
    model1 <- seer_model_helper(mset, exposure=FALSE, interact=FALSE)
    model2 <- seer_model_helper(mset, exposure=TRUE, interact=FALSE)
    model3 <- seer_model_helper(mset, exposure=TRUE, interact=TRUE)
    models <- bind_rows(model1, model2, model3)
    print(models)
    return(seer_model_helper(dset, exposure=TRUE, interact=TRUE, fitonly=TRUE))
}
fit <- seer_model_table(sset)

##################################################
# Visualize observed and fitted prostate cancer mortality
##################################################
seer_plot <- function(dset, fit){
    dset <- dset %>% mutate(Race=relevel(Race, ref='Black'),
                            Observed=1e5*Deaths/Population,
                            Fitted=1e5*fitted(fit)/Population)
    theme_update(legend.position='none')
    gg <- ggplot(dset)
    gg <- gg+geom_point(aes(x=Year, y=Observed, colour=Race),
                        size=2)
    gg <- gg+geom_line(aes(x=Year, y=Fitted, colour=Race),
                       size=1)
    gg <- gg+geom_text(data=dset %>% filter(Year == 2005),
                       aes(x=Year, y=Fitted, colour=Race, label=Race),
                       vjust=-2,
                       size=6)
    gg <- gg+scale_x_continuous(name='\nYear of death')
    gg <- gg+scale_y_continuous(name='Death rate per 100,000 men\n',
                                limits=c(0, 160),
                                breaks=seq(0, 150, by=50),
                                expand=c(0, 0))
    gg <- gg+scale_colour_manual(name='', values=c(White=muted('#56B4E9'), Black=muted('seagreen')))
    print(gg)
}
seer_plot(sset, fit)

