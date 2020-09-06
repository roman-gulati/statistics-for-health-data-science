##################################################
# Create figures and tables for Chapter 7
##################################################
library(tidyverse)
library(xtable)
library(grid)
library(scales)
library(viridis)
library(rsample)

if(Sys.getenv('RSTUDIO') == '1')
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('shared.R')

set.seed(1111)

##################################################
# Health expenditures data
##################################################
grab_and_curate_meps <- function(year){
  dset <- grab_meps(year)
  dset <- dset %>% select(Arthritis=ARTHDX,
                          Age=AGE17X,
                          Sex=SEX,
                          Race=RACETHX,
                          Cost=TOTEXP17,
                          Diabetes=DIABDX)
  dset <- dset %>% mutate(Arthritis=factor(Arthritis,
                                           levels=c(-9, -8, -7, -1, 2, 1),
                                           labels=c('not ascertained',
                                                    'don\'t know',
                                                    'refused',
                                                    'inapplicable',
                                                    'no',
                                                    'yes')),
                          Sex=factor(Sex,
                                     levels=c(1, 2),
                                     labels=c('male', 'female')),
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
                                          levels=c(-9, -8, -7, -1, 2, 1),
                                          labels=c('not ascertained',
                                                   'don\'t know',
                                                   'refused',
                                                   'inapplicable',
                                                   'no',
                                                   'yes')),
                          Anycost=Cost > 0,
                          Logcost=log10(Cost+1))
  cat('Remove',
      dset %>% filter(Age <= 17) %>% nrow(),
      'age 17 or younger\n')
  dset <- dset %>% filter(Age > 17)
  dset %>% droplevels()
}
#pset <- grab_and_curate_meps(2017)

##################################################
# Select a sample of 2000 records
##################################################
#sset <- pset %>% sample_n(2000)

##################################################
# Visualize histograms
##################################################
two_panel_plot <- function(dset1,
                           dset2,
                           arg1,
                           arg2,
                           fun,
                           ...,
                           filename=NULL,
                           ext='pdf',
                           saveit=FALSE){
    panel1 <- fun(dset1, arg1, ...)
    panel2 <- fun(dset2, arg2, ...)
    Layout <- grid.layout(ncol=2,
                          nrow=1,
                          heights=unit(5, 'null'))
    if(saveit){
        filename <- paste(filename, ext, sep='.')
        get(ext)(here('figures', filename),
                 height=5,
                 width=10)
        on.exit(graphics.off())
    }
    grid.newpage()
    pushViewport(viewport(layout=Layout))
    print(panel1, vp=viewport(layout.pos.col=1, layout.pos.row=1))
    print(panel2, vp=viewport(layout.pos.col=2, layout.pos.row=1))
}

cost_histogram <- function(dset, panel_title, nbreaks=30, xmax=6){
    hdat <- hist(dset[['Logcost']], breaks=nbreaks, plot=FALSE)
    hset <- with(hdat, data.frame(mids, density))
    cat('Median of costs in', tolower(panel_title), ':       ', median(dset$Cost), '\n')
    cat('Median of log(costs+1) in', tolower(panel_title), ':', median(dset$Logcost), '\n')
    gg_theme()
    gg <- ggplot(hset)
    gg <- gg+ggtitle(paste0(panel_title, ' (n=', nrow(dset), ')'))
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+scale_x_continuous(name=parse(text='Total~medical~expenditures~(log[10])'),
                                limits=c(0, xmax+0.02),
                                breaks=seq(0, xmax, by=1),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Density\n',
                                limits=c(0, 1),
                                expand=c(0, 0))
    gg
}
#two_panel_plot(pset,
#               sset,
#               arg1='Population',
#               arg2='Sample',
#               fun=cost_histogram,
#               filename='07-cost_histogram',
#               saveit=TRUE)

##################################################
# Bootstrapping the median
##################################################
median_cost_histogram <- function(dset, panel_title, B=5000, sample_size=2000, nbreaks=30){
    original_size <- nrow(dset)
    if(tolower(panel_title) == 'population')
        dset <- dset %>% sample_n(size=sample_size, replace=TRUE)
    rset <- bootstraps(dset, times=B)
    ests <- map_dbl(rset$splits, function(x) median(as.data.frame(x)$Cost))
    hdat <- hist(ests, breaks=nbreaks, plot=FALSE)
    hset <- with(hdat, data.frame(mids, density))
    cat('Range of median costs in', tolower(panel_title), ':', min(ests), max(ests), '\n')
    cat('SD of median costs in', tolower(panel_title), ':', sd(ests), '\n')
    gg_theme()
    gg <- ggplot(hset)
    gg <- gg+ggtitle(paste0(panel_title, ' (n=', original_size, ')'))
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+scale_x_continuous(name='Median medical expenditures ($)',
                                limits=c(800, 1520),
                                breaks=seq(800, 1500, by=100))
    gg <- gg+scale_y_continuous(name='Density\n',
                                limits=c(0, 0.008),
                                expand=c(0, 0))
    gg
}
#two_panel_plot(pset,
#               sset,
#               arg1='Population',
#               arg2='Sample',
#               fun=median_cost_histogram,
#               filename='07-median_cost_histogram',
#               saveit=TRUE)

##################################################
# Illustrate bootstrap methods to study the probability
# of no medical costs
##################################################
format_model <- function(fit, dset){
    cset <- coef(summary(fit))
    cset <- data.frame(Predictor=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% transmute(Predictor=sub('(Intercept)', 'Intercept', Predictor, fixed=TRUE),
                               Predictor=gsub(paste0('(', paste(names(dset), collapse='|'), ')'), '\\1 = ', Predictor),
                               Predictor=sub('Age =', 'Age', Predictor),
                               Coefficient=Estimate,
                               'Maximum likelihood SE'=`Std. Error`)
    cset
}

simple_histogram <- function(dset, panel_title, xlabel, nbreaks=30){
    hdat <- hist(dset, breaks=nbreaks, plot=FALSE)
    hset <- with(hdat, data.frame(mids, density))
    gg_theme()
    gg <- ggplot(hset)
    gg <- gg+ggtitle(panel_title)
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+scale_x_continuous(name=xlabel,
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Density\n',
                                #limits=c(0, 0.001),
                                #breaks=seq(0, 0.001, by=0.0002),
                                #labels=label_number(0.0001),
                                expand=c(0, 0))
    gg
}

bootstrap_no_cost_analysis <- function(dset, B=5000, alpha=0.05, saveit=FALSE){
    # logistic regression for no medical costs
    dset <- dset %>% mutate(Nocost=!Anycost)
    fit <- glm(Nocost~Age+Sex,
               data=dset,
               family=binomial(link='logit'))
    fset <- format_model(fit, dset)
    # bootstrap SEs for regression coefficients
    rset <- bootstraps(dset, times=B)
    ests <- map(rset$splits, function(x){
                    rfit <- glm(Nocost~Age+Sex, 
                                data=as.data.frame(x),
                                family=binomial(link='logit'))
                    coef(rfit)
               })
    bset <- do.call(rbind, ests)
    bses <- apply(bset, 2, sd)
    fset <- fset %>% mutate('Bootstrap SE'=bses)
    if(saveit){
        filename <- '07-bootstrap_no_cost.tex'
        Caption <- 'Fitted logistic regression of no medical expenditures based on a sample of 2000 persons age 18 years or older in the MEPS 2017 data with maximum likelihood and bootstrap standard errors (SEs) based on 5000 bootstrap samples'
        print(xtable(fset,
                     digits=4,
                     align='llccc',
                     label='tab:bootstrap_no_cost',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              math.style.negative=TRUE,
              hline.after=0)
    }
    # visualize regression coefficients
    bset <- as_tibble(bset)
    two_panel_plot(bset$Age,
                   bset$Sexfemale,
                   arg1='Coefficient for age',
                   arg2='Coefficient for sex',
                   fun=simple_histogram,
                   xlabel='Coefficient values',
                   filename='07-no_cost_coefficient_histogram',
                   saveit=saveit)
    # marginal additive effect for sex in the sample
    wset <- dset %>% mutate(Sex='female')
    mset <- dset %>% mutate(Sex='male')
    dset <- dset %>% mutate(Female=predict(fit, newdata=wset, type='response'),
                            Male=predict(fit, newdata=mset, type='response'),
                            MAE=Female-Male)
    cat('Marginal additive effect for sex:', mean(dset$MAE), '\n')
    # bootstrap standard error for the marginal additive effect for sex
    pred <- map(rset$splits, function(x){
                    rfit <- glm(Nocost~Age+Sex, 
                                data=as.data.frame(x),
                                family=binomial(link='logit'))
                    pset <- dset %>% mutate(Female=predict(rfit, newdata=wset, type='response'),
                                            Male=predict(rfit, newdata=mset, type='response'),
                                            MAE=Female-Male)
                    mean(pset$MAE)
                   })
    maes <- do.call(c, pred)
    cat('SE for marginal additive effect for sex:', sd(maes), '\n')
    gg <- simple_histogram(maes, '', 'Marginal additive effect for sex')
    if(saveit)
        ggsave(plot=gg,
               file=here('figures', '07-marginal_additive_effect_for_sex.pdf'),
               height=5,
               width=10)
    # predict probability of no costs for hypothetical group
    pred <- map(rset$splits, function(x){
                    rfit <- glm(Nocost~Age+Sex, 
                                data=as.data.frame(x),
                                family=binomial(link='logit'))
                    predict(rfit, type='response')
                   })
    pset <- do.call(rbind, pred)
    pest <- rowMeans(pset)
    cat('Bootstrap SE for predicted proportion with no costs:', sd(pest), '\n')
    # calculate bootstrap-t and bootstrap percentiles confidence intervals
    mult <- qnorm(1-alpha/2, lower.tail=TRUE)
    fset <- fset %>% rename(Parameter='Predictor',
                            Estimate='Coefficient')
    fset <- fset %>% mutate('Bootstrap-t'=paste0('(',
                                                 sprintf('%5.3f', Estimate-mult*`Bootstrap SE`),
                                                 ', ',
                                                 sprintf('%5.3f', Estimate+mult*`Bootstrap SE`),
                                                 ')'),
                            'Bootstrap percentile'=paste0('(',
                                                 sprintf('%5.3f', apply(bset, 2, quantile, alpha/2)),
                                                 ', ',
                                                 sprintf('%5.3f', apply(bset, 2, quantile, 1-alpha/2)),
                                                 ')'))
    fset <- fset %>% mutate(Estimate=sprintf('%5.3f', Estimate),
                            'Bootstrap SE'=sprintf('%5.3f', `Bootstrap SE`))
    fset <- fset %>% select(-contains('Maximum likelihood'))
    ests <- map_dbl(rset$splits, function(x) median(as.data.frame(x)$Cost))
    mcis <- tibble(Parameter='Median cost',
                   Estimate=median(dset$Cost),
                   'Bootstrap SE'=sd(ests),
                   'Bootstrap-t'=paste0('(',
                                        sprintf('%4.0f', Estimate-mult*sd(ests)),
                                        ', ',
                                        sprintf('%4.0f', Estimate+mult*sd(ests)),
                                        ')'),
                   'Bootstrap percentile'=paste0('(',
                                                 sprintf('%4.0f', quantile(ests, alpha/2)),
                                                 ', ',
                                                 sprintf('%4.0f', quantile(ests, 1-alpha/2)),
                                                 ')'))
    mcis <- mcis %>% mutate(Estimate=sprintf('%4.0f', Estimate),
                            'Bootstrap SE'=sprintf('%2.0f', `Bootstrap SE`))
    pcis <- tibble(Parameter='Probability of no costs',
                   Estimate=mean(pest),
                   'Bootstrap SE'=sd(pest),
                   'Bootstrap-t'=paste0('(',
                                        sprintf('%5.3f', Estimate-mult*sd(pest)),
                                        ', ',
                                        sprintf('%5.3f', Estimate+mult*sd(pest)),
                                        ')'),
                   'Bootstrap percentile'=paste0('(',
                                                 sprintf('%5.3f', quantile(pest, alpha/2)),
                                                 ', ',
                                                 sprintf('%5.3f', quantile(pest, 1-alpha/2)),
                                                 ')'))
    pcis <- pcis %>% mutate(Estimate=sprintf('%5.3f', Estimate),
                            'Bootstrap SE'=sprintf('%5.3f', `Bootstrap SE`))
    cset <- bind_rows(fset, mcis, pcis)
    cset <- cset %>% mutate(Estimate=sub('-', '$-$', Estimate),
                            'Bootstrap-t'=gsub('-', '$-$', `Bootstrap-t`),
                            'Bootstrap percentile'=gsub('-', '$-$', `Bootstrap percentile`))
    if(saveit){
        filename <- '07-bootstrap_confidence_intervals.tex'
        Caption <- 'Bootstrap-t and bootstrap percentile 95\\% confidence intervals for logistic regression coefficients and median medical expenditures based on a sample of 2000 persons age 18 years or older in the MEPS 2017 data'
        print(xtable(cset,
                     digits=0,
                     align='llcccc',
                     label='tab:bootstrap_confidence_intervals',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.text.function=identity,
              math.style.negative=TRUE,
              hline.after=c(0, nrow(fset)))
    }
}
#bootstrap_no_cost_analysis(sset, saveit=TRUE)

##################################################
# Bootstrap SE for marginal adjusted additive effect
# for sex from a two-part model
##################################################
bootstrap_two_part_analysis <- function(dset, B=1000, saveit=FALSE){
    dset <- dset %>% filter(Diabetes %in% c('no', 'yes'))
    dset <- dset %>% mutate(ScaledAge=Age-18)
    bfit <- glm(Anycost~ScaledAge+Sex+Race+Diabetes,
                family=binomial(link='logit'),
                data=dset)
    gfit <- glm(Cost~ScaledAge+Sex+Race+Diabetes,
                family=Gamma(link='log'),
                data=dset,
                subset=Anycost)
    bpred0 <- predict(bfit, newdata=dset %>% mutate(Diabetes='no'), type='response')
    bpred1 <- predict(bfit, newdata=dset %>% mutate(Diabetes='yes'), type='response')
    gpred0 <- predict(gfit, newdata=dset %>% mutate(Diabetes='no'), type='response')
    gpred1 <- predict(gfit, newdata=dset %>% mutate(Diabetes='yes'), type='response')
    mean0 <- mean(bpred0*gpred0)
    mean1 <- mean(bpred1*gpred1)
    mae <- mean1-mean0
    cat('Marginal additive effect for diabetes:', round(mae), '\n')
    rset <- bootstraps(dset, times=B)
    pred <- map(rset$splits, function(x){
                    rbfit <- glm(Anycost~ScaledAge+Sex+Race+Diabetes,
                                 family=binomial(link='logit'),
                                 data=as.data.frame(x))
                    rgfit <- glm(Cost~ScaledAge+Sex+Race+Diabetes,
                                 family=Gamma(link='log'),
                                 data=as.data.frame(x),
                                 subset=Anycost)
                    nset <- dset %>% mutate(Diabetes='no')
                    yset <- dset %>% mutate(Diabetes='yes')
                    bpred0 <- predict(rbfit, newdata=nset, type='response')
                    bpred1 <- predict(rbfit, newdata=yset, type='response')
                    gpred0 <- predict(rgfit, newdata=nset, type='response')
                    gpred1 <- predict(rgfit, newdata=yset, type='response')
                    mean(bpred1*gpred1)-mean(bpred0*gpred0)
                   })
    maes <- do.call(c, pred)
    cat('SE for marginal adjusted additive effect for diabetes:', sd(maes), '\n')
    gg <- simple_histogram(maes, '', 'Marginal additive effect for diabetes')
    if(saveit)
        ggsave(plot=gg,
               file=here('figures', '07-marginal_additive_effect_for_diabetes.pdf'),
               height=5,
               width=10)
}
#bootstrap_two_part_analysis(pset, saveit=TRUE)

##################################################
# Test differences in distributions of medical
# expenditures between persons with diabetes vs
# persons with arthritis
##################################################
obs_exp_table <- function(dset, group){
    chi2 <- chisq.test(dset[[group]], dset[['Strata']])
    dset %>% with(data.frame(Group=group,
                             Strata=levels(Strata),
                             Observed=chi2$observed['yes', ],
                             Expected=chi2$expected['yes', ],
                             stringsAsFactors=FALSE,
                             row.names=NULL))
}

diabetes_vs_arthritis_histogram <- function(dset, panel_title, value, nbreaks=30){
    hdat <- hist(dset, breaks=nbreaks, plot=FALSE)
    hset <- with(hdat, data.frame(mids, density))
    gg_theme()
    gg <- ggplot(hset)
    gg <- gg+ggtitle(panel_title)
    gg <- gg+geom_bar(aes(x=mids, y=density), stat='identity')
    gg <- gg+geom_vline(xintercept=value)
    gg <- gg+scale_x_continuous(name='Q statistic',
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Density\n',
                                limits=c(0, 0.2),
                                expand=c(0, 0))
    gg
}

diabetes_vs_arthritis_analysis <- function(dset, B=5000, alpha=0.05, saveit=FALSE){
    tset <- dset %>% filter(xor(Arthritis == 'yes', Diabetes == 'yes'))
    tset <- tset %>% droplevels()
    cat(tset %>% filter(Diabetes == 'yes') %>% nrow(), 'persons have diabetes\n')
    cat(tset %>% filter(Arthritis == 'yes') %>% nrow(), 'persons have arthritis\n')
    tset <- tset %>% mutate(Strata=cut(Cost,
                                       breaks=c(0, 500, 1000, 5000, 10000, Inf),
                                       label=c('$<\\$500$',
                                               '\\$500--\\$9999',
                                               '\\$1000--\\$4999',
                                               '\\$5000--\\$99999',
                                               '$\\ge\\$10000$'),
                                       right=FALSE))
    diab <- obs_exp_table(tset, 'Diabetes')
    arth <- obs_exp_table(tset, 'Arthritis')
    gset <- bind_rows(diab, arth)
    gset <- gset %>% reshape::melt(id.vars=c('Group', 'Strata'))
    gset <- gset %>% rename(Statistic='variable')
    # calculate sum of squared differences
    cset <- gset %>% group_by(Strata, Group)
    cset <- cset %>% summarize(Exp=value[Statistic == 'Expected'],
                               Obs=value[Statistic == 'Observed'],
                               Term=(Exp-Obs)^2/Exp)
    cset <- cset %>% ungroup()
    cest <- cset %>% with(sum(Term))
    cat('Sum of squared differences:', cest, '\n')
    # format table of observed and expected counts
    gset <- gset %>% mutate(Strata=factor(Strata, levels=levels(tset$Strata)),
                            value=ifelse(Statistic == 'Observed',
                                         sprintf('%2.0f', value),
                                         sprintf('%5.2f', value)))
    gset <- gset %>% reshape::cast(Statistic+Group~Strata, value='value')
    gset <- gset %>% mutate(Statistic=factor(Statistic, levels=c('Observed', 'Expected')),
                            Group=factor(Group, levels=c('Diabetes', 'Arthritis')))
    gset <- gset %>% arrange(Statistic, Group)
    if(saveit){
        filename <- '07-diabetes_vs_arthritis.tex'
        Caption <- 'Observed and expected counts of persons with prior diagnosis of diabetes or arthritis based on a sample of 2000 persons from the MEPS 2017 data'
        print(xtable(gset,
                     digits=0,
                     align='lllccccc',
                     label='tab:diabetes_vs_arthritis',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.text.function=identity,
              math.style.negative=TRUE,
              hline.after=c(0, 2))
    }
    # bootstrap p-value
    group <- tset %>% with(Diabetes)
    strata <- tset %>% with(Strata)
    scounts <- table(strata)
    nstrata <- length(scounts)
    npersons <- length(group)
    bootstrap_statistics <- vector(mode='numeric', length=B)
    permutation_statistics <- vector(mode='numeric', length=B)
    for (b in 1:B){
        rstrata <- sample(x=nstrata, size=npersons, replace=TRUE, prob=scounts)
        bootstrap_statistics[b] <- chisq.test(rstrata, group)$statistic
        rstrata <- sample(x=strata)
        permutation_statistics[b] <- chisq.test(rstrata, group)$statistic
    }
    bootstrap_pvalue <- mean(bootstrap_statistics >= cest)
    permutation_pvalue <- mean(permutation_statistics >= cest)
    cat('Bootstrap p-value:', bootstrap_pvalue, '\n')
    cat('Permutation p-value:', permutation_pvalue, '\n')
    two_panel_plot(bootstrap_statistics,
                   permutation_statistics,
                   arg1='Bootstrap samples',
                   arg2='Permutations',
                   fun=diabetes_vs_arthritis_histogram,
                   value=cest,
                   filename='07-diabetes_vs_arthritis_histogram',
                   saveit=saveit)
}
#diabetes_vs_arthritis_analysis(sset, saveit=TRUE)

