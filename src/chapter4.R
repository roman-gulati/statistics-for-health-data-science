##################################################
# Chapter 4 examples
##################################################
library(tidyverse)
library(scales)
library(grid)
library(viridis)
library(margins)
library(pROC)
library(nnet)

source('grab_nhanes.R')

##################################################
# Grab NHANES data for initial explorations
##################################################
iset <- bind_rows(grab_and_curate_nhanes(1999, minage=20, maxage=59),
                  grab_and_curate_nhanes(2015, minage=20, maxage=59))
iset <- iset %>% mutate(BMICAT=cut(BMI,
                                   breaks=c(0, 18.5, 25, 30, 100),
                                   labels=c('Underweight',
                                            'Normal',
                                            'Overweight',
                                            'Obese'),
                                   right=FALSE),
                        OBESE=factor(BMI >= 30,
                                     levels=c('FALSE', 'TRUE'),
                                     labels=c('Not obese', 'Obese')))

##############################################
# Two-way table for obese persons by year
##############################################
obese_frequency_table <- function(dset){
    dset <- dset %>% mutate(YEAR=sub('-', '--', YEAR))
    fset <- dset %>% group_by(YEAR, OBESE)
    fset <- fset %>% summarize(Count=n())
    fset <- fset %>% group_by(YEAR)
    fset <- fset %>% mutate(Percent=sprintf('%4.1f%%', 100*Count/sum(Count)),
                            Content=paste0(Count, ' (', Percent, ')'),
                            Total=paste(sum(Count), '(100.0%)'))
    fset <- fset %>% select(-Count, -Percent)
    fset <- fset %>% spread(OBESE, Content)
    fset <- fset %>% select(Year=YEAR, 'Not obese', Obese, Total)
    return(fset)
}
obese_frequency_table(iset)

##################################################
# Compare OR versus RR
##################################################
or_versus_rr_table <- function(){
    dset <- tibble(X0=rep(c(0.01, 0.05, 0.10), each=2),
                   X1=c(0.02, 0.04, 0.10, 0.20, 0.20, 0.40))
    dset <- dset %>% mutate(RR=X1/X0,
                            O0=X0/(1-X0),
                            O1=X1/(1-X1),
                            OR=O1/O0)
    dset <- dset %>% select('$P(Y=1 \\mid X=0)$'=X0,
                            '$P(Y=1 \\mid X=1)$'=X1,
                            RR,
                            OR)
    return(dset)
}

##################################################
# Proportion of obese persons by age 
##################################################
count_scatterplot <- function(dset){
    theme_update(legend.position='none')
    gg <- ggplot(dset, aes(x=AGE, y=OBESE))
    gg <- gg+geom_count(color='gray')
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_discrete(name='')
    print(gg)
}
count_scatterplot(iset)

proportion_scatterplot <- function(dset){
    dset <- dset %>% mutate(OBESE=as.integer(OBESE)-1)
    gset <- dset %>% group_by(AGE)
    gset <- gset %>% summarize(Proportion=sum(OBESE)/n())
    gg <- ggplot(gset, aes(x=AGE, y=Proportion))
    gg <- gg+geom_point(color='gray')
    gg <- gg+geom_smooth(data=dset,
                         aes(x=AGE, y=OBESE),
                         method='glm',
                         size=1,
                         color=muted('purple'),
                         method.args=list(family='binomial'),
                         se=FALSE)
    gg <- gg+geom_smooth(method='lm',
                         size=1,
                         color=muted('seagreen'),
                         se=FALSE)
    gg <- gg+geom_text(data=gset %>% filter(AGE == max(AGE)),
                       aes(x=AGE,
                           y=0.4,
                           label='Logistic regression'),
                       colour=muted('purple'),
                       angle=5,
                       hjust=1,
                       vjust=-4,
                       size=5)
    gg <- gg+geom_text(data=gset %>% filter(AGE == max(AGE)),
                       aes(x=AGE,
                           y=0.4,
                           label='Linear regression'),
                       colour=muted('seagreen'),
                       angle=5,
                       hjust=1,
                       vjust=2,
                       size=5)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Percent of obese persons',
                                limits=c(0, 1),
                                breaks=seq(0, 1, by=0.2),
                                labels=percent_format(accuracy=1),
                                expand=c(0, 0))
    print(gg)
}
proportion_scatterplot(iset)

##################################################
# Logistic regression
##################################################
format_model <- function(fit, dset){
    cset <- coef(summary(fit))
    cset <- data.frame(Predictor=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% transmute(Predictor=sub('(Intercept)',
                                             'Intercept',
                                             Predictor,
                                             fixed=TRUE),
                               Predictor=gsub(paste0('(',
                                                     paste(names(dset), collapse='|'),
                                                     ')'),
                                              '\\1 = ',
                                              Predictor),
                               Predictor=sub(' = $', '', Predictor),
                               Predictor=sub('2015-2016', '2015--2016', Predictor),
                               Coefficient=Estimate,
                               SE=`Std. Error`,
                               'P-value'=`Pr(>|z|)`)
    return(cset)
}

logistic_regression_table <- function(dset){
    fit <- glm(OBESE~YEAR,
               family=binomial(link='logit'),
               data=dset)
    fset <- format_model(fit, dset)
    return(fset)
}
logistic_regression_table(iset)

##################################################
# Univariate summary table
##################################################
univariate_summary <- function(dset, varname, reflevel){
    varsym <- sym(varname)
    uset <- dset %>% group_by(!!varsym)
    uset <- uset %>% summarize(Total=n(),
                               N.obese=sum(OBESE == 'Obese'),
                               P.obese=N.obese/Total,
                               Percent.obese=sprintf('%4.1f%%', 100*P.obese),
                               Content.obese=paste0(N.obese, ' (', Percent.obese, ')'),
                               Odds.obese=P.obese/(1-P.obese))
    uset <- uset %>% ungroup()
    uset <- uset %>% mutate(Variable=paste(varname, !!varsym, sep=' = '),
                            OR=Odds.obese/Odds.obese[!!varsym == reflevel],
                            OR=ifelse(!!varsym == reflevel, 'reference', sprintf('%4.2f', OR)))
    uset <- uset %>% select(Variable,
                            'Sample, n'=Total,
                            'Obese, n (%)'=Content.obese,
                            'Odds ratio'=OR)
    return(uset)
}

univariate_table <- function(dset){
    dset <- dset %>% mutate(YEAR=sub('-', '--', YEAR))
    sset <- univariate_summary(dset, 'SEX', 'Male')
    rset <- univariate_summary(dset, 'RACE', 'Non-Hispanic White')
    yset <- univariate_summary(dset, 'YEAR', '1999--2000')
    uset <- bind_rows(sset, rset, yset)
    return(uset)
}
univariate_table(iset)

##################################################
# Multivariate logistic regression
##################################################
multivariate_table <- function(dset, alpha=0.05){
    dset <- dset %>% mutate(AGE=AGE-min(AGE))
    fit <- glm(OBESE~SEX+RACE+AGE+YEAR,
               family=binomial(link='logit'),
               data=dset)
    fset <- format_model(fit, dset)
    mult <- qnorm(alpha/2, lower.tail=FALSE)
    fset <- fset %>% mutate(OR=ifelse(Predictor == 'Intercept',
                                      NA,
                                      sprintf('%4.2f', exp(Coefficient))),
                            '95% CI'=ifelse(Predictor == 'Intercept',
                                      NA,
                                      paste0('(',
                                             sprintf('%4.2f', exp(Coefficient-mult*SE)),
                                             ',',
                                             sprintf('%4.2f', exp(Coefficient+mult*SE)),
                                             ')')))
    return(fit)
}
fit <- multivariate_table(iset)

##################################################
# Predicted probabilities of obesity by sex and year
# restricted to RACE = Other/Mixed
##################################################
prediction_plot <- function(fit, dset){
    minage <- min(dset$AGE)
    dset <- dset %>% mutate(AGE=AGE-minage)
    nset <- dset %>% with(expand.grid(SEX=c('Male', 'Female'),
                                      RACE='Other/Mixed',
                                      AGE=seq(min(AGE), max(AGE)),
                                      YEAR=c('1999-2000', '2015-2016')))
    nset <- nset %>% mutate(PHAT=predict(fit, newdata=nset, 'response'),
                            AGE=AGE+minage)
    nset <- nset %>% group_by(SEX, AGE)
    gg <- ggplot(nset)
    gg <- gg+geom_line(aes(x=AGE,
                           y=PHAT,
                           colour=SEX,
                           linetype=YEAR),
                       size=1)
    gg <- gg+geom_text(data=nset %>% filter(AGE == 40),
                       aes(x=AGE,
                           y=PHAT,
                           colour=SEX,
                           label=paste(SEX, YEAR, sep=', '),
                           vjust=ifelse(YEAR == '1999-2000', 2, -1)),
                       angle=5,
                       size=5)
    gg <- gg+scale_colour_viridis(discrete=TRUE, begin=0.2, end=0.5)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Predicted probability of obesity',
                                limits=c(0, 0.5),
                                breaks=seq(0, 0.5, by=0.1),
                                labels=percent_format(accuracy=1),
                                expand=c(0, 0))
    print(gg)
}
prediction_plot(fit, iset)

##################################################
# Compare differences, RRs, and ORs for obesity by
# sex and year restricted to RACE = Other/Mixed
##################################################
prediction_comparison <- function(fit, dset, comparison=FALSE){
    minage <- min(dset$AGE)
    dset <- dset %>% mutate(AGE=AGE-minage)
    nset <- dset %>% with(expand.grid(SEX=c('Male', 'Female'),
                                      RACE='Other/Mixed',
                                      AGE=seq(min(AGE), max(AGE)),
                                      YEAR=c('1999-2000', '2015-2016')))
    nset <- nset %>% mutate(PHAT=predict(fit, newdata=nset, 'response'),
                            AGE=AGE+minage)
    nset <- nset %>% group_by(SEX, AGE)
    nset <- nset %>% mutate(Difference=diff(PHAT),
                            'Risk ratio'=PHAT[YEAR == '2015-2016']/PHAT[YEAR == '1999-2000'],
                            'Odds ratio'=PHAT[YEAR == '2015-2016']/(1-PHAT[YEAR == '2015-2016'])/
                                        (PHAT[YEAR == '1999-2000']/(1-PHAT[YEAR == '1999-2000'])))
    nset <- nset %>% reshape2::melt(measure.vars=c('Difference', 'Risk ratio', 'Odds ratio'))
    nset <- nset %>% rename(Measure='variable', Value='value')
    gg <- ggplot(nset)
    gg <- gg+geom_line(aes(x=AGE,
                           y=Value,
                           colour=SEX,
                           linetype=YEAR),
                       size=1)
    gg <- gg+facet_grid(Measure~., scales='free_y')
    gg <- gg+scale_colour_viridis(discrete=TRUE, begin=0.2, end=0.5)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Predicted effect on obesity',
                                expand=c(0, 0))
    print(gg)
}
prediction_comparison(fit, iset)

##################################################
# Marginal additive effect of year
##################################################
marginal_additive_effect <- function(fit, dset, variable){
    mset <- margins(fit, data=dset, variables=variable)
    print(summary(mset))
}
marginal_additive_effect(fit, iset, 'YEAR')

lm_check <- function(dset, vars){
    dset <- dset %>% mutate(AGE=AGE-min(AGE),
                            OBESE=as.integer(OBESE)-1)
    lmod <- lm(OBESE~SEX+RACE+AGE+YEAR, data=dset)
    ests <- cbind(coef(lmod), confint(lmod))
    print(ests[grepl(vars, rownames(ests)), ])
}
lm_check(iset, 'YEAR')

##################################################
# AIC and BIC for all 8 possible models with YEAR
##################################################
aicbic_row <- function(vset, dset, response='OBESE'){
    predictor <- with(vset, paste(AGE, RACE, SEX, YEAR, sep='+'))
    predictor <- gsub('\\+\\+*', '\\+', predictor)
    equation <- paste(response, predictor, sep='~')
    fit <- glm(as.formula(equation),
               family=binomial(link='logit'),
               data=dset)
    vset <- vset %>% mutate(AIC=AIC(fit), BIC=BIC(fit))
    return(vset)
}

aicbic_table <- function(dset){
    dset <- dset %>% mutate(AGE=AGE-min(AGE))
    vset <- expand.grid(AGE=c('', 'AGE'),
                        RACE=c('', 'RACE'),
                        SEX=c('', 'SEX'),
                        stringsAsFactors=FALSE)
    vset <- vset %>% mutate(YEAR='YEAR',
                            Model=seq(nrow(vset)))
    vset <- vset %>% group_by(Model)
    vset <- vset %>% do(aicbic_row(., dset))
    vset <- vset %>% select(Model, SEX, RACE, AGE, YEAR, AIC, BIC)
    vset <- vset %>% mutate(SEX=sub('SEX', '+', SEX),
                            RACE=sub('RACE', '+', RACE),
                            AGE=sub('AGE', '+', AGE),
                            YEAR=sub('YEAR', '+', YEAR))
    vset[vset == ''] <- '--'
    return(vset)
}
aicbic_table(iset)

##################################################
# Hosmer-Lemeshow comparison of observed and predicted
##################################################
hosmer_lemeshow_table <- function(fit, dset, ngroups=10){
    dset <- dset %>% mutate(AGE=AGE-min(AGE))
    pset <- dset %>% mutate(Predicted=predict(fit, type='response'))
    pset <- pset %>% arrange(Predicted)
    pset <- pset %>% mutate(Group=cut(Predicted,
                                      breaks=quantile(Predicted, probs=seq(0, 1, by=1/ngroups)),
                                      include.lowest=TRUE,
                                      right=FALSE),
                            Group=as.integer(Group))
    gset <- pset %>% group_by(Group)
    gset <- gset %>% summarize('Obs.obese'=sum(OBESE == 'Obese'),
                               'Obs.nonobese'=sum(OBESE == 'Not obese'),
                               'Pred.obese'=sum(Predicted),
                               'Pred.nonobese'=sum(1-Predicted),
                               'Res.obese'=(Obs.obese-Pred.obese)/sqrt(Pred.obese),
                               'Res.nonobese'=(Obs.nonobese-Pred.nonobese)/sqrt(Pred.nonobese))
    gset <- gset %>% ungroup()
    cat('\tHosmer-Lemeshow goodness-of-fit test\n',
        'Statistic:', statistic <- gset %>% with(sum(Res.obese^2)+sum(Res.nonobese^2)),
        '\n',
        'P-value:  ', p=1-pchisq(statistic, ngroups-2),
        '\n')
    names(gset)[-1] <- sapply(strsplit(names(gset[-1]), '\\.'),
                              function(x) paste0('$\\text{', x[1], '}_\\text{', x[2], '}$'))
    names(gset) <- sub('nonobese', 'non-obese', names(gset))
}
hosmer_lemeshow_table(fit, iset)

##################################################
# Prediction accuracy 
##################################################
prediction_accuracy_table <- function(fit, dset, threshold=0.5){
    dset <- dset %>% mutate(AGE=AGE-min(AGE))
    pset <- dset %>% mutate(Predicted=predict(fit, type='response'),
                            PREDICTED=Predicted > threshold)
    fset <- pset %>% group_by(OBESE, PREDICTED)
    fset <- fset %>% summarize(Count=n())
    fset <- fset %>% group_by(OBESE)
    fset <- fset %>% mutate(Percent=sprintf('%4.1f%%', 100*Count/sum(Count)),
                            Content=paste0(Count, ' (', Percent, ')'),
                            Total=paste(sum(Count), '(100.0%)'))
    fset <- fset %>% ungroup()
    fset <- fset %>% mutate(OBESE=factor(OBESE,
                                         levels=c('Not obese', 'Obese'),
                                         labels=c('$\\text{Obs}_\\text{non-obese}$',
                                                  '$\\text{Obs}_\\text{obese}$')),
                            Content=sub('%', '\\%', Content, fixed=TRUE),
                            Total=sub('%', '\\%', Total, fixed=TRUE))
    fset <- fset %>% select(-Count, -Percent)
    fset <- fset %>% spread(PREDICTED, Content)
    fset <- fset %>% select(Status=OBESE,
                            '$\\text{Pred}_\\text{non-obese}$'='FALSE',
                            '$\\text{Pred}_\\text{obese}$'='TRUE',
                            Total)
    return(fset)
}
prediction_accuracy_table(fit, iset)

##################################################
# ROC curve
##################################################
roc_curve_plot <- function(fit, dset){
    dset <- dset %>% mutate(AGE=AGE-min(AGE))
    pset <- dset %>% mutate(Predicted=predict(fit, type='response'))
    afit <- roc(OBESE~Predicted, pset, direction='>')
    aset <- with(afit, tibble(thresholds,
                              sensitivities,
                              specificities))
    aset <- aset %>% arrange(desc(specificities), sensitivities)
    aset <- aset %>% mutate(floored_thresholds=signif(thresholds, 0.1))
    tset <- aset %>% group_by(floored_thresholds)
    tset <- tset %>% filter(is.finite(thresholds) &
                            thresholds == min(thresholds),
                            !near(floored_thresholds, 0.3))
    tset <- tset %>% mutate(thresholds=sprintf('%3.1f', thresholds),
                            thresholds=sub('^', 'T=', thresholds))
    gg <- ggplot(data=aset)
    gg <- gg+geom_abline(intercept=0,
                         slope=1,
                         colour='gray',
                         size=0.5)
    gg <- gg+geom_step(aes(x=sensitivities,
                           y=1-specificities),
                       size=1)
    gg <- gg+geom_text(data=tset,
                       aes(x=sensitivities,
                           y=1-specificities,
                           label=thresholds),
                       hjust=1.5,
                       size=6)
    gg <- gg+scale_x_continuous('False positive rate (1-specificity)',
                                labels=percent_format(),
                                limits=c(0, 1),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous('True positive rate (sensitivity)',
                                labels=percent_format(),
                                limits=c(0, 1),
                                expand=c(0, 0))
    print(gg)
    cat('AUC:', 1-auc(afit), '\n')
}
roc_curve_plot(fit, iset)

##################################################
# Histograms of obese and non-obese persons
##################################################
histogram_plot <- function(fit, dset){
    dset <- dset %>% mutate(AGE=AGE-min(AGE))
    pset <- dset %>% mutate(Predicted=predict(fit, type='response'))
    gg <- ggplot(pset)
    gg <- gg+geom_histogram(aes(x=Predicted,
                                fill=OBESE,
                                color=OBESE),
                            position='identity',
                            alpha=0.4,
                            bins=40)
    gg <- gg+scale_fill_viridis(name='', discrete=TRUE, begin=0.8, end=0.2)
    gg <- gg+scale_colour_viridis(name='', discrete=TRUE, begin=0.8, end=0.2)
    gg <- gg+scale_x_continuous('Predicted probability',
                                labels=percent_format(),
                                limits=c(0, 1.025),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous('Frequency',
                                limits=c(0, 500),
                                expand=c(0, 0))
    print(gg)
}
histogram_plot(fit, iset)

##################################################
# Frequencies of BMI categories by year
##################################################
univariate_relative_risk <- function(dset, risk){
    rsym <- sym(risk)
    rset <- dset %>% select(BMICAT, all_of(risk))
    rset <- rset %>% distinct(BMICAT, !!rsym)
    rset <- rset %>% mutate(!!rsym:=sprintf('%4.2f', !!rsym),
                            !!rsym:=sub('-', '$-$', !!rsym))
    rset <- rset %>% spread(BMICAT, risk)
    uset <- data.frame(YEAR=switch(risk,
                                   Diff='Difference',
                                   RR='Risk ratio',
                                   OR='Odds ratio',
                                   ORM='Multi-category OR'),
                       rset,
                       stringsAsFactors=FALSE)
    return(uset)
}

bmi_frequency_table <- function(dset){
    dset <- dset %>% mutate(YEAR=sub('-', '--', YEAR))
    fset <- dset %>% group_by(YEAR, BMICAT)
    fset <- fset %>% summarize(Count=n())
    fset <- fset %>% group_by(YEAR)
    fset <- fset %>% mutate(Proportion=Count/sum(Count),
                            Percent=sprintf('%4.1f%%', 100*Count/sum(Count)),
                            Content=paste0(Count, ' (', Percent, ')'),
                            Content=sub('%', '\\\\%', Content),
                            Total=paste(sum(Count), '(100.0\\\\%)'))
    fset <- fset %>% group_by(BMICAT)
    fset <- fset %>% mutate(Diff=Proportion[YEAR == '2015--2016']-Proportion[YEAR == '1999--2000'],
                            RR=Proportion[YEAR == '2015--2016']/Proportion[YEAR == '1999--2000'],
                            ODDS=Proportion/(1-Proportion),
                            OR=ODDS[YEAR == '2015--2016']/ODDS[YEAR == '1999--2000'])
    cset <- fset %>% select(YEAR, BMICAT, Content, Total)
    cset <- cset %>% spread(BMICAT, Content)
    dset <- fset %>% univariate_relative_risk('Diff')
    rset <- fset %>% univariate_relative_risk('RR')
    oset <- fset %>% univariate_relative_risk('OR')
    mset <- fset %>% select(YEAR, BMICAT, Proportion)
    mset <- mset %>% group_by(YEAR)
    mset <- mset %>% mutate(ODDSM=Proportion/Proportion[BMICAT == 'Normal'])
    mset <- mset %>% group_by(BMICAT)
    mset <- mset %>% mutate(ORM=ODDSM[YEAR == '2015--2016']/ODDSM[YEAR == '1999--2000'])
    mset <- mset %>% univariate_relative_risk('ORM')
    xset <- bind_rows(cset, dset, rset, oset, mset)
    xset <- xset %>% select(Year=YEAR, Underweight, Normal, Overweight, Obese)
    return(xset)
}
bmi_frequency_table(iset)

##################################################
# BMI proportions by age 
##################################################
bmi_age_scatterplot <- function(dset){
    gset <- dset %>% group_by(AGE, BMICAT)
    gset <- gset %>% summarize(N=n())
    gset <- gset %>% group_by(AGE)
    gset <- gset %>% mutate(Proportion=N/sum(N))
    gset <- gset %>% ungroup(AGE)
    gg <- ggplot(gset)
    gg <- gg+geom_point(aes(x=AGE,
                            y=Proportion,
                            color=BMICAT,
                            shape=BMICAT),
                        size=2)
    for(level in levels(dset$BMICAT)){
        lset <- gset %>% filter(BMICAT == level)
        gg <- gg+geom_smooth(data=lset,
                             aes(x=AGE,
                                 y=Proportion,
                                 color=BMICAT),
                             method='lm',
                             se=FALSE)
        aset <- lset %>% filter(between(AGE, 45, 54))
        aset <- aset %>% summarize(AGE=mean(AGE),
                                   BMICAT=unique(BMICAT),
                                   Proportion=mean(Proportion))
        gg <- gg+geom_text(data=aset,
                           aes(x=AGE,
                               y=Proportion,
                               label=BMICAT,
                               color=BMICAT),
                           vjust=-1,
                           size=5)
    }
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Percent of persons',
                                labels=percent_format(accuracy=1),
                                limits=c(0, 0.5),
                                breaks=seq(0, 0.5, by=0.1),
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(discrete=TRUE, begin=0.6, end=0)
    print(gg)
}
bmi_age_scatterplot(iset)

##################################################
# Compare separate logistic regressions with
# multinomial regression
##################################################
format_multinomial <- function(fit, dset){
    cset <- coef(summary(fit))
    cset <- data.frame(BMICAT=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% gather(Predictor, Estimate, -BMICAT)
    sset <- summary(fit)$standard.errors
    sset <- data.frame(BMICAT=rownames(sset), sset, check.names=FALSE)
    sset <- sset %>% gather(Predictor, SE, -BMICAT)
    mset <- full_join(cset, sset, by=c('BMICAT', 'Predictor'))
    mset <- mset %>% transmute(BMICAT=as.character(BMICAT),
                               Predictor=sub('(Intercept)',
                                             'Intercept',
                                             Predictor,
                                             fixed=TRUE),
                               Predictor=gsub(paste0('(',
                                                     paste(names(dset), collapse='|'),
                                                     ')'),
                                              '\\1 = ',
                                              Predictor),
                               Predictor=sub(' = $', '', Predictor),
                               Coef.=Estimate,
                               SE,
                               P=2*pnorm(abs(Coef./SE), lower.tail=FALSE),
                               OR=exp(Coef.))
    return(mset)
}

format_logistic <- function(rset, dset, ref='Normal'){
    fit <- glm(BMICAT~SEX+RACE+AGE+YEAR,
               family=binomial(link='logit'),
               data=dset %>% filter(BMICAT %in% c(as.character(rset$BMICAT), ref)))
    fset <- format_model(fit, dset)
    fset <- fset %>% transmute(Predictor,
                               Coef.=Coefficient,
                               SE,
                               P=`P-value`,
                               OR=exp(Coef.))
    return(fset)
}

logistic_multinomial_table <- function(dset, reference='Normal'){
    dset <- dset %>% mutate(AGE=AGE-min(AGE),
                            YEAR=sub('-', '--', YEAR),
                            BMICAT=relevel(BMICAT, ref=reference))
    fit <- multinom(BMICAT~SEX+RACE+AGE+YEAR, data=dset)
    fset <- format_multinomial(fit, dset)
    lset <- dset %>% distinct(BMICAT)
    lset <- lset %>% filter(BMICAT != reference)
    lset <- lset %>% group_by(BMICAT)
    lset <- lset %>% do(format_logistic(., dset, ref=reference))
    xset <- full_join(lset,
                      fset,
                      by=c('BMICAT', 'Predictor'),
                      suffix=c(' (LR)', ' (MR)'))
    xset <- xset %>% ungroup()
    xset <- xset %>% mutate(BMICAT=ifelse(Predictor == 'Intercept', BMICAT, ''))
    xset <- xset %>% rename(Category='BMICAT')
    return(fit)
}
fit <- logistic_multinomial_table(iset)

##################################################
# Marginal additive effect from multinomial regression
##################################################
predict_at_level <- function(fit, dset, at){
    pset <- prediction(fit, data=dset, at=at)
    pset <- pset %>% summarize(Underweight=mean(`Pr(Underweight)`),
                               Normal=mean(`Pr(Normal)`),
                               Overweight=mean(`Pr(Overweight)`),
                               Obese=mean(`Pr(Obese)`))
    avar <- unlist(at)
    varname <- paste0('$\\text{', names(avar), '}=\\text{', avar, '}$')
    pset <- data.frame(Category=names(pset), Predicted=t(pset), row.names=NULL)
    pset %>% rename(!!sym(varname):='Predicted')
}

multinomial_additive_effect <- function(fit, dset){
    dset <- dset %>% mutate(AGE=AGE-min(AGE),
                            YEAR=sub('-', '--', YEAR))
    pold <- predict_at_level(fit, dset, list(YEAR='1999--2000'))
    pnew <- predict_at_level(fit, dset, list(YEAR='2015--2016'))
    pset <- full_join(pold, pnew, by='Category')
    pset <- pset %>% mutate('Marginal additive effect'=pset[, 3]-pset[, 2])
    if(FALSE){
        # check using margins package
        cats <- dset %>% with(levels(BMICAT))
        mset <- sapply(cats, function(x) mean(margins(fit,
                                                      variables=variable,
                                                      data=dset,
                                                      category=x)$dydx))
        print(mset)
    }
    return(pset)
}
multinomial_additive_effect(fit, iset)

