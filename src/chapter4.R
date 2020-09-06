##################################################
# Create figures and tables for Chapter 4
##################################################
library(tidyverse)
library(scales)
library(grid)
library(viridis)
library(xtable)
library(here)
library(margins)
library(pROC)
library(nnet)

if(Sys.getenv('RSTUDIO') == '1')
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('shared.R')

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
obese_frequency_table <- function(dset, saveit=FALSE){
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
    if(saveit){
        filename <- '04-obese_year.tex'
        Caption <- 'Frequency of obese and non-obese persons aged 20--59 years in NHANES sample data for 1999--2000 and 2015--2016.'
        print(xtable(fset,
                     digits=0,
                     align='lc|cc|c',
                     label='tab:obese_year',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              hline.after=0)
    }
}
#obese_frequency_table(iset, saveit=TRUE)

##################################################
# Compare OR versus RR
##################################################
or_versus_rr_table <- function(saveit=FALSE){
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
    if(saveit){
        filename <- '04-or_versus_rr.tex'
        Caption <- 'Risk ratio (RR) and odds ratio (OR) for selected probabilities of positive outcomes.'
        print(xtable(dset,
                     digits=2,
                     align='lcccc',
                     label='tab:or_versus_rr',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.colnames.function=identity,
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              hline.after=0)
    }
}
#or_versus_rr_table(saveit=TRUE)

##################################################
# Proportion of obese persons by age 
##################################################
count_scatterplot <- function(dset, ext='pdf', saveit=FALSE){
    gg_theme()
    gg <- ggplot(dset, aes(x=AGE, y=OBESE))
    gg <- gg+geom_count(color='gray')
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_discrete(name='')
    print(gg)
    if(saveit){
        filename <- '04-bmi_age'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
}
#count_scatterplot(iset, saveit=TRUE)

proportion_scatterplot <- function(dset, ext='pdf', saveit=FALSE){
    dset <- dset %>% mutate(OBESE=as.integer(OBESE)-1)
    gset <- dset %>% group_by(AGE)
    gset <- gset %>% summarize(Proportion=sum(OBESE)/n())
    gg_theme()
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
    if(saveit){
        filename <- '04-proportion_obese'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
}
#proportion_scatterplot(iset, saveit=TRUE)

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
                               'P-value'=format_pvalue(`Pr(>|z|)`))
    cset
}

logistic_regression_table <- function(dset, saveit=FALSE){
    fit <- glm(OBESE~YEAR,
               family=binomial(link='logit'),
               data=dset)
    fset <- format_model(fit, dset)
    if(saveit){
        filename <- '04-obese_year_regression.tex'
        Caption <- 'Logistic regression of obese persons aged 20--59 years in NHANES sample data for 1999--2000 and 2015--2016.'
        print(xtable(fset,
                     digits=3,
                     align='llccc',
                     label='tab:obese_year_regression',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              hline.after=0)
    }
}
#logistic_regression_table(iset, saveit=TRUE)

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
    uset %>% select(Variable,
                    'Sample, n'=Total,
                    'Obese, n (%)'=Content.obese,
                    'Odds ratio'=OR)
}

univariate_table <- function(dset, saveit=FALSE){
    dset <- dset %>% mutate(YEAR=sub('-', '--', YEAR))
    sset <- univariate_summary(dset, 'SEX', 'Male')
    rset <- univariate_summary(dset, 'RACE', 'Non-Hispanic White')
    yset <- univariate_summary(dset, 'YEAR', '1999--2000')
    uset <- bind_rows(sset, rset, yset)
    if(saveit){
        filename <- '04-univariate_summary.tex'
        Caption <- 'Univariate summaries of selected categorical variables for obese persons aged 20--59 years in NHANES sample data in 1999--2000 and 2015--2016.'
        print(xtable(uset,
                     digits=2,
                     align='llccc',
                     label='tab:univariate_summary',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.colnames.function=sanitize,
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              hline.after=c(0, cumsum(c(nrow(sset),
                                        nrow(rset)))))
    }
}
#univariate_table(iset, saveit=TRUE)

##################################################
# Multivariate logistic regression
##################################################
multivariate_table <- function(dset, alpha=0.05, saveit=FALSE){
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
    if(saveit){
        filename <- '04-multivariate_regression.tex'
        Caption <- 'Multivariate logistic regression of being obese for persons age 20--59 years in NHANES sample data in 1999--2000 and 2015--2016 showing each estimated coefficient ($\\hat{\\beta}$), standard error (SE), P-value for the hypothesis of no association ($\\hat{\\beta}=0$), odds ratio ($\\text{OR}=\\exp(\\hat{\\beta})$) and 95\\% confidence interval (CI) for the OR ($\\exp(\\hat{\\beta} \\pm 1.96\\times \\text{SE})$).'
        print(xtable(fset,
                     digits=2,
                     align='llccccc',
                     label='tab:multivariate_regression',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.colnames.function=sanitize,
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              hline.after=0)
    }
    fit
}
#fit <- multivariate_table(iset, saveit=TRUE)

##################################################
# Predicted probabilities of obesity by sex and year
# restricted to RACE = Other/Mixed
##################################################
prediction_plot <- function(fit, dset, ext='pdf', saveit=FALSE){
    minage <- min(dset$AGE)
    dset <- dset %>% mutate(AGE=AGE-minage)
    nset <- dset %>% with(expand.grid(SEX=c('Male', 'Female'),
                                      RACE='Other/Mixed',
                                      AGE=seq(min(AGE), max(AGE)),
                                      YEAR=c('1999-2000', '2015-2016')))
    nset <- nset %>% mutate(PHAT=predict(fit, newdata=nset, 'response'),
                            AGE=AGE+minage)
    nset <- nset %>% group_by(SEX, AGE)
    gg_theme()
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
    if(saveit){
        filename <- '04-predicted_probabilities'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
}
#prediction_plot(fit, iset, saveit=TRUE)

##################################################
# Compare differences, RRs, and ORs for obesity by
# sex and year restricted to RACE = Other/Mixed
##################################################
prediction_comparison <- function(fit, dset, ext='pdf', comparison=FALSE, saveit=FALSE){
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
    gg_theme()
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
    if(saveit){
        filename <- '04-prediction_comparison'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
}
#prediction_comparison(fit, iset, saveit=TRUE)

##################################################
# Marginal additive effect of year
##################################################
marginal_additive_effect <- function(fit, dset, variable){
    mset <- margins(fit, data=dset, variables=variable)
    print(summary(mset))
}
#marginal_additive_effect(fit, iset, 'YEAR')

lm_check <- function(dset, vars){
    dset <- dset %>% mutate(AGE=AGE-min(AGE),
                            OBESE=as.integer(OBESE)-1)
    lmod <- lm(OBESE~SEX+RACE+AGE+YEAR, data=dset)
    ests <- cbind(coef(lmod), confint(lmod))
    print(ests[grepl(vars, rownames(ests)), ])
}
#lm_check(iset, 'YEAR')

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
    vset %>% mutate(AIC=AIC(fit), BIC=BIC(fit))
}

aicbic_table <- function(dset, saveit=FALSE){
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
    if(saveit){
        filename <- '04-aic_bic.tex'
        Caption <- 'AIC and BIC for logistic regressions with year and all possible combinations of sex, race and age for persons age 20--59 years in NHANES sample data for years 1999--2000 and 2015--2016.'
        print(xtable(vset,
                     align='lccccccc',
                     digits=1,
                     label='tab:aic_bic',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              include.rownames=FALSE,
              hline.after=0)
    }
}
#aicbic_table(iset, saveit=TRUE)

##################################################
# Hosmer-Lemeshow comparison of observed and predicted
##################################################
hosmer_lemeshow_table <- function(fit, dset, ngroups=10, saveit=FALSE){
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
    if(saveit){
        filename <- '04-hosmer_lemeshow.tex'
        Caption <- 'Observed (Obs) and predicted (Pred) numbers of obese and non-obese persons by decile of predicted probabilities and corresponding Pearson residuals (Res) based on the model in Table~\\ref{tab:multivariate_regression}.'
        print(xtable(gset,
                     align='lccccccc',
                     digits=1,
                     label='tab:hosmer_lemeshow',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              sanitize.colnames.function=identity,
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              include.rownames=FALSE,
              hline.after=0)
    }
}
#hosmer_lemeshow_table(fit, iset, saveit=TRUE)

##################################################
# Prediction accuracy 
##################################################
prediction_accuracy_table <- function(fit, dset, threshold=0.5, saveit=FALSE){
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
    if(saveit){
        filename <- '04-prediction_accuracy.tex'
        Caption <- 'Frequency of observed (Obs) and predicted (Pred) obese and non-obese persons based on the model in Table~\\ref{tab:multivariate_regression} using a threshold of $T=0.5$ for predictions.'
        print(xtable(fset,
                     digits=0,
                     align='lc|cc|c',
                     label='tab:prediction_accuracy',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.colnames.function=identity,
              sanitize.text.function=identity,
              math.style.negative=TRUE,
              hline.after=0)
    }
}
#prediction_accuracy_table(fit, iset, saveit=TRUE)

##################################################
# ROC curve
##################################################
roc_curve_plot <- function(fit, dset, ext='pdf', saveit=FALSE){
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
    gg_theme(aspect.ratio=1)
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
    if(saveit){
        filename <- '04-roc_curve'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               here('figures', filename),
               height=5,
               width=8)
    }
    cat('AUC:', 1-auc(afit), '\n')
}
#roc_curve_plot(fit, iset, saveit=TRUE)

##################################################
# Histograms of obese and non-obese persons
##################################################
histogram_plot <- function(fit, dset, ext='pdf', saveit=FALSE){
    dset <- dset %>% mutate(AGE=AGE-min(AGE))
    pset <- dset %>% mutate(Predicted=predict(fit, type='response'))
    gg_theme(legend.position=c(0.8, 0.8))
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
    if(saveit){
        filename <- '04-prediction_histograms'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               here('figures', filename),
               height=5,
               width=10)
    }
}
#histogram_plot(fit, iset, saveit=TRUE)

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
    data.frame(YEAR=switch(risk,
                           Diff='Difference',
                           RR='Risk ratio',
                           OR='Odds ratio',
                           ORM='Multi-category OR'),
               rset,
               stringsAsFactors=FALSE)
}

bmi_frequency_table <- function(dset, saveit=FALSE){
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
    #xset <- xset %>% select(Year=YEAR, Underweight, Normal, Overweight, Obese, Total)
    xset <- xset %>% select(Year=YEAR, Underweight, Normal, Overweight, Obese)
    if(saveit){
        filename <- '04-bmi_year.tex'
        Caption <- 'Frequency of persons in BMI categories aged 20--59 years in NHANES sample data for 1999--2000 and 2015--2016 and measures of relative change.'
        print(xtable(xset,
                     digits=0,
                     #align='lc|cccc|c',
                     align='lc|cccc',
                     label='tab:bmi_year',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.text.function=identity,
              math.style.negative=TRUE,
              hline.after=c(0, nrow(cset), nrow(xset)-1))
    }
}
#bmi_frequency_table(iset, saveit=TRUE)

##################################################
# BMI proportions by age 
##################################################
bmi_age_scatterplot <- function(dset, ext='pdf', saveit=FALSE){
    gset <- dset %>% group_by(AGE, BMICAT)
    gset <- gset %>% summarize(N=n())
    gset <- gset %>% group_by(AGE)
    gset <- gset %>% mutate(Proportion=N/sum(N))
    gset <- gset %>% ungroup(AGE)
    gg_theme()
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
    if(saveit){
        filename <- '04-multicategory_bmi_age'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               here('figures', filename),
               height=5,
               width=10)
    }
}
#bmi_age_scatterplot(iset, saveit=TRUE)

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
                               P=format_pvalue(2*pnorm(abs(Coef./SE), lower.tail=FALSE)),
                               OR=exp(Coef.))
    mset
}

format_logistic <- function(rset, dset, ref='Normal'){
    fit <- glm(BMICAT~SEX+RACE+AGE+YEAR,
               family=binomial(link='logit'),
               data=dset %>% filter(BMICAT %in% c(as.character(rset$BMICAT), ref)))
    fset <- format_model(fit, dset)
    fset %>% transmute(Predictor,
                       Coef.=Coefficient,
                       SE,
                       P=`P-value`,
                       OR=exp(Coef.))
}

logistic_multinomial_table <- function(dset, reference='Normal', saveit=FALSE){
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
    if(saveit){
        filename <- '04-logistic_multinomial.tex'
        Caption <- 'Separate logistic regressions (LR) and multinomial regression (MR) of BMI categories for persons age 20--59 years in NHANES sample data for 1999--2000 and 2015--2016.'
        print(xtable(xset,
                     digits=2,
                     align='lll|cccc|cccc',
                     label='tab:logistic_multinomial',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.text.function=sanitize,
              math.style.negative=TRUE,
              hline.after=c(0, 8, 16))
    }
    fit
}
#fit <- logistic_multinomial_table(iset, saveit=TRUE)

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

multinomial_additive_effect <- function(fit, dset, saveit=FALSE){
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
    if(saveit){
        filename <- '04-multinomial_additive_effect.tex'
        Caption <- 'The marginal additive effect of YEAR calculated as the mean difference of predicted probabilities if all persons were in the indicated BMI category in 1999--2000 and in 2015--2016 based on the multinomial model in Table~\\ref{tab:logistic_multinomial}.'
        Label <- 'multinomial_additive_effect'
        print(xtable(pset,
                     align='llccc',
                     digits=3,
                     label='tab:multinomial_additive_effect',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              include.rownames=FALSE,
              sanitize.text.function=identity,
              math.style.negative=TRUE,
              hline.after=0)
    }
}
#multinomial_additive_effect(fit, iset, saveit=TRUE)

