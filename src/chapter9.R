##################################################
# Create figures and tables for Chapter 9
##################################################
library(tidyverse)
library(xtable)
library(scales)
library(survey)

if(Sys.getenv('RSTUDIO') == '1')
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('shared.R')

options(survey.lonely.psu='adjust')

set.seed(12345)

##################################################
# Load MEPS 2017 data
##################################################
grab_and_curate_meps <- function(year, inscope_only=FALSE){
    dset <- grab_meps(year, inscope_only=inscope_only)
    suffix <- substr(year, 3, 4)
    inscope <- paste0('INSCOP', suffix)
    dset <- dset %>% select(Age=AGE17X,
                            Sex=SEX,
                            Race=RACETHX,
                            Diabetes=DIABDX,
                            Health=RTHLTH31,
                            Total=TOTEXP17,
                            Clusters=VARPSU,
                            Strata=VARSTR,
                            Weight=PERWT17F,
                            Inscope=!!sym(inscope))
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
                                            labels=c('No', 'Yes')),
                            Health=factor(Health,
                                          levels=seq(5),
                                          labels=c('Excellent',
                                                   'Very good',
                                                   'Good',
                                                   'Fair',
                                                   'Poor')))
    dset
}
#dset <- grab_and_curate_meps(2017)

##################################################
# Age distributions and proportion of persons in
# excellent health using stratified versus simple
# random sampling
##################################################
format_estimates <- function(vec, fmt='%4.2f'){
    paste0(sprintf(fmt, mean(vec)),
           ' (',
           sprintf(fmt, quantile(vec, prob=0.025)),
           ', ',
           sprintf(fmt, quantile(vec, prob=0.975)),
           ')')
}

stratified_vs_simple_outcome <- function(dset, b=100, size=1000, fmt='%4.2f'){
    # define age strata and outcome
    dset <- dset %>% mutate(Agegroup=cut(Age,
                                         breaks=c(0, 25, 55, Inf),
                                         right=FALSE),
                            Excellent=Health == 'Excellent')
    # define weights for age strata
    prop <- dset %>% with(c(prop.table(table(Agegroup))))
    # simple random sampling
    simprs_once <- function(dset, size){
        y <- dset %>% sample_n(size)
        with(y, mean(Excellent, na.rm=TRUE))
    }
    simprs_many <- function(times, dset, size){
        replicate(times, simprs_once(dset=dset, size=size))
    }
    simprs <- simprs_many(b, dset, size)
    # age-stratified sampling in each bootstrap sample
    strars_once <- function(dset, size){
        y <- dset %>% slice(stratsample(dset$Agegroup, size*prop))
        with(y, mean(Excellent, na.rm=TRUE))
    }
    strars_many <- function(times, dset, size){
        replicate(times, strars_once(dset=dset, size=size))
    }
    strars <- strars_many(b, dset, size)
    sset <- tibble(Design=c('Simple random sampling',
                            'Age-stratified sampling'),
                   Proportion=c(format_estimates(simprs, fmt=fmt),
                                format_estimates(strars, fmt=fmt)))
    print(sset)
}
#stratified_vs_simple_outcome(dset)

stratified_vs_simple_distribution <- function(dset,
                                              b=100,
                                              size=1000,
                                              fmt='%4.2f',
                                              saveit=FALSE){
    # define age strata and outcome
    dset <- dset %>% mutate(Agegroup=cut(Age,
                                         breaks=c(0, 25, 55, Inf),
                                         labels=c('<25', '25-54', '>54'),
                                         right=FALSE))
    # define weights for age strata
    prop <- dset %>% with(c(prop.table(table(Agegroup))))
    # simple random sampling
    simprs_once <- function(dset, size){
        y <- dset %>% sample_n(size)
        as.data.frame(prop.table(table(y$Agegroup, dnn='Agegroup')))
    }
    simprs <- data.frame(Sample=seq(b))
    simprs <- simprs %>% group_by(Sample)
    simprs <- simprs %>% do(simprs_once(dset, size))
    # age-stratified sampling in each bootstrap sample
    strars_once <- function(dset, size){
        y <- dset %>% slice(stratsample(dset$Agegroup, size*prop))
        as.data.frame(prop.table(table(y$Agegroup, dnn='Agegroup')))
    }
    strars <- data.frame(Sample=seq(b))
    strars <- strars %>% group_by(Sample)
    strars <- strars %>% do(strars_once(dset, size))
    # merge and format table
    mset <- full_join(simprs,
                      strars,
                      by=c('Sample', 'Agegroup'),
                      suffix=c('.simprs', '.strars'))
    mset <- mset %>% group_by(Agegroup)
    mset <- mset %>% summarize('Simple random sampling'=format_estimates(Freq.simprs, fmt),
                               'Age-stratified sampling'=format_estimates(Freq.strars, fmt))
    mset <- mset %>% rename(Age='Agegroup')
    mset <- mset %>% mutate(Age=sub('25-54', '25--54', Age))
    if(saveit){
        filename <- '09-stratified_vs_simple.tex'
        Caption <- paste('Proportions and 95\\% confidence intervals of MEPS 2017 participants in various age groups selected by repeated simple random sampling or repeated age-stratified sampling; based on', b, 'replicates of', size, 'participants each.')
        print(xtable(mset,
                     digits=0,
                     align='llcc',
                     label='tab:stratified_vs_simple',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              sanitize.text.function=sanitize,
              include.rownames=FALSE,
              hline.after=0)
    }
}
#stratified_vs_simple_distribution(dset, saveit=TRUE)

##################################################
# Analyze costs associated with prior diagnosis of
# diabetes using MEPS data 
##################################################
format_model <- function(fit, type, fmt='%4.2f'){
    cset <- coef(summary(fit))
    cset <- data.frame(Predictor=rownames(cset),
                       cset,
                       check.names=FALSE,
                       row.names=NULL)
    cset <- cset %>% filter(Predictor == 'DiabetesYes')
    cset <- cset %>% transmute(Model=type,
                               Coefficient=Estimate,
                               '95% CI'=paste0('(',
                                               sprintf(fmt, Estimate-1.96*`Std. Error`),
                                               ', ',
                                               sprintf(fmt, Estimate+1.96*`Std. Error`),
                                               ')'),
                               'P-value'=sub('(.*)', '$\\1 $', format_pvalue(`Pr(>|t|)`)))
    cset
}

diabetes_analysis <- function(dset, saveit=FALSE){
    cat('Adults in sample:',
        adults <- dset %>% filter(Age > 17) %>% nrow(),
        '\n')
    cat('Persons with diabetes in sample:',
        diabetics <- dset %>% filter(Age > 17, Diabetes == 'Yes') %>% nrow(),
        paste0('(',
               sprintf('%3.1f%%', 100*diabetics/adults),
               ')'),
        '\n')
    # specify survey designs
    weighted_design <- svydesign(id=~Clusters,
                                 strata=~Strata,
                                 weights=~Weight,
                                 data=dset,
                                 nest=TRUE)
    unweighted_design <- svydesign(id=~Clusters,
                                   strata=~Strata,
                                   weights=NULL,
                                   data=dset,
                                   nest=TRUE)
    # number of persons with prior diagnosis of diabetes
    persons <- svytotal(~Diabetes,
                        design=subset(weighted_design, Inscope == 1 & Age > 17),
                        na.rm=TRUE)
    persons <- as.data.frame(persons)
    cat('Persons with prior diabetes diagnosis:',
        persons['DiabetesYes', 'total'],
        paste0('(',
               sprintf('%3.1f%%',
                       100*persons['DiabetesYes', 'total']/sum(persons[, 'total'])),
               ')'),
        '\n')
    # estimate incremental total expenditures
    equation <- Total~Age+Sex+Race+Diabetes
    naive_fit <- lm(as.formula(equation), data=dset, subset=Inscope == 1 & Age > 17)
    unweighted_fit <- svyglm(as.formula(equation),
                             design=subset(unweighted_design, Inscope == 1 & Age > 17))
    weighted_fit <- svyglm(as.formula(equation),
                           design=subset(weighted_design, Inscope == 1 & Age > 17))
    fset <- bind_rows(format_model(naive_fit, 'Naive'),
                      format_model(unweighted_fit, 'Unweighted'),
                      format_model(weighted_fit, 'Weighted'))
    if(saveit){
        filename <- '09-incremental_cost_diabetes.tex'
        Caption <- 'Estimated incremental medical costs for persons age 18 years or older in the US using a naive linear regression and unweighted and weighted regressions that account for the MEPS 2017 survey design, each adjusted for age, sex and race/ethnicity.'
        print(xtable(fset,
                     digits=2,
                     align='llccc',
                     label='tab:incremental_cost_diabetes',
                     caption=Caption),
              file=here('tables', filename),
              table.placement='!ht',
              caption.placement='bottom',
              sanitize.colnames.function=sanitize,
              sanitize.text.function=identity,
              include.rownames=FALSE,
              hline.after=0)
    }
    fset <- fset %>% mutate(Total=Coefficient*persons['DiabetesYes', 'total'])
    print(fset)
}
#diabetes_analysis(dset, saveit=TRUE)

