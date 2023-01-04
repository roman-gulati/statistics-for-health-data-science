##################################################
# Chapter 8 examples
##################################################
library(tidyverse)
library(scales)
library(MatchIt)
library(viridis)

source('grab_meps.R')

set.seed(12345)

##################################################
# Kidney stone example
##################################################
kset <- tibble(Treatment=factor(c('Open surgery',
                                  'Open surgery',
                                  'PCNL',
                                  'PCNL')),
               Size=factor(c('Small ($<$2 cm)',
                             'Large ($\\geq$2 cm)',
                             'Small ($<$2 cm)',
                             'Large ($\\geq$2 cm)'),
                           levels=c('Small ($<$2 cm)',
                                    'Large ($\\geq$2 cm)')),
               Patients=c(87, 263, 270, 80),
               Successes=c(81, 192, 234, 55),
               Total=Successes/Patients)

kidney_example_table <- function(dset){
    tset <- dset %>% group_by(Treatment) %>% summarize(Patients=sum(Patients),
                                                       Successes=sum(Successes))
    sset <- dset %>% group_by(Size) %>% summarize(Patients=sum(Patients),
                                                  Successes=sum(Successes))
    oset <- dset %>% summarize(Patients=sum(Patients),
                               Successes=sum(Successes))
    dtab <- bind_rows(dset, tset, sset, oset)
    dtab <- dtab %>% transmute(Treatment=fct_explicit_na(Treatment, na_level='Total'),
                               Size=fct_explicit_na(Size, na_level='Total'),
                               Total=paste(paste(Successes, Patients, sep='/'),
                                           paste0('(',
                                                  sprintf('%4.1f%%', 100*Successes/Patients),
                                                  ')')))
    dtab <- dtab %>% spread(Size, Total)
    return(dtab)
}
kidney_example_table(kset)

stratified_kidney_analysis <- function(dset){
    sset <- dset %>% group_by(Size)
    sset <- sset %>% summarize(Proportion=sum(Patients)/sum(dset$Patients),
                               Difference=-diff(Total))
    cat('Proportions and success rates by kidney size:\n')
    print(sset)
    sset <- full_join(dset, sset)
    sset <- sset %>% mutate(Weighted=Proportion*Total)
    sset <- sset %>% group_by(Treatment)
    sset <- sset %>% summarize(Mean=sum(Weighted))
    cat('Weighted success rates by treatment group:\n')
    print(sset)
    sset <- sset %>% ungroup()
    sset <- sset %>% summarize(Effect=-diff(Mean))
    cat('Difference in weighted success rates:\n')
    print(sset)
    fit <- lm(Total~Treatment+Size, data=dset)
    cat('Fit from linear regression:\n')
    print(coef(fit))
    dset <- dset %>% select(-Total)
    eset <- dset[rep(1:nrow(dset), dset$Patients), ]
    eset <- eset %>% group_by(Treatment, Size)
    eset <- eset %>% mutate(Successes=c(rep(1, unique(Successes)),
                                        rep(0, unique(Patients)-unique(Successes))))
    fit <- glm(Successes~Treatment+Size, data=eset, family=binomial(link='logit'))
    cat('Fit from logistic regression:\n')
    print(exp(coef(fit)))
}
stratified_kidney_analysis(kset)

matched_kidney_analysis <- function(dset, nreps=100){
    dset <- dset %>% select(-Total)
    # expand dset from counts to individual records
    eset <- dset[rep(1:nrow(dset), dset$Patients), ]
    eset <- eset %>% group_by(Treatment, Size)
    eset <- eset %>% mutate(Successes=c(rep(1, unique(Successes)),
                                        rep(0, unique(Patients)-unique(Successes))))
    # determine smaller number of patients in each treatment group
    nset <- dset %>% select(Treatment, Target=Patients)
    nset <- nset %>% group_by(Treatment)
    nset <- nset %>% filter(Target == min(Target))
    # merge individual records with sample target size
    mset <- full_join(eset, nset, by='Treatment')
    for(index in 1:nreps){
        # sample target size of individual records
        rset <- mset %>% group_by(Treatment, Size)
        rset <- rset %>% sample_n(size=unique(Target))
        # calculate success rates for each treatment
        rset <- rset %>% group_by(Treatment)
        rset <- rset %>% summarize(Rate=sum(Successes)/n())
        rset <- rset %>% mutate(Rep=index)
        if(index == 1)
            sset <- rset
        else
            sset <- bind_rows(sset, rset)
    }
    sset <- sset %>% group_by(Treatment)
    sset <- sset %>% summarize(Rate=mean(Rate))
    cat('Matched analysis success rates:\n')
    print(sset)
    cat('Difference in success rates:', with(sset, -diff(Rate)), '\n')
}
matched_kidney_analysis(kset)

weighted_kidney_analysis <- function(dset){
    dset <- dset %>% mutate(Treatment=as.character(Treatment), Total=NULL)
    dset <- dset %>% group_by(Size)
    # estimate probability of open surgery
    oset <- dset %>% summarize(Ratio=paste(Patients[Treatment == 'Open surgery'],
                                           sum(Patients),
                                           sep='/'),
                               Probability=Patients[Treatment == 'Open surgery']/sum(Patients))
    oset <- oset %>% mutate(Treatment='Open surgery')
    # estimate probability of PCNL
    pset <- dset %>% summarize(Ratio=paste(Patients[Treatment == 'PCNL'],
                                           sum(Patients),
                                           sep='/'),
                               Probability=Patients[Treatment == 'PCNL']/sum(Patients))
    pset <- pset %>% mutate(Treatment='PCNL')
    # merge with aggregate counts
    tset <- bind_rows(oset, pset)
    sset <- full_join(dset, tset, by=c('Size', 'Treatment'))
    # calculate weighted numbers of patients and successes
    sset <- sset %>% mutate(Weighted.Patients=Patients/Probability,
                            Weighted.Successes=Successes/Probability)
    # calculate success rates and treatment difference
    sset <- sset %>% group_by(Treatment)
    sset <- sset %>% summarize(Weighted.Rate=sum(Weighted.Successes)/sum(Weighted.Patients))
    cat('Weighted analysis success rates:\n')
    print(sset)
    cat('Difference in success rates:', with(sset, -diff(Weighted.Rate)), '\n')
}
weighted_kidney_analysis(kset)

##################################################
# Arthritis and health expenditures example
##################################################
exclude <- function(dset,
                    variable,
                    missing_codes=c('not ascertained',
                                    'don\'t know',
                                    'refused',
                                    'inapplicable')){
    varsym <- sym(variable)
    cat('Remove',
        dset %>% filter(!!varsym %in% missing_codes) %>% nrow(),
        'missing',
        variable,
        '\n')
    dset <- dset %>% filter(!(!!varsym %in% missing_codes))
    return(dset)
}

grab_and_curate_meps <- function(year){
    dset <- grab_meps(year)
    dset <- dset %>% select(Arthritis=ARTHDX,
                            Age=AGE17X,
                            Sex=SEX,
                            Cost=TOTEXP17,
                            Pain=ADPAIN42,
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
                            Pain=factor(Pain,
                                        levels=c(-9, -1, 1, 2, 3, 4, 5),
                                        labels=c('not ascertained',
                                                 'inapplicable',
                                                 'not at all',
                                                 'a little bit',
                                                 'moderately',
                                                 'quite a bit',
                                                 'extremely')),
                            Diabetes=factor(Diabetes,
                                            levels=c(-9, -8, -7, -1, 2, 1),
                                            labels=c('not ascertained',
                                                     'don\'t know',
                                                     'refused',
                                                     'inapplicable',
                                                     'no',
                                                     'yes')),
                            Anycost=Cost > 0)
    cat('Remove',
        dset %>% filter(Age <= 17) %>% nrow(),
        'age 17 or younger\n')
    dset <- dset %>% filter(Age > 17)
    dset <- dset %>% exclude('Arthritis')
    dset <- dset %>% exclude('Pain')
    dset <- dset %>% exclude('Diabetes')
    dset <- dset %>% droplevels()
    return(dset)
}
dset <- grab_and_curate_meps(2017)

format_model <- function(fit, lab=NULL){
    if(is.null(lab))
        lab <- switch(family(fit)$family,
                      binomial=' (logistic)',
                      gaussian=' (linear)')
    cset <- coef(summary(fit))
    cset <- data.frame(Variable=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% transmute(Variable=sub('(Intercept)', 'Intercept', Variable, fixed=TRUE),
                               Variable=gsub(paste0('(', paste(names(dset), collapse='|'), ')'), '\\1 = ', Variable),
                               Variable=sub('Age =', 'Age', Variable),
                               Variable=sub('Pain =', 'Pain', Variable),
                               Coef.=Estimate,
                               SE=`Std. Error`)
    names(cset)[-1] <- sub('$', lab, names(cset)[-1])
    return(cset)
}

stratified_arthritis_analysis <- function(dset,
                                          logistic_response='Anycost',
                                          linear_response='Cost',
                                          predictor=c('Arthritis',
                                                      'Sex',
                                                      'Diabetes',
                                                      'Age')){
    linear_predictor <- paste(predictor, collapse='+')
    linear_equation <- paste(linear_response, linear_predictor, sep='~')
    logistic_equation <- paste(logistic_response, linear_predictor, sep='~')
    logistic_fit <- glm(as.formula(logistic_equation),
                        data=dset,
                        family=binomial(link='logit'))
    linear_fit <- lm(as.formula(linear_equation),
                     data=dset,
                     subset=Anycost)
    linear_formatted <- format_model(linear_fit)
    logistic_formatted <- format_model(logistic_fit)
    fset <- bind_cols(logistic_formatted,
                      linear_formatted %>% select(-Variable))
    yset <- dset %>% mutate(Arthritis='yes')
    yset <- yset %>% mutate(Probability=predict(logistic_fit,
                                                newdata=yset,
                                                type='response'),
                            Amount=predict(linear_fit,
                                           newdata=yset,
                                           type='response'),
                            Scaled=Probability*Amount)
    nset <- dset %>% mutate(Arthritis='no')
    nset <- nset %>% mutate(Probability=predict(logistic_fit,
                                                newdata=nset,
                                                type='response'),
                            Amount=predict(linear_fit,
                                           newdata=nset,
                                           type='response'),
                            Scaled=Probability*Amount)
    cat('Mean health expenditure if all persons had arthritis (conditional): $',
        ycon <- yset %>% with(mean(Scaled)),
        '\n',
        sep='')
    cat('Mean health expenditure if no persons had arthritis (conditional): $',
        ncon <- nset %>% with(mean(Scaled)),
        '\n',
        sep='')
    cat('Difference from conditional analysis: $', ycon-ncon, '\n', sep='')
    cat('Mean health expenditure if all persons had arthritis (unconditional): $',
        yunc <- dset %>% filter(Arthritis == 'yes') %>% with(mean(Cost)),
        '\n',
        sep='')
    cat('Mean health expenditure if no persons had arthritis (unconditional): $',
        nunc <- dset %>% filter(Arthritis == 'no') %>% with(mean(Cost)),
        '\n',
        sep='')
    cat('Difference from conditional analysis: $', yunc-nunc, '\n', sep='')
    uset <- dset %>% group_by(Arthritis)
    uset <- uset %>% summarize(Cost=mean(Cost))
    uset <- uset %>% ungroup()
    uset <- uset %>% mutate(Arthritis=as.character(Arthritis),
                            Method='No conditioning')
    sset <- bind_rows(yset, nset)
    sset <- sset %>% group_by(Arthritis)
    sset <- sset %>% summarize(Cost=mean(Scaled))
    sset <- sset %>% mutate(Method='Stratifying')
    bset <- bind_rows(uset, sset)
    return(bset)
}
stratified_arthritis_analysis(dset)

matched_arthritis_analysis <- function(dset, max_age_difference=10){
    # matchit requires 0/1 response
    dset <- dset %>% mutate(Arthritis=as.integer(Arthritis)-1)
    mfit <- matchit(Arthritis~Age+Diabetes+Sex,
                    data=dset,
                    exact=c('Sex', 'Diabetes'))
    if(max_age_difference > 0){
        cases <- match.data(mfit, group='treat')
        controls <- match.data(mfit, group='control')
        cases <- cases %>% arrange(Sex, Diabetes, Age)
        controls <- controls %>% arrange(Sex, Diabetes, Age)
        age_difference <- abs(cases$Age-controls$Age)
        cases <- cases %>% filter(age_difference <= max_age_difference)
        controls <- controls %>% filter(age_difference <= max_age_difference)
        cset <- bind_rows(cases, controls)
    } else {
        cset <- match.data(mfit, group='all')
    }
    cset <- cset %>% group_by(Arthritis)
    sset <- cset %>% summarize(Cost=mean(Cost))
    sset <- sset %>% mutate(Arthritis=ifelse(Arthritis == 0, 'no', 'yes'))
    cat('Mean health expenditures in matched groups:\n')
    print(sset)
    cat('Difference in mean health expenditures rates:', with(sset, diff(Cost)), '\n')
    sset <- sset %>% mutate(Method='Matching')
    return(sset)
}
matched_arthritis_analysis(dset)

stratified_propensity_analysis <- function(dset, nstrata=4){
    sset <- dset %>% mutate(Strata=cut_interval(Propensity, n=nstrata))
    sset <- sset %>% group_by(Strata)
    pset <- sset %>% summarize(Proportion=n()/nrow(sset))
    sset <- sset %>% do(format_model(lm(Cost~Arthritis, data=.), lab=''))
    mset <- sset %>% full_join(pset, by='Strata')
    mset <- mset %>% mutate(Weighted=Proportion*Coef.)
    mset <- mset %>% group_by(Variable)
    mset <- mset %>% summarize(Coef.=sum(Weighted))
    mset <- mset %>% rename(Arthritis='Variable')
    mset <- mset %>% transmute(Arthritis=sub('Intercept', 'no', Arthritis),
                               Arthritis=sub('Arthritis = yes', 'yes', Arthritis),
                               Cost=ifelse(Arthritis == 'yes', sum(Coef.), Coef.))
    cat('Mean health expenditures across', nstrata, 'strata:\n')
    print(mset)
    cat('Difference in mean health expenditures rates:', with(mset, -diff(Cost)), '\n')
    mset <- mset %>% mutate(Method='Stratifying')
    return(mset)
}

weighted_propensity_analysis <- function(dset){
    uset <- dset %>% do(format_model(lm(Cost~Arthritis, data=.), lab=''))
    uset <- uset %>% mutate(Method='No conditioning')
    wset <- dset %>% mutate(Weight=ifelse(Arthritis == 'yes',
                                          1/Propensity,
                                          1/(1-Propensity)))
    wset <- wset %>% do(format_model(lm(Cost~Arthritis, data=., weights=Weight), lab=''))
    wset <- wset %>% mutate(Method='Weighting')
    sset <- bind_rows(uset, wset)
    sset <- sset %>% rename(Arthritis='Variable')
    sset <- sset %>% group_by(Method)
    sset <- sset %>% mutate(Arthritis=sub('Intercept', 'no', Arthritis),
                            Arthritis=sub('Arthritis = yes', 'yes', Arthritis),
                            Cost=ifelse(Arthritis == 'yes', sum(Coef.), Coef.))
    sset <- sset %>% ungroup()
    sset <- sset %>% transmute(Arthritis, Cost, Method)
    return(sset)
}

propensity_score_estimation <- function(dset){
    fit <- glm(Arthritis~Sex+Diabetes+Age,
               data=dset,
               family=binomial(link='logit'))
    pset <- dset %>% mutate(Propensity=fitted(fit))
    gg <- ggplot(pset)
    gg <- gg+geom_histogram(aes(x=Propensity,
                                y=..density..,
                                colour=Arthritis,
                                fill=Arthritis),
                            binwidth=0.05,
                            alpha=0.3,
                            position='identity')
    gg <- gg+theme(legend.position=c(0.95, 0.95),
                   legend.justification=c('right', 'top'))
    gg <- gg+scale_x_continuous(name='Propensity score',
                                limits=c(0, 1.02),
                                breaks=seq(0, 1, by=0.2),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Density',
                                limits=c(0, 4),
                                breaks=seq(0, 4, by=1),
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(name='Arthritis', discrete=TRUE, begin=0.75, end=0.25)
    gg <- gg+scale_fill_viridis(name='Arthritis', discrete=TRUE, begin=0.75, end=0.25)
    print(gg)
    return(pset)
}
pset <- propensity_score_estimation(dset)

propensity_score_analysis <- function(dset, nstrata=10){
    sset <- stratified_propensity_analysis(dset, nstrata=nstrata)
    mset <- matched_arthritis_analysis(dset)
    wset <- weighted_propensity_analysis(dset)
    tset <- bind_rows(sset, mset, wset)
    tset <- tset %>% spread(Arthritis, Cost)
    tset <- tset %>% mutate(Method=factor(Method, levels=c('No conditioning',
                                                           'Stratifying',
                                                           'Matching',
                                                           'Weighting')))
    tset <- tset %>% arrange(Method)
    tset <- tset %>% transmute(Method,
                               'Cost (no arthritis)'=dollar_format(accuracy=1)(no),
                               'Cost (arthritis)'=dollar_format(accuracy=1)(yes),
                               'Cost difference'=dollar_format(accuracy=1)(yes-no))
    return(tset)
}
propensity_score_analysis(pset)

mediation_analysis <- function(dset,
                               response='Cost',
                               mediator='Pain',
                               cause='Arthritis',
                               confounders=c('Sex',
                                             'Diabetes',
                                             'Age')){
    dset <- dset %>% filter(Anycost)
    dset <- dset %>% mutate(Pain=as.integer(Pain))
    predictor1 <- paste(c(cause, confounders), collapse='+')
    predictor2 <- paste(c(cause, confounders), collapse='+')
    predictor3 <- paste(c(cause, confounders, mediator), collapse='+')
    model1 <- paste(response, predictor1, sep='~')
    model2 <- paste(mediator, predictor2, sep='~')
    model3 <- paste(response, predictor3, sep='~')
    fit1 <- lm(as.formula(model1), data=dset)
    fit2 <- lm(as.formula(model2), data=dset)
    fit3 <- lm(as.formula(model3), data=dset)
    formatted1 <- format_model(fit1, ' (1)')
    formatted2 <- format_model(fit2, ' (2)')
    formatted3 <- format_model(fit3, ' (3)')
    fset <- full_join(formatted1, formatted2, by='Variable')
    fset <- full_join(fset, formatted3, by='Variable')
    total_effect <- coef(fit1)['Arthritisyes']
    effect_on_pain <- coef(fit2)['Arthritisyes']
    effect_of_pain <- coef(fit3)['Pain']
    indirect_effect <- effect_on_pain*effect_of_pain
    direct_effect <- total_effect-indirect_effect
    cat('Total effect: $', total_effect, '\n', sep='')
    cat('Indirect effect: $', indirect_effect, '\n', sep='')
    cat('Direct effect: $', direct_effect, '\n', sep='')
    cat('Indirect/Total: ', indirect_effect/total_effect, '\n', sep='')
    fit3_interaction <- update(fit3, .~.+Arthritis:Pain)
    cat('Model (3) after adding cause-mediator interaction:\n')
    return(coef(summary(fit3_interaction)))
}
mediation_analysis(dset)

