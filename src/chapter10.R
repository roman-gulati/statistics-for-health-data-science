##################################################
# Chapter 10 examples
##################################################
library(tidyverse)
library(scales)
library(grid)
library(viridis)
library(glmnet)
library(rpart)
library(rpart.plot)
library(randomForest)

set.seed(1233)

source('grab_meps.R')
source('grab_nhanes.R')

##################################################
# Load or download/extract/save longitudinal MEPS
##################################################
grab_longitudinal_meps <- function(panel){
    code <- c(23, 35, 48, 58, 65, 71, 80, 86, 98, 106, 114, 122, 130, 139, 148, 156, 164, 172, 183, 193, 202)
    code <- paste0('h', code)
    stopifnot(between(panel, 1, 21))
    year_code <- code[panel]
    archive_filename <- paste0(year_code, 'ssp.zip')
    remote_sas_archive <- paste0('https://meps.ahrq.gov/mepsweb/data_files/pufs/',
                                 archive_filename)
    local_sas_archive <- paste0(year_code, 'ssp.zip')
    local_rdata <- paste(year_code, 'Rdata', sep='.')
    if(file.exists(local_rdata)){
        load(file=local_rdata)
    } else {
        download.file(remote_sas_archive, local_sas_archive, mode='wb')
        dset <- foreign::read.xport(unzip(local_sas_archive))
        local_sas_file <- sub('ssp.zip$', '.ssp', local_sas_archive)
        file.copy(local_sas_file, local_rdata, overwrite=TRUE)
        unlink(local_sas_file)
        save(dset, file=local_rdata)
    }
    dset <- dset %>% filter(INSCOPY1 == 1, INSCOPY2 == 1)
    return(dset)
}

##################################################
# Curate MEPS data from a longitudinal panel with
# predictor variables from a specified survey year
# or round
##################################################
grab_and_curate_meps <- function(panel=21, survey_round=3, survey_year=1){
    stopifnot(between(survey_round, 1, 5))
    dset <- grab_longitudinal_meps(panel)
    variables <- c(Sex='SEX',
                   Race='RACEV1X',
                   Education='EDUCYR',
                   Age=paste0('AGEY', survey_year, 'X'),
                   Married=paste0('MARRYY', survey_year, 'X'),
                   Region=paste0('REGIONY', survey_year),
                   Poverty=paste0('POVCATY', survey_year),
                   Insurance=paste0('INSCOVY', survey_year),
                   BP.diagnosis=paste0('HIBPDXY', survey_year),
                   MI.diagnosis=paste0('MIDXY', survey_year),
                   CHD.diagnosis=paste0('CHDDXY', survey_year),
                   OHD.diagnosis=paste0('OHRTDXY', survey_year),
                   Angina.diagnosis=paste0('CHDDXY', survey_year),
                   Asthma.diagnosis=paste0('ASTHDXY', survey_year),
                   Arthritis.diagnosis=paste0('ARTHDXY', survey_year),
                   Cancer.diagnosis=paste0('CANCERY', survey_year),
                   Cholesterol.diagnosis=paste0('CHOLDXY', survey_year),
                   Diabetes.diagnosis=paste0('DIABDXY', survey_year),
                   Emphysema.diagnosis=paste0('EMPHDXY', survey_year),
                   Stroke.diagnosis=paste0('STRKDXY', survey_year),
                   Perceived.health=paste0('RTHLTH', survey_round),
                   Mental.health=paste0('MNHLTH', survey_round),
                   Inpatient.year1='IPTEXPY1',
                   Inpatient.year2='IPTEXPY2',
                   Outpatient.year1='OPTEXPY1',
                   Outpatient.year2='OPTEXPY2',
                   Total.year1='TOTEXPY1',
                   Total.year2='TOTEXPY2')
    # some variables are only available for certain rounds
    if(survey_round %in% c(3, 5)){
        variables <- c(variables, Physical.exercise=paste0('PHYEXE', survey_round))
        if(survey_round == 3)
            variables <- c(variables, BMI=paste0('BMINDX', survey_round))
    } else if (survey_round %in% c(2, 4)){
        variables <- c(variables, Current.smoker=paste0('ADSMOK', survey_round))
    }
    dset <- dset %>% select(all_of(variables))
    dset <- dset %>% mutate(Sex=factor(Sex,
                                       levels=c(1, 2),
                                       labels=c('Male', 'Female')),
                            Race=factor(Race,
                                        levels=c(1, 2, 3, 4, 6),
                                        labels=c('White',
                                                 'Black',
                                                 'Amer Indian/Alasaka Native',
                                                 'Asian/Native Hawaiian/Pacific Islander',
                                                 'Mixed')),
                            Education=ifelse(between(Education, 0, 17), Education, NA),
                            Married=factor(Married,
                                           levels=c(1, 2, 3, 4, 5, 6),
                                           labels=c('Married',
                                                    'Widowed',
                                                    'Divorced',
                                                    'Separated',
                                                    'Never married',
                                                    'Inapplicable')),
                            Region=factor(Region,
                                          levels=c(1, 2, 3, 4),
                                          labels=c('Northeast',
                                                   'South',
                                                   'Midwest',
                                                   'West')),
                            Poverty=factor(Poverty,
                                           levels=c(1, 2, 3, 4, 5),
                                           labels=c('Poor',
                                                    'Near poor',
                                                    'Low income',
                                                    'Middle income',
                                                    'High income')),
                            Insurance=factor(Insurance,
                                             levels=c(1, 2, 3),
                                             labels=c('Any private',
                                                      'Public only',
                                                      'Uninsured')),
                            BP.diagnosis=factor(BP.diagnosis,
                                                levels=c(2, 1),
                                                labels=c('No', 'Yes')),
                            MI.diagnosis=factor(MI.diagnosis,
                                                levels=c(2, 1),
                                                labels=c('No', 'Yes')),
                            CHD.diagnosis=factor(CHD.diagnosis,
                                                 levels=c(2, 1),
                                                 labels=c('No', 'Yes')),
                            OHD.diagnosis=factor(OHD.diagnosis,
                                                 levels=c(2, 1),
                                                 labels=c('No', 'Yes')),
                            Angina.diagnosis=factor(Angina.diagnosis,
                                                    levels=c(2, 1),
                                                    labels=c('No', 'Yes')),
                            Asthma.diagnosis=factor(Asthma.diagnosis,
                                                    levels=c(2, 1),
                                                    labels=c('No', 'Yes')),
                            Arthritis.diagnosis=factor(Arthritis.diagnosis,
                                                       levels=c(2, 1),
                                                       labels=c('No', 'Yes')),
                            Cancer.diagnosis=factor(Cancer.diagnosis,
                                                    levels=c(2, 1),
                                                    labels=c('No', 'Yes')),
                            Cholesterol.diagnosis=factor(Cholesterol.diagnosis,
                                                         levels=c(2, 1),
                                                         labels=c('No', 'Yes')),
                            Diabetes.diagnosis=factor(Diabetes.diagnosis,
                                                      levels=c(2, 1),
                                                      labels=c('No', 'Yes')),
                            Emphysema.diagnosis=factor(Emphysema.diagnosis,
                                                       levels=c(2, 1),
                                                       labels=c('No', 'Yes')),
                            Stroke.diagnosis=factor(Stroke.diagnosis,
                                                    levels=c(2, 1),
                                                    labels=c('No', 'Yes')),
                            Perceived.health=factor(Perceived.health,
                                                    levels=seq(5),
                                                    labels=c('Excellent',
                                                             'Very good',
                                                             'Good',
                                                             'Fair',
                                                             'Poor')),
                            Mental.health=factor(Mental.health,
                                                 levels=seq(5),
                                                 labels=c('Excellent',
                                                          'Very good',
                                                          'Good',
                                                          'Fair',
                                                          'Poor')),
                            Anyhosp2016=Inpatient.year1 > 0,
                            Anyhosp2017=Inpatient.year2 > 0,
                            Edu9plus=Education >= 9,
                            Edu12plus=Education >= 12,
                            Edu13plus=Education >= 13,
                            Edu16plus=Education >= 16)
    if('Physical.exercise' %in% names(dset))
        dset <- dset %>% mutate(Physical.exercise=factor(Physical.exercise,
                                                         levels=c(2, 1),
                                                         labels=c('No', 'Yes')))
    if('BMI' %in% names(dset))
        dset <- dset %>% mutate(BMI=ifelse(between(BMI, 12.5, 83.2), BMI, NA))
    if('Current.smoker' %in% names(dset))
        dset <- dset %>% mutate(Current.smoker=factor(Current.smoker,
                                                      levels=c(2, 1),
                                                      labels=c('No', 'Yes')))
    cat('Remove',
        dset %>% filter(Age <= 17) %>% nrow(),
        'age 17 or younger\n')
    dset <- dset %>% filter(Age > 17)
    return(dset)
}

##################################################
# Illustrate data science concepts using NHANES 2015
##################################################
iset <- grab_and_curate_nhanes(2015, minage=0, maxage=80)

##################################################
# Define age groups of various sizes
##################################################
append_agegroups <- function(dset, ngroups=c(2, 5, 10, 20, 40, 79)){
    for(n in ngroups)
        iset <- iset %>% mutate(!!sym(paste0('AGEGROUPS', n)):=cut_interval(AGE, n=n))
    return(iset)
}
iset <- append_agegroups(iset)

##################################################
# Slice out training and testing data
##################################################
iset.train <- iset %>% slice(1:1000)
iset.test <- iset %>% slice(1001:2000)

##################################################
# Fit linear models for BMI given age group and
# collect predictions for training and testing data
##################################################
predictors <- grep('^AGEGROUPS', names(iset), value=TRUE)
fits.train <- lapply(predictors, function(x) lm(paste('BMI', x, sep='~'), data=iset.train))
fits.test <- lapply(predictors, function(x) lm(paste('BMI', x, sep='~'), data=iset.test))

predict_bmi <- function(fits, dset, predictors){
    plist <- lapply(fits, predict, newdata=dset)
    pset <- do.call(cbind, plist)
    pset <- with(dset, data.frame(BMI, AGE, pset))
    names(pset) <- c('BMI', 'AGE', predictors)
    pset <- pset %>% gather(Predictor, Prediction, -BMI, -AGE)
    pset <- pset %>% mutate(Predictor=as.factor(sub('^AGEGROUPS', '', Predictor)))
    return(pset)
}
pset.train <- predict_bmi(fits.train, iset.train, predictors)
pset.test <- predict_bmi(fits.test, iset.test, predictors)

##################################################
# Visualize selected fitted models
##################################################
tidy_predictions <- function(dset, correction=0.5){
    pset <- dset %>% select(Predictor, AGE, Prediction)
    pset <- pset %>% unique()
    pset <- pset %>% mutate(AGE=AGE-correction)
    pset <- pset %>% group_by(Predictor)
    pset <- pset %>% ungroup()
    xset <- pset %>% filter(AGE == max(AGE))
    xset <- xset %>% mutate(AGE=AGE+correction*2)
    pset <- bind_rows(pset, xset)
    pset <- pset %>% arrange(Predictor, AGE, Prediction)
    return(pset)
}

bmi_by_agegroup_plot <- function(train, test=NA){
    if(!is.null(dim(test))){
        rset <- test
        train.alpha <- 0.5
    } else {
        rset <- train
        train.alpha <- 1
    }
    gg <- ggplot(rset)
    gg <- gg+geom_point(aes(x=AGE, y=BMI), color='gray', size=1)
    if(!is.null(dim(test)))
        gg <- gg+geom_step(data=tidy_predictions(test),
                           aes(x=AGE,
                               y=Prediction,
                               colour=Predictor),
                           size=1)
    gg <- gg+geom_step(data=tidy_predictions(train),
                       aes(x=AGE,
                           y=Prediction,
                           colour=Predictor),
                       alpha=train.alpha,
                       size=1)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(0, 81),
                                breaks=seq(0, 80, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_colour_manual(values=c('2'='#5DC863FF',
                                          '5'='#2FB47CFF',
                                          '10'='#1E9C89FF',
                                          '20'='#25848EFF',
                                          '40'='#2F6C8EFF',
                                          '79'='#3B528BFF',
                                          'CART'='purple'))
    print(gg)
}
bmi_by_agegroup_plot(pset.train %>% filter(Predictor %in% c(2, 79)))
bmi_by_agegroup_plot(pset.train %>% filter(Predictor == 2),
                     pset.test %>% filter(Predictor == 2))
bmi_by_agegroup_plot(pset.train %>% filter(Predictor == 79),
                     pset.test %>% filter(Predictor == 79))

##################################################
# Calculate mean squared errors
##################################################
calculate_mse <- function(fits, dset, predictors, dataset){
    dset <- predict_bmi(fits, dset, predictors)
    dset <- dset %>% mutate(Error=BMI-Prediction)
    dset <- dset %>% group_by(Predictor)
    mset <- dset %>% summarize(Dataset=dataset,
                               MSE=mean(Error^2)/var(BMI))
    return(mset)
}

mse_by_agegroup_size <- function(fits, train, test, predictors){
    mse.train <- calculate_mse(fits, train, predictors, 'Training')
    mse.test <- calculate_mse(fits, test, predictors, 'Testing')
    mset <- bind_rows(mse.train, mse.test)
    mset <- mset %>% mutate(Predictor=factor(as.integer(as.character(Predictor))))
    gg <- ggplot(mset)
    gg <- gg+geom_point(aes(x=Predictor,
                            y=MSE,
                            colour=Predictor,
                            alpha=Dataset),
                        size=4)
    gg <- gg+scale_x_discrete(name='Number of age groups')
    gg <- gg+scale_y_continuous(name='Standardized mean squared error',
                                limits=c(0.5, 0.9),
                                breaks=seq(0, 1, by=0.1),
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(discrete=TRUE, begin=0.75, end=0.25)
    gg <- gg+scale_alpha_manual(values=c('Training'=0.4, 'Testing'=1))
    print(gg)
}
mse_by_agegroup_size(fits.train, iset.train, iset.test, predictors)

##################################################
# Visualize optimal fitted model
##################################################
bmi_by_agegroup_plot(pset.test %>% filter(Predictor == 10))

##################################################
# K-fold cross-validation
##################################################
cv_mse <- function(cv_vars, dset, k=10){
    test.index <- sample(nrow(dset), nrow(dset)/k)
    test.set <- dset %>% slice(test.index)
    train.set <- dset %>% slice(-test.index)
    fit <- lm(paste('BMI', cv_vars$Predictor, sep='~'), data=train.set)
    test.set <- test.set %>% mutate(Prediction=predict(fit, newdata=test.set))
    test.set <- test.set %>% mutate(Error=BMI-Prediction)
    test.set <- test.set %>% summarize(MSE=mean(Error^2)/var(BMI))
    return(test.set)
}

cv_mse_by_agegroup_size <- function(dset){
    cvset <- expand.grid(Predictor=predictors, Replicate=seq(10))
    cvset <- cvset %>% group_by(Predictor, Replicate)
    cvset <- cvset %>% do(cv_mse(., dset))
    cvset <- cvset %>% ungroup(Predictor)
    cvset <- cvset %>% mutate(Predictor=as.integer(sub('^AGEGROUPS', '', Predictor)))
    cvset <- cvset %>% mutate(Predictor=factor(Predictor))
    meset <- cvset %>% group_by(Predictor)
    meset <- meset %>% summarize(Mean=mean(MSE))
    gg <- ggplot(cvset)
    gg <- gg+geom_point(aes(x=Predictor,
                            y=MSE,
                            colour=Predictor),
                        alpha=0.4,
                        size=2)
    gg <- gg+geom_point(data=meset,
                        aes(x=Predictor,
                            y=Mean,
                            colour=Predictor),
                        size=5)
    gg <- gg+scale_x_discrete(name='Number of age groups')
    gg <- gg+scale_y_continuous(name='Standardized mean squared error',
                                limits=c(0, 1),
                                breaks=seq(0, 1, by=0.2),
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(discrete=TRUE, begin=0.75, end=0.25)
    print(gg)
}
cv_mse_by_agegroup_size(iset.train)

##################################################
# Make design and prediction matrices to show
# LASSO and ridge regression
##################################################
obs.ages <- iset.train$AGE
uni.ages <- sort(unique(obs.ages))
gmat <- outer(obs.ages, uni.ages, '<')
pmat <- upper.tri(matrix(0, length(uni.ages), length(uni.ages)))

##################################################
# LASSO
##################################################
cv_lasso <- function(dset, gmat, pmat, uni.ages){
    cv.lasso.fit <- cv.glmnet(x=gmat, y=dset$BMI, alpha=1, family='gaussian')
    pred.lasso <- tibble(AGE=uni.ages,
                         Predictor='LASSO',
                         Prediction=as.vector(predict(cv.lasso.fit, newx=pmat, s='lambda.min')))
    return(pred.lasso)
}
pred.lasso <- cv_lasso(iset.train, gmat, pmat, uni.ages)

##################################################
# Ridge regression
##################################################
cv_ridge <- function(dset, gmat, pmat, uni.ages){
    cv.ridge.fit <- cv.glmnet(x=gmat, y=iset.train$BMI, alpha=0, family='gaussian')
    pred.ridge <- tibble(AGE=uni.ages,
                         Predictor='Ridge',
                         Prediction=as.vector(predict(cv.ridge.fit, newx=pmat, s='lambda.min')))
    return(pred.ridge)
}
pred.ridge <- cv_ridge(dset, gmat, pmat, uni.ages)

##################################################
# Visualize standard and regularized predictions
##################################################
bmi_by_age_standard_and_regularized_plot <- function(pred.standard,
                                                     pred.lasso,
                                                     pred.ridge){
    sset <- tidy_predictions(pred.standard)
    sset <- sset %>% mutate(Predictor=paste(Predictor, 'age groups'))
    theme_update(legend.position=c(0.15, 0.85))
    gg <- ggplot(pred.standard)
    gg <- gg+geom_point(aes(x=AGE, y=BMI), color='gray', size=1)
    gg <- gg+geom_step(data=sset,
                       aes(x=AGE,
                           y=Prediction,
                           colour=Predictor),
                       alpha=0.4,
                       size=1)
    gg <- gg+geom_step(data=pred.lasso,
                       aes(x=AGE,
                           y=Prediction,
                           colour=Predictor),
                       size=1)
    gg <- gg+geom_step(data=pred.ridge,
                       aes(x=AGE,
                           y=Prediction,
                           colour=Predictor),
                       size=1)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(0, 81),
                                breaks=seq(0, 80, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_colour_manual(name='Model',
                                 values=c('2 age groups'='#5DC863FF',
                                          '5 age groups'='#2FB47CFF',
                                          '10 age groups'='#1E9C89FF',
                                          '20 age groups'='#25848EFF',
                                          '40 age groups'='#2F6C8EFF',
                                          '79 age groups'='#3B528BFF',
                                          'LASSO'='#7AD151FF',
                                          'Ridge'='#483480FF'))
    print(gg)
}
bmi_by_age_standard_and_regularized_plot(pset.train %>% filter(Predictor == 79),
                                         pred.lasso,
                                         pred.ridge)

##################################################
# Read longitudinal MEPS data
##################################################
mset <- grab_and_curate_meps(panel=21, survey_round=2)

##################################################
# Restrict to candidate predictor variables
##################################################
predictors <- c('Age', 'Sex', 'Race', 'Married', 'Poverty',
                'Edu9plus', 'Edu12plus', 'Edu13plus', 'Edu16plus',
                'Insurance',
                'Angina.diagnosis',
                'Asthma.diagnosis',
                'Arthritis.diagnosis',
                'Cancer.diagnosis',
                'CHD.diagnosis',
                'Diabetes.diagnosis',
                'Emphysema.diagnosis',
                'MI.diagnosis',
                'BP.diagnosis',
                'Cholesterol.diagnosis',
                'OHD.diagnosis',
                'Stroke.diagnosis',
                'Perceived.health',
                'Mental.health',
                'Anyhosp2016')
cat('Keep',
    mset %>% drop_na(any_of(predictors)) %>% nrow(),
    'with non-missing data\n')
mset <- mset %>% drop_na(any_of(predictors))

##################################################
# Remove intercept and empty variable then
# fit standard, LASSO, and ridge regression
##################################################
fit_standard_regularized <- function(dset, predictors){
    mmat <- model.matrix(~., dset %>% select(all_of(predictors)))
    takeout <- which(colnames(mmat) %in% c('MarriedInapplicable', '(Intercept)'))
    mmat <- mmat[, -takeout]
    fit.standard <- glmnet(x=mmat, y=dset$Anyhosp2017, lambda=0, family='binomial')
    fit.lasso <- cv.glmnet(x=mmat, y=dset$Anyhosp2017, alpha=1, family='binomial')
    fit.ridge <- cv.glmnet(x=mmat, y=dset$Anyhosp2017, alpha=0, family='binomial')
    fset <- tibble(Predictor=rownames(coef(fit.standard)),
                   Standard=as.vector(coef(fit.standard)),
                   LASSO=as.vector(coef(fit.lasso)),
                   Ridge=as.vector(coef(fit.ridge)))
    return(fset)
}
fset <- fit_standard_regularized(mset, predictors)

##################################################
# Show standard and regularized regression results
##################################################
standard_and_regularized_barplot <- function(dset){
    pset <- dset %>% gather(Model, Coefficient, -Predictor)
    pset <- pset %>% filter(Predictor != '(Intercept)')
    pset <- pset %>% mutate(Predictor=factor(Predictor, levels=unique(Predictor)),
                            Model=factor(Model, levels=c('Standard', 'LASSO', 'Ridge')))
    theme_update(legend.position=c(0.8, 0.2),
                 axis.text.x=element_blank())
    gg <- ggplot(pset)
    gg <- gg+geom_hline(aes(yintercept=0))
    gg <- gg+geom_bar(aes(x=Predictor,
                          y=Coefficient,
                          fill=Model),
                      stat='identity',
                      position='dodge')
    gg <- gg+scale_x_discrete(name='Predictor')
    gg <- gg+scale_y_continuous(name='Coefficient',
                                limits=c(-1, 1),
                                breaks=seq(-1, 1, by=0.2),
                                expand=c(0, 0))
    gg <- gg+scale_fill_viridis(discrete=TRUE, option='cividis', begin=0.9, end=0.1)
    print(gg)
}
standard_and_regularized_barplot(fset)

standard_and_regularized_coefficients <- function(fset, mset){
    fset <- fset %>% mutate(Predictor=sub('(Intercept)', 'Intercept', Predictor, fixed=TRUE),
                            Predictor=gsub(paste0('(', paste(names(mset), collapse='|'), ')'), '\\1 = ', Predictor),
                            Predictor=sub('Age =', 'Age', Predictor),
                            Predictor=sub('[.]', ' ', Predictor),
                            Predictor=sub('([0-9]*)plus = TRUE', ' = \\1+', Predictor),
                            Predictor=sub('Edu', 'Education', Predictor),
                            Predictor=sub('CHD', 'Coronary heart disease', Predictor),
                            Predictor=sub('MI', 'Heart attack', Predictor),
                            Predictor=sub('BP', 'High blood pressure', Predictor),
                            Predictor=sub('Cholesterol', 'High cholesterol', Predictor),
                            Predictor=sub('OHD', 'Other heart disease', Predictor),
                            Predictor=sub('Anyhosp2016 = TRUE', 'Hospitalized in 2016 = Yes', Predictor))
    return(fset)
}
standard_and_regularized_coefficients(fset, mset)

##################################################
# Use CART to predict BMI by age
##################################################
tree_mse <- function(dset, splits=NA){
    if(all(is.na(splits)))
        predictors <- 1
    else {
        predictors <- paste('AGE', splits, sep='<')
        predictors <- paste0('I(', predictors, ')')
        predictors <- paste(predictors, collapse='+')
    }
    equation <- paste('BMI', predictors, sep='~')
    fit <- lm(equation, data=dset)
    return(summary(fit)$sigma^2)
}

bmi_tree_analysis <- function(dset){
    fit.tree <- rpart(BMI~AGE, data=dset, method='anova')
    #printcp(fit.tree)
    rpart.plot(fit.tree, tweak=1.05)
    mse0 <- tree_mse(dset, splits=NA)
    mse1 <- tree_mse(dset, splits=13)
    mse2 <- tree_mse(dset, splits=c(13, 26))
    mse3 <- tree_mse(dset, splits=c(13, 26, 71))
    mses <- c(mse0, mse1, mse2, mse3)
    ratios <- 1-lead(mses)/mses
    cat('Reductions in MSEs:', ratios, '\n')
    pset <- dset %>% mutate(Predictor='CART', Prediction=predict(fit.tree))
    bmi_by_agegroup_plot(pset)
}
bmi_tree_analysis(iset.train)

##################################################
# Use CART to predict total expenditures
##################################################
tree_error_plot <- function(dset){
    gg <- ggplot(dset)
    gg <- gg+geom_hline(data=dset %>% filter(xerror+xstd == min(xerror+xstd)),
                        aes(yintercept=xerror+xstd),
                        colour='orange')
    gg <- gg+geom_pointrange(aes(x=nsplit+1,
                                 y=xerror,
                                 ymin=xerror-xstd,
                                 ymax=xerror+xstd))
    gg <- gg+scale_x_continuous(name='Size of tree',
                                limits=c(1, 55),
                                breaks=seq(1, 55, by=2),
                                minor_breaks=seq(2, 54, by=2))
    gg <- gg+scale_y_continuous(name='Relative error',
                                limits=c(0.7, 1.1),
                                breaks=seq(0.7, 1.1, by=0.1),
                                expand=c(0, 0))
    print(gg)
}

expenditure_tree_analysis <- function(dset, predictors){
    dset <- dset %>% select(all_of(predictors), Total.year1, Total.year2)
    dset <- dset %>% rename('Total expenditures in 2016'='Total.year1',
                            'Perceived health'='Perceived.health')
    cost.tree <- rpart(Total.year2~., data=dset, method='anova')
    #printcp(cost.tree)
    #rpart.plot(cost.tree)
    # explore greater depth
    cost.tree2 <- rpart(Total.year2~.,
                        data=dset,
                        method='anova',
                        control=rpart.control(cp=0.001, maxdepth=20))
    #plotcp(cost.tree2)
    tree_error_plot(as.data.frame(cost.tree2$cptable))
    best.cp <- cost.tree2$cptable[which.min(cost.tree2$cptable[, 'xerror']), 'CP']
    # constrain depth
    cost.tree3 <- rpart(Total.year2~.,
                        data=dset,
                        method='anova',
                        control=rpart.control(cp=best.cp))
    rpart.plot(cost.tree3)
    rpart.plot(cost.tree3, tweak=1.05)
    # calculate summaries within prediction strata
    pset <- dset %>% mutate(Prediction=predict(cost.tree3))
    pset <- pset %>% select('Total expenditures in 2016',
                            'Perceived health',
                            Age,
                            Prediction)
    pset <- pset %>% group_by(Prediction)
    pset <- pset %>% summarize('Total expenditures in 2016'=paste(dollar_format()(min(`Total expenditures in 2016`)),
                                                                  dollar_format()(max(`Total expenditures in 2016`)), sep='-'),
                               'Perceived health'=paste(sort(unique(`Perceived health`)), collapse=', '),
                               Age=paste(min(Age), max(Age), sep='-'))
    pset <- pset %>% mutate('Perceived health'=ifelse(grepl('Very good, Good', `Perceived health`), 'Any', `Perceived health`))
    pset <- pset %>% select('Total expenditures in 2016',
                            'Perceived health',
                            Age,
                            Prediction)
    # unfortunately, summaries miss gaps between groups
    # manually construct table
    print(cost.tree3)
    pset <- tibble('Costs in 2016'=c('$0--$2,127',
                                     '$2,128--$11,004',
                                     '$11,005--$40,221',
                                     '$11,005--$40,221',
                                     '>$40,221',
                                     '$40,222--$72,906',
                                     '>$72,906',
                                     '>$72,906'),
                   'Perceived health'=c('Any',
                                        'Any',
                                        'Excellent, Very good',
                                        'Good, Fair, Poor',
                                        'Any',
                                        'Any',
                                        'Any',
                                        'Any'),
                   'Insurance'=c('Any',
                                 'Any',
                                 'Any',
                                 'Any',
                                 'Any private, Uninsured',
                                 'Public only',
                                 'Public only',
                                 'Public only'),
                   Age=c('Any',
                         'Any',
                         'Any',
                         'Any',
                         'Any',
                         'Any',
                         '18--53',
                         '54--85'),
                   'Predicted in 2017'=c('$2,064',
                                         '$7,895',
                                         '$12,186',
                                         '$22,149',
                                         '$29,884',
                                         '$39,288',
                                         '$35,261',
                                         '$95,723'))
    return(pset)
}
expenditure_tree_analysis(mset, predictors)

##################################################
# Use CART to predict hospitalization. Illustrate
# classification tree with splits even though this
# isn't optimal
##################################################
hospitalization_tree_analysis <- function(dset, predictors){
    dset <- dset %>% select(all_of(predictors), Total.year1, Anyhosp2017)
    dset <- dset %>% mutate(Edu12plus=factor(Edu12plus,
                                             levels=c('FALSE', 'TRUE'),
                                             labels=c('No', 'Yes')))
    dset <- dset %>% rename('Total expenditures in 2016'='Total.year1',
                            'Perceived health'='Perceived.health',
                            'Angina diagnosis'='Angina.diagnosis',
                            'Asthma diagnosis'='Asthma.diagnosis',
                            'Heart attack diagnosis'='MI.diagnosis',
                            'Education 12+'='Edu12plus')
    hosp.tree <- rpart(Anyhosp2017~.,
                       data=dset,
                       method='class',
                       control=rpart.control(minsplit=20, cp=0))
    printcp(hosp.tree)
    #plotcp(hosp.tree)
    #rpart.plot(hosp.tree)
    # explore greater depth
    hosp.tree2 <- rpart(Anyhosp2017~.,
                        data=dset,
                        method='class',
                        control=rpart.control(minsplit=20, cp=0.0032))
    plotcp(hosp.tree2, col='blue', lty='dashed')
    rpart.plot(hosp.tree2)
    #tree_error_plot(as.data.frame(hosp.tree2$cptable))
    best.cp <- hosp.tree2$cptable[which.min(hosp.tree2$cptable[, 'xerror']), 'CP']
    # constrain depth
    hosp.tree3 <- rpart(Anyhosp2017~.,
                        data=dset,
                        method='class',
                        control=rpart.control(cp=best.cp))
    #rpart.plot(hosp.tree3)
    pset <- dset %>% mutate(Prediction=predict(hosp.tree2, type='class'))
    ptab <- pset %>% with(prop.table(table(Anyhosp2017, Prediction)))
    ptab <- ptab %>% as.data.frame()
    ptab <- ptab %>% filter(Anyhosp2017 != Prediction)
    cat('Misclassification error:', ptab %>% with(sum(Freq)), '\n')
    rpart.plot(hosp.tree2, tweak=1.05)
}
hospitalization_tree_analysis(mset, predictors)

##################################################
# Examine random forests
##################################################
forest_importance_plot <- function(rfobject){
    iset <- importance(rfobject)
    iset <- data.frame(Predictor=rownames(iset),
                       iset,
                       row.names=NULL,
                       check.names=FALSE)
    iset <- iset %>% rename('PCMSE'='%IncMSE')
    iset <- iset %>% arrange(PCMSE)
    iset <- iset %>% mutate(Predictor=sub('[.]', ' ', Predictor),
                            Predictor=sub('Total year1', 'Total expenditures in 2016', Predictor),
                            Predictor=sub('Edu(.*)plus', 'Education = \\1+', Predictor),
                            Predictor=sub('CHD', 'Coronary heart disease', Predictor),
                            Predictor=sub('MI', 'Heart attack', Predictor),
                            Predictor=sub('BP', 'High blood pressure', Predictor),
                            Predictor=sub('Cholesterol', 'High cholesterol', Predictor),
                            Predictor=sub('OHD', 'Other heart disease', Predictor),
                            Predictor=sub('Anyhosp2016', 'Hospitalized in 2016', Predictor),
                            Predictor=factor(Predictor, levels=Predictor))
    theme_update(legend.position='none',
                 axis.text.y=element_text(size=10))
    gg <- ggplot(iset)
    gg <- gg+geom_point(aes(x=PCMSE/100,
                            y=Predictor,
                            colour=Predictor),
                        size=4)
    gg <- gg+scale_x_continuous(name='Increase in mean squared error',
                                label=percent_format(accuracy=1))
    gg <- gg+scale_colour_viridis(discrete=TRUE)
    gg <- gg+scale_y_discrete(name='')
    print(gg)
}

expenditure_forest_analysis <- function(dset, predictors, p=0.7){
    dset <- dset %>% select(all_of(predictors), Total.year1, Total.year2)
    total.rows <- nrow(dset)
    train.rows <- sample(total.rows, total.rows*p)
    test.rows <- setdiff(seq(total.rows), train.rows)
    dset.train <- dset %>% slice(train.rows)
    dset.test <- dset %>% slice(test.rows)
    cost.forest <- randomForest(Total.year2~.,
                                data=dset.train,
                                ntree=500,
                                importance=TRUE)
    dset.test <- dset.test %>% mutate(Forest=predict(cost.forest, newdata=dset.test))
    forest.mse <- dset.test %>% summarize(MSE=mean((Total.year2-Forest)^2))
    forest.mse <- forest.mse %>% unlist()
    cat('Forest MSE:', forest.mse, '\n')
    forest_importance_plot(cost.forest)
    # grow large regression tree
    cost.tree <- rpart(Total.year2~.,
                       data=dset.train,
                       method='anova',
                       control=rpart.control(cp=0.001))
    #plotcp(cost.tree)
    #printcp(cost.tree)
    #tree_error_plot(as.data.frame(cost.tree$cptable))
    best.cp <- cost.tree$cptable[which.min(cost.tree$cptable[, 'xerror']), 'CP']
    # best tree
    cost.tree2 <- rpart(Total.year2~.,
                        data=dset.train,
                        method='anova',
                        control=rpart.control(cp=best.cp))
    #rpart.plot(cost.tree2)
    dset.test <- dset.test %>% mutate(Tree=predict(cost.tree2, newdata=dset.test))
    tree.mse <- dset.test %>% summarize(MSE=mean((Total.year2-Tree)^2))
    tree.mse <- tree.mse %>% unlist()
    cat('Tree MSE:', tree.mse, '\n')
    cat('Improvement in forest over tree:', (tree.mse-forest.mse)/tree.mse, '\n')
}
expenditure_forest_analysis(mset, predictors)

hospitalization_forest_analysis <- function(dset, predictors, p=0.7){
    dset <- dset %>% select(all_of(predictors), Total.year1, Anyhosp2017)
    dset <- dset %>% mutate(Anyhosp2017=factor(Anyhosp2017))
    total.rows <- nrow(dset)
    train.rows <- sample(total.rows, total.rows*p)
    test.rows <- setdiff(seq(total.rows), train.rows)
    dset.train <- dset %>% slice(train.rows)
    dset.test <- dset %>% slice(test.rows)
    hosp.forest <- randomForest(Anyhosp2017~.,
                                data=dset.train,
                                ntree=500,
                                importance=TRUE)
    dset.test <- dset.test %>% mutate(Forest=predict(hosp.forest, newdata=dset.test))
    ftab <- dset.test %>% with(prop.table(table(Anyhosp2017, Forest)))
    ftab <- ftab %>% as.data.frame()
    ftab <- ftab %>% filter(Anyhosp2017 != Forest)
    cat('Misclassification error for forest:', ftab %>% with(sum(Freq)), '\n')
    # grow large classification tree
    hosp.tree <- rpart(Anyhosp2017~.,
                       data=dset.train,
                       method='class',
                       control=rpart.control(cp=0.001))
    #printcp(hosp.tree)
    #tree_error_plot(as.data.frame(hosp.tree$cptable))
    best.cp <- hosp.tree$cptable[which.min(hosp.tree$cptable[, 'xerror']), 'CP']
    # best tree
    hosp.tree2 <- rpart(Anyhosp2017~.,
                        data=dset.train,
                        method='class',
                        control=rpart.control(cp=best.cp))
    dset.test <- dset.test %>% mutate(Tree=predict(hosp.tree2, newdata=dset.test, typ='class'))
    ttab <- dset.test %>% with(prop.table(table(Anyhosp2017, Tree)))
    ttab <- ttab %>% as.data.frame()
    ttab <- ttab %>% filter(Anyhosp2017 != Tree)
    cat('Misclassification error for tree:', ttab %>% with(sum(Freq)), '\n')
}
hospitalization_forest_analysis(mset, predictors)

