##################################################
# Chapter 3 examples
##################################################
library(tidyverse)
library(scales)
library(grid)
library(viridis)
library(quantreg)

set.seed(1233)

source('grab_nhanes.R')

##################################################
# Visualize variability around a regression line
##################################################
regression_variability_plot <- function(dset, type){
    dset <- dset %>% mutate(p=predict(lm(y~x)))
    panel_title <- switch(type, y='Variability around E(Y)', regression='Variability around E(Y|X)')
    gg <- ggplot(dset)
    gg <- gg+ggtitle(panel_title)
    if(type == 'y'){
        gg <- gg+geom_segment(aes(x=x, y=y, xend=x, yend=mean(y)),
                              colour='gray',
                              size=0.75)
        gg <- gg+geom_hline(aes(yintercept=mean(y)),
                            size=0.75)
    } else {
        gg <- gg+geom_segment(aes(x=x, y=y, xend=x, yend=p),
                              colour='gray',
                              size=0.75)
        gg <- gg+geom_line(aes(x=x, y=p),
                           size=0.75)
    }
    gg <- gg+geom_point(aes(x=x, y=y), size=2)
    gg <- gg+labs(x='\nx', y='y\n')
    print(gg)
}
dset <- tibble(x=seq(10), y=2*x+rnorm(10))
regression_variability_plot(dset, 'y')
regression_variability_plot(dset, 'regression')

##################################################
# Initial exploration: age 20-59 years
##################################################
iset <- bind_rows(grab_and_curate_nhanes(1999, minage=20, maxage=59),
                  grab_and_curate_nhanes(2015, minage=20, maxage=59))

##############################################
# A single binary covariate
##############################################
bmi_density_plot <- function(dset){
    gg <- ggplot(dset)
    gg <- gg+geom_histogram(aes(BMI,
                                y=..density..,
                                colour=YEAR,
                                fill=YEAR),
                            binwidth=1,
                            alpha=0.3,
                            position='identity')
    gg <- gg+theme(legend.position=c(0.95, 0.95),
                   legend.justification=c('right', 'top'))
    gg <- gg+scale_x_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Density',
                                limits=c(0, 0.08),
                                breaks=seq(0, 0.08, by=0.02),
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(name='', discrete=TRUE, begin=0.75, end=0.25)
    gg <- gg+scale_fill_viridis(name='', discrete=TRUE, begin=0.75, end=0.25)
    print(gg)
    # Overall summary
    oset <- dset %>% summarize(N=n(),
                               Mean=mean(BMI),
                               Var=var(BMI),
                               SD=sd(BMI))
    # Year-specific summary
    yset <- dset %>% group_by(YEAR)
    yset <- yset %>% summarize(N=n(),
                               Mean=mean(BMI),
                               Var=var(BMI),
                               SD=sd(BMI))
    # Calculate within-group variance
    wset <- yset %>% summarize(WGVar=sum(Var*N)/sum(N))
    # Alternative calculation
    yset <- dset %>% group_by(YEAR)
    yset <- yset %>% mutate(BMISCALED=scale(BMI, scale=FALSE))
    yset <- yset %>% ungroup()
    wset <- yset %>% summarize(WGVar=var(BMISCALED))
    # Proportion of variance explained
    with(oset, (Var-wset$WGVar)/Var)
}
bmi_density_plot(iset)

##################################################
# Age scatterplot
##################################################
age_scatterplot <- function(dset, correction=0.5){
    gset <- dset %>% group_by(AGE)
    gset <- gset %>% summarize(BMI=mean(BMI))
    gg <- ggplot(dset, aes(x=AGE, y=BMI))
    gg <- gg+geom_point(color='gray', size=1)
    gg <- gg+geom_segment(data=gset,
                          aes(x=AGE-correction,
                              y=BMI,
                              xend=AGE+correction,
                              yend=BMI),
                          size=1)
    gg <- gg+geom_smooth(method='lm', color='black', se=FALSE, size=1)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    print(gg)
    # Compare simple averages and fitted regression values
    fit <- lm(BMI~AGE, data=dset)
    fset <- dset %>% distinct(AGE)
    fset <- fset %>% arrange(AGE)
    fset <- fset %>% mutate(REG=predict(fit, newdata=tibble(AGE)))
    fset <- full_join(gset, fset, by='AGE')
    fset <- bind_cols(fset %>% filter(between(AGE, 20, 29)),
                      fset %>% filter(between(AGE, 30, 39)),
                      fset %>% filter(between(AGE, 40, 49)),
                      fset %>% filter(between(AGE, 50, 59)))
    names(fset) <- gsub('[.0-9]*$', '', names(fset))
    return(fset)
}
age_scatterplot(iset)

##################################################
# Visualize BMI trends by age, sex, and race/ethnicity
##################################################
bmi_trends_plot <- function(dset){
    dset <- dset %>% mutate(AGEGROUP=sub('≤', '""<=', AGEGROUP),
                            AGEGROUP=sub('<=', '""<=', AGEGROUP),
                            AGEGROUP=sub('>', '"">', AGEGROUP),
                            RACE=factor(RACE, labels=sub(' ', '~', levels(RACE))),
                            YEAR=sub('-', '-\n', YEAR))
    mset <- dset %>% group_by(SEX, AGEGROUP, RACE, YEAR)
    mset <- mset %>% summarize(BMI=mean(BMI))
    gg <- ggplot(dset)
    gg <- gg+geom_dotplot(aes(x=YEAR,
                              y=BMI,
                              colour=interaction(SEX, AGEGROUP),
                              fill=interaction(SEX, AGEGROUP)),
                          binaxis='y',
                          binwidth=1.5,
                          stackdir='centerwhole',
                          alpha=0.4)
    gg <- gg+geom_line(data=mset,
                       aes(x=YEAR,
                           y=BMI,
                           group=interaction(SEX, RACE)),
                       size=0.75)
    gg <- gg+facet_grid(RACE~SEX+AGEGROUP, labeller=label_parsed)
    gg <- gg+scale_x_discrete(name='')
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=20),
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(name='', discrete=TRUE, begin=0.75, end=0.25)
    gg <- gg+scale_fill_viridis(name='', discrete=TRUE, begin=0.75, end=0.25)
    print(gg)
}
bmi_trends_plot(iset)

##################################################
# QQ normal plot to assess normality
##################################################
qqnormal_plot <- function(dset){
    gg <- ggplot(dset, aes(sample=Residuals))
    gg <- gg+ggtitle('QQ-normal plot')
    gg <- gg+stat_qq(geom='point',
                     position='jitter',
                     alpha=0.1,
                     size=2)
    gg <- gg+geom_hline(yintercept=0, colour='darkgray', linetype='dashed')
    gg <- gg+scale_x_continuous(name='Theoretical quantiles')
    gg <- gg+scale_y_continuous(name='Residual quantiles',
                                limits=c(-20, 40),
                                expand=c(0, 0))
    print(gg)
}

##################################################
# Residuals vs fitted plot to assess linearity
# and constant variance
##################################################
resvsfit_plot <- function(dset){
    gg <- ggplot(dset)
    gg <- gg+ggtitle('Residuals-versus-fitted-values plot')
    gg <- gg+geom_point(aes_string(x='Fitted',
                                   y='Residuals'),
                        position='jitter',
                        alpha=0.1,
                        size=2)
    gg <- gg+geom_hline(yintercept=0, colour='darkgray', linetype='dashed')
    gg <- gg+scale_x_continuous(name='Fitted')
    gg <- gg+scale_y_continuous(name='Residuals',
                                limits=c(-20, 40),
                                expand=c(0, 0))
    print(gg)
}

##################################################
# Fit linear regression with interaction
##################################################
linear_model <- function(dset,
                         alpha=0.05,
                         include.interaction=FALSE,
                         continuous.age=FALSE,
                         scaled.age=FALSE){
    if(scaled.age)
        dset <- dset %>% mutate(AGE=AGE/10)
    if(include.interaction){
        if(continuous.age){
            equation <- 'BMI~SEX+RACE+AGE*YEAR'
            Label <- 'bmi_trends_continuous_age'
        } else {
            equation <- 'BMI~SEX+RACE+AGEGROUP*YEAR'
            Label <- 'bmi_trends_categorical_age'
        }
    } else {
        if(continuous.age){
            equation <- 'BMI~SEX+RACE+AGE+YEAR'
            Label <- 'bmi_trends_continuous_age_no_interaction'
        } else {
            equation <- 'BMI~SEX+RACE+AGEGROUP+YEAR'
            Label <- 'bmi_trends_categorical_age_no_interaction'
        }
    }
    fit <- lm(as.formula(equation), data=dset)
    cset <- coef(summary(fit))
    cset <- data.frame(Variable=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% transmute(Variable=sub('(Intercept)', 'Intercept', Variable, fixed=TRUE),
                               Variable=gsub(paste0('(', paste(names(dset), collapse='|'), ')'), '\\1 = ', Variable),
                               Variable=sub('GROUP', '', Variable),
                               Variable=sub('.*:.*', 'AGE-YEAR interaction', Variable),
                               Coefficient=Estimate,
                               SE=`Std. Error`,
                               # '95% lower'=Estimate-qnorm(1-alpha/2)*`Std. Error`,
                               # '95% upper'=Estimate+qnorm(1-alpha/2)*`Std. Error`,
                               'P-value'=`Pr(>|t|)`)
    if(continuous.age)
        cset <- cset %>% mutate(Variable=sub('AGE =', 'AGE', Variable))
    if(!include.interaction & !continuous.age){
        rset <- dset %>% mutate(Fitted=fitted(fit),
                                Residuals=resid(fit))
        qqnormal_plot(rset)
        resvsfit_plot(rset)
        # white women young early
        wwye <- matrix(c(1, 1, 0, 0, 0, 0, 0, 0), 1)
        wwye.ht <- multcomp::glht(fit, linfct=wwye)
        confint(wwye.ht)
        # Mexican women old early
        mwoe <- matrix(c(1, 1, 0, 1, 0, 0, 1, 0), 1)
        mwoe.ht <- multcomp::glht(fit, linfct=mwoe)
        confint(mwoe.ht)
        # race using F-test and LRT
        fit.norace <- update(fit, .~.-RACE)
        ftest.race <- anova(fit, fit.norace, test='F')
        lrt.race <- anova(fit, fit.norace, test='LRT')
        cat('Unadjusted R-squared:', summary(fit)$r.squared, '\n')
        cat('Adjusted R-squared:  ', summary(fit)$adj.r.squared, '\n')
    }
    lstats <- tibble(Age=ifelse(continuous.age, 'Continuous', 'Categorical'),
                     Interaction=ifelse(include.interaction, 'Yes', 'No'),
                     Parameters=length(coef(fit))+1,
                     'Log likelihood'=as.numeric(logLik(fit)),
                     AIC=AIC(fit),
                     BIC=BIC(fit))
    return(lstats)
}
cat_ni <- linear_model(iset, include.interaction=FALSE, continuous.age=FALSE)
con_ni <- linear_model(iset, include.interaction=FALSE, continuous.age=TRUE)
cat_wi <- linear_model(iset, include.interaction=TRUE, continuous.age=FALSE)
con_wi <- linear_model(iset, include.interaction=TRUE, continuous.age=TRUE)
lset <- bind_rows(cat_ni, con_ni, cat_wi, con_wi)

##################################################
# Load and filter NHANES demographic, BMI, PE data
##################################################
grab_and_curate_nhanes2 <- function(year, minage=0, maxage=84, minbmi=0, maxbmi=Inf){
    bset <- grab_nhanes('bmx', year)
    dset <- grab_nhanes('demo', year)
    paset <- grab_nhanes('paq',year)
    dset <- dset %>% select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)
    bset <- bset %>% select(SEQN, BMXBMI)
    paset <- paset %>% select(SEQN,PAQ650)
    rset <- right_join(dset, bset, by='SEQN')
    mset <- right_join(rset, paset, by='SEQN')
    mset <- mset %>% rename(ID='SEQN',
                            SEX='RIAGENDR',
                            AGE='RIDAGEYR',
                            RACE='RIDRETH1',
                            BMI='BMXBMI',
                            PE='PAQ650')
    mset <- mset %>% mutate(SEX=factor(SEX,
                                       levels=seq(2),
                                       labels=c('Male', 'Female')),
                            AGEGROUP=factor(AGE <= 50,
                                            levels=c('TRUE', 'FALSE'),
                                            labels=c('<=50', '>50')),
                            AGEDECADE=cut(AGE,
                                          breaks=seq(0, 90, by=10),
                                          labels=paste(seq(0, 80, by=10),
                                                       seq(10, 90, by=10),
                                                       sep='-'),
                                          right=FALSE),
                            RACE=factor(RACE,
                                        levels=c(3, 4, 1, 2, 5),
                                        labels=c('Non-Hispanic White',
                                                 'Non-Hispanic Black',
                                                 'Mexican American',
                                                 'Other Hispanic',
                                                 'Other/Mixed')),
                            PE=ifelse(PE != 1 & PE != 2, NA, PE),
                            PE=factor(PE,
                                      levels=c(2, 1),
                                      labels=c('No', 'Yes')),
                            YEAR=paste(year, year+1, sep='-'))
    mset <- mset %>% filter(between(AGE, minage, maxage),
                            between(BMI, minbmi, maxbmi))
    return(mset)
}

##################################################
# Example of confounding
##################################################
pe_analysis <- function(dset){
    fit0 <- lm(AGE~PE, data=dset)
    cat('AGE~PE:', coef(fit0)['PEYes'], '\n')
    fit1 <- lm(BMI~AGE, data=dset)
    cat('BMI~AGE:', coef(fit1)['AGE'], '\n')
    fit2 <- lm(BMI~PE, data=dset)
    cat('BMI~PE:', coef(fit2)['PEYes'], '\n')
    fit3 <- lm(BMI~PE+AGE, data=dset)
    cat('BMI~PE+AGE:', coef(fit3)['PEYes'], '\n')
}
peset <- grab_and_curate_nhanes2(2015, minage=20, maxage=59)
pe_analysis(peset)

############################################
# Quantile regression
############################################
simple_quantile_regression_plot <- function(dset, correction=0.5){
    q25 <- rq(BMI~AGE, data=dset, tau=0.25)
    q75 <- rq(BMI~AGE, data=dset, tau=0.75)
    pset <- tibble(AGE=seq(20, 59))
    pset <- pset %>% mutate(Predict25=predict(q25, newdata=pset),
                            Predict75=predict(q75, newdata=pset))
    gg <- ggplot(dset)
    gg <- gg+geom_rect(data=pset,
                       aes(xmin=AGE-correction,
                           ymin=Predict25,
                           xmax=AGE+correction,
                           ymax=Predict75),
                       fill='yellow',
                       alpha=0.3)
    gg <- gg+geom_point(aes(x=AGE, y=BMI),
                        color='gray')
    gg <- gg+geom_abline(intercept=coef(q25)[1],
                         slope=coef(q25)[2],
                         size=1)
    gg <- gg+geom_abline(intercept=coef(q75)[1],
                         slope=coef(q75)[2],
                         size=1)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(19, 61),
                                breaks=seq(20, 60, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    print(gg)
}
simple_quantile_regression_plot(iset)

multivariate_quantile_regression <- function(dset, tau=0.5){
    lab <- switch(as.character(tau), '0.25'=' (1st)', '0.75'=' (3rd)')
    fit <- rq(BMI~SEX+RACE+AGE+YEAR, data=dset, tau=tau)
    cset <- coef(summary(fit))
    cset <- data.frame(Variable=rownames(cset), cset, check.names=FALSE)
    cset <- cset %>% transmute(Variable=sub('(Intercept)', 'Intercept', Variable, fixed=TRUE),
                               Variable=gsub(paste0('(', paste(names(dset), collapse='|'), ')'), '\\1 = ', Variable),
                               Variable=sub('AGE =', 'AGE', Variable),
                               Variable=sub('GROUP', '', Variable),
                               Variable=sub('.*:.*', 'AGE-YEAR interaction', Variable),
                               Coef.=Value,
                               SE=`Std. Error`,
                               P=`Pr(>|t|)`)
    names(cset)[-1] <- sub('$', lab, names(cset)[-1])
    return(cset)
}

quartile_regression_table <- function(dset){
    q25 <- multivariate_quantile_regression(dset, tau=0.25)
    q75 <- multivariate_quantile_regression(dset, tau=0.75)
    fset <- bind_cols(q25, q75 %>% select(-Variable))
    return(fset)
}
quartile_regression_table(iset)

############################################
# Nonparametric regression - kernel methods
# sample observations for age 2-84 years for
# pedagogical and numerical reasons
############################################
dset <- bind_rows(grab_and_curate_nhanes(1999),
                  grab_and_curate_nhanes(2015))
sset <- dset %>% sample_n(size=1000, replace=FALSE)

regression_subgroup_plot <- function(dset, group){
    group_sym <- sym(group)
    rset <- dset %>% group_by(!!group_sym)
    rset <- rset %>% summarize(Xmin=min(AGE),
                               Xmax=max(AGE),
                               Mean=mean(BMI))
    panel_title <- switch(group,
                          AGEGROUP='Age ≤50 vs >50 years',
                          AGEDECADE='Age 0-10, ..., 80-90 years')
    gg <- ggplot(dset)
    gg <- gg+ggtitle(panel_title)
    gg <- gg+geom_point(aes(x=AGE, y=BMI),
                        color='gray')
    gg <- gg+geom_segment(data=rset,
                          aes(x=Xmin,
                              y=Mean,
                              xend=Xmax,
                              yend=Mean),
                          size=1)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(0, 86),
                                breaks=seq(0, 80, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    print(gg)
}
regression_subgroup_plot(sset, 'AGEGROUP')
regression_subgroup_plot(sset, 'AGEDECADE')

############################################
# Manually calculate uniform kernel estimate
# across ages using bandwith 10
############################################
halfwidth <- 5
uset <- sset %>% mutate(Lower=AGE-halfwidth, Upper=AGE+halfwidth)
uset <- uset %>% group_by(Lower, Upper)
uset <- uset %>% mutate(BMIUNIFORM=mean(BMI))

uniform_kernel_plot <- function(dset,
                                example_age=42,
                                halfwidth=5,
                                correction=0.5){
    gg <- ggplot(dset)
    gg <- gg+annotate(geom='rect',
                      xmin=example_age-halfwidth-correction,
                      xmax=example_age+halfwidth+correction,
                      ymin=-Inf,
                      ymax=Inf,
                      fill='red',
                      alpha=0.2)
    gg <- gg+geom_point(aes(x=AGE, y=BMI),
                        color='gray')
    gg <- gg+geom_step(aes(x=AGE, y=BMIUNIFORM),
                       size=1)
    gg <- gg+geom_point(aes(x=example_age, y=unique(BMIUNIFORM[AGE == example_age])),
                        colour='darkred',
                        size=4)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(0, 86),
                                breaks=seq(0, 80, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    print(gg)
}
uniform_kernel_plot(uset)

############################################
# Manually calculate normal kernel estimate
############################################
bandwidth <- 3
nset <- sset %>% within({WEIGHTS <- outer(AGE, AGE, function(x, y) dnorm(x-y, sd=bandwidth))
                         BMINORMAL <- apply(WEIGHTS, 1, function(x) sum(BMI*x)/sum(x))
                         rm(WEIGHTS)})

normal_kernel_plot <- function(dset,
                               example_age=42,
                               bandwidth=3,
                               correction=0.5){
    aset <- dset %>% distinct(AGE)
    aset <- aset %>% filter(AGE <= example_age)
    ages <- with(aset, unique(AGE))
    gg <- ggplot(dset)
    for(age in ages){
        example_diff <- abs(example_age-age)
        if(example_diff %in% c(0, 1, 2))
            cat('Diff:', example_diff,
                'Weight:', dnorm(example_diff, sd=bandwidth), '\n')
        gg <- gg+annotate(geom='rect',
                          xmin=example_age-example_diff-correction,
                          xmax=example_age+example_diff+correction,
                          ymin=-Inf,
                          ymax=Inf,
                          fill='red',
                          alpha=dnorm(example_diff, sd=bandwidth))
    }
    gg <- gg+geom_point(aes(x=AGE, y=BMI),
                        color='gray')
    gg <- gg+geom_line(aes(x=AGE, y=BMINORMAL),
                       size=1)
    gg <- gg+geom_point(aes(x=example_age, y=unique(BMINORMAL[AGE == example_age])),
                        colour='darkred',
                        size=4)
    gg <- gg+scale_x_continuous(name='Age (years)',
                                limits=c(0, 86),
                                breaks=seq(0, 80, by=10),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name=parse(text='Body~mass~index~(kg/m^2)'),
                                limits=c(9, 71),
                                breaks=seq(10, 70, by=10),
                                expand=c(0, 0))
    print(gg)
}
normal_kernel_plot(nset)

