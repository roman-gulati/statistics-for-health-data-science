##################################################
# Load or download, extract, read, and save NHANES
##################################################
library(tidyverse)

grab_nhanes <- function(type, year){
    years <- seq(1999, 2015, by=2)
    stopifnot(year %in% years)
    code <- letters[seq_along(years)]
    code[1] <- ''
    names(code) <- years
    year_span <- paste(year, year+1, sep='-')
    remote_sas_file <- paste0('https://wwwn.cdc.gov/nchs/nhanes/',
                              year_span,
                              '/',
                              type,
                              ifelse(year > min(years),
                                     paste0('_', code[[as.character(year)]]),
                                     ''),
                              '.xpt')
    local_sas_file <- paste0(type, '_', year_span, '.xpt')
    local_rdata <- paste0(type, '_', year_span, '.Rdata')
    if(file.exists(local_rdata)){
        load(file=local_rdata)
    } else {
        download.file(remote_sas_file, local_sas_file, mode='wb')
        dset <- foreign::read.xport(local_sas_file)
        save(dset, file=local_rdata)
    }
    return(dset)
}

##################################################
# Load and filter NHANES demographic and BMI data
##################################################
grab_and_curate_nhanes <- function(year, minage=0, maxage=84, minbmi=0, maxbmi=Inf){
    dset <- grab_nhanes('demo', year)
    bset <- grab_nhanes('bmx', year)
    dset <- dset %>% select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)
    bset <- bset %>% select(SEQN, BMXBMI)
    mset <- right_join(dset, bset, by='SEQN')
    mset <- mset %>% rename(ID='SEQN',
                            SEX='RIAGENDR',
                            AGE='RIDAGEYR',
                            RACE='RIDRETH1',
                            BMI='BMXBMI')
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
                            YEAR=paste(year, year+1, sep='-'))
    mset <- mset %>% filter(between(AGE, minage, maxage),
                            between(BMI, minbmi, maxbmi))
    return(mset)
}

