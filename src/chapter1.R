##################################################
# Create figures for Chapter 1
##################################################
library(dplyr)
library(readr)
library(ggplot2)
library(grid)
library(viridis)

source('shared.R')

##################################################
# Visualize SEER incidence of cervical, melanoma,
# and prostate cancer
##################################################
seer_incidence <- function(filename, ext='pdf', saveit=FALSE){
    dset <- read_csv(here('data', filename), col_types='ciddd')
    dset <- dset %>% rename(Site='Site recode ICD-O-3/WHO 2008 (PCM)',
                            Year='Year of diagnosis (1975-2017)',
                            Rate='Age-Adjusted Rate')
    gg_theme()
    gg <- ggplot(data=dset)
    gg <- gg+geom_line(aes(Year, Rate, colour=Site),
                       size=1)
    gg <- gg+geom_text(data=dset %>% filter(Year == 2007),
                       aes(Year,
                           Rate,
                           label=Site,
                           colour=Site),
                       vjust=-1.5,
                       size=7)
    gg <- gg+scale_colour_viridis(discrete=TRUE, option='cividis')
    gg <- gg+scale_x_continuous(name='',
                                limits=c(1974, 2018),
                                breaks=seq(1975, 2015, by=5),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Diagnoses per 100,000\n',
                                limits=c(0, 100),
                                breaks=seq(0, 100, by=20),
                                expand=c(0, 0))
    print(gg)
    if(saveit){
        filename <- '01-seer_incidence'
        filename <- paste(filename, ext, sep='.')
        ggsave(plot=gg,
               file=here('figures', filename),
               height=5,
               width=10)
    }
}
#seer_incidence('inc_cervix_melanoma_prostate_1975-2017.csv', saveit=TRUE)

##################################################
# Load MEPS 2017 data
##################################################
grab_and_format_meps <- function(year){
    dset <- grab_meps(2017)
    dset <- dset %>% mutate(SCALEDTOTEXP17=TOTEXP17/10000,
                            LOGTOTEXP17=log(TOTEXP17+1))
    dset <- dset %>% select(SCALEDTOTEXP17, LOGTOTEXP17)
    return(dset)
}
#dset <- grab_and_format_meps(2017)

##################################################
# Visualize histogram of total expenditures
##################################################
totexp_panel <- function(dset, varname, xmax, xstep, ymax, ystep){
    if(varname == 'SCALEDTOTEXP17'){
        xlabel <- '\nTotal expenditure ($10,000)'
        bpoints <- c(seq(0, xmax, length=15), max(dset[['SCALEDTOTEXP17']]))
    } else {
        xlabel <- '\nLog(Total expenditure+$1)'
        bpoints <- 20
    }
    hdat <- hist(dset[[varname]], breaks=bpoints, plot=FALSE)
    hset <- with(hdat, data.frame(mids, counts))
    hset <- hset %>% mutate(scaledcounts=counts/1e3)
    gg_theme()
    gg <- ggplot(hset)
    gg <- gg+geom_bar(aes(x=mids, y=scaledcounts), stat='identity')
    gg <- gg+scale_x_continuous(name=xlabel,
                                breaks=seq(0, xmax, by=xstep),
                                limits=c(0, xmax))
    gg <- gg+scale_y_continuous(name='Participants (1,000)\n',
                                breaks=seq(0, ymax, by=ystep),
                                limits=c(0, ymax))
    return(gg)
}

totexp_plot <- function(dset, ext='pdf', saveit=FALSE){
    spanel <- totexp_panel(dset, 'SCALEDTOTEXP17', xmax=10, xstep=1, ymax=35, ystep=5)
    lpanel <- totexp_panel(dset, 'LOGTOTEXP17', xmax=16, xstep=2, ymax=6, ystep=1)
    Layout <- grid.layout(ncol=2, nrow=1, heights=unit(5, 'null'))
    if(saveit){
        filename <- '01-natural_and_transformed_expenditures'
        filename <- paste(filename, ext, sep='.')
        get(ext)(here('figures', filename),
                 height=5,
                 width=10)
        on.exit(graphics.off())
    }
    grid.newpage()
    pushViewport(viewport(layout=Layout))
    print(spanel, vp=viewport(layout.pos.col=1, layout.pos.row=1))
    print(lpanel, vp=viewport(layout.pos.col=2, layout.pos.row=1))
}
#totexp_plot(dset, saveit=TRUE)

