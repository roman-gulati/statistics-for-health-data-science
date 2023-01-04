##################################################
# Load or download, extract, read, and save MEPS
##################################################
grab_meps <- function(year, inscope_only=TRUE){
    code <- paste0('h', seq(12, 201, by=9))
    names(code) <- seq(1996, 2017)
    stopifnot(year %in% names(code))
    year_code <- code[[as.character(year)]]
    local_rdata <- paste(year_code, 'Rdata', sep='.')
    if(file.exists(local_rdata)){
        load(file=local_rdata)
    } else {
        if(year < 2017){
            archive_filename <- paste0(year_code, 'ssp.zip')
            remote_archive <- paste0('https://meps.ahrq.gov/mepsweb/data_files/pufs/',
                                     archive_filename)
            local_archive <- paste0(year_code, 'ssp.zip')
            old_timeout <- getOption('timeout')
            options(timeout=max(300, old_timeout))
            download.file(remote_archive, local_archive, mode='wb')
            options(timeout=old_timeout)
            dset <- foreign::read.xport(unzip(local_archive))
            local_file <- sub('ssp.zip$', '.ssp', local_archive)
        } else {
            archive_filename <- paste0(year_code, 'dta.zip')
            remote_archive <- paste0('https://meps.ahrq.gov/data_files/pufs/',
                                     year_code,
                                     '/',
                                     archive_filename)
            local_archive <- paste0(year_code, 'dta.zip')
            old_timeout <- getOption('timeout')
            options(timeout=max(300, old_timeout))
            download.file(remote_archive, local_archive, mode='wb')
            options(timeout=old_timeout)
            dset <- haven::read_dta(unzip(local_archive))
            local_file <- sub('dta.zip$', '.dta', local_archive)
        }
        file.copy(local_file, local_rdata, overwrite=TRUE)
        unlink(local_file)
        save(dset, file=local_rdata)
    }
    if(inscope_only){
        suffix <- substr(year, 3, 4)
        inscope <- paste0('INSCOP', suffix)
        dset <- dset[dset[[inscope]] == 1, ]
    }
    return(dset)
}

