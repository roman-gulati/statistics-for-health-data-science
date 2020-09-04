##################################################
# Load or download, extract, read, and save MEPS
##################################################
grab_meps <- function(year, inscope_only=TRUE){
    code <- paste0('h', seq(12, 201, by=9))
    names(code) <- seq(1996, 2017)
    stopifnot(year %in% names(code))
    year_code <- code[[as.character(year)]]
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
    if(inscope_only){
        suffix <- substr(year, 3, 4)
        inscope <- paste0('INSCOP', suffix)
        dset <- dset[dset[[inscope]] == 1, ]
    }
    return(dset)
}

