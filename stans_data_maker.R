library(tidyverse)
library(data.table)
library(pbapply)

pmppm <- function(mass, ppm=4){
  if(mass<200){
    as.numeric(mass)+(c(-ppm, ppm)*200/1000000)
  } else {
    c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
  }
}
grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  
  return(all_data)
}

# RML file import suggestions --------------------------------------------------

## Consider using something like the below code I added in- you always want to try and avoid 
# desktop-specific paths for data input. You can use whatever works for you, especially as you're
# going through your Google Drive, but just consider something beyond a normal file path.

# An example of how to use the below function and loop is in the Targeted Pipeline project on the 
# Github.
RemoveCsv <- function(full.filepaths) {
  # Gathers all files in given directory and drops the csv extension.
  #
  # Args
  #   full.filepaths: list of files in a directory matching given patterns.
  #
  # Returns
  #   no.path: list of files, with filepath and csv extension removed.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

filenames <- RemoveCsv(list.files(path = "data_raw", pattern = file.pattern))

for (i in filenames) {
  filepath <- file.path("data_raw", paste(i,".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = TRUE))
}
# --------------------------------------------------

ms_file_paths <- dir("(G:\My Drive\FalkorFactor\mzMLs\pos)",
                     pattern = ".mzML",
                     full.names = TRUE) %>%
  normalizePath() %>% # This is a cool function, haven't heard of it before.
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .)) %>%
  `[`(!grepl("DDApos", x = .)) %>%
  `[`(!grepl("180205", x = .)) %>%
  `[`(grepl("Std", x=.))

all_stans <- read.csv("falkor_stans.csv") %>%
  mutate(mz=gsub("Â", "", .$mz)) %>%
  mutate(rt=as.numeric(rt)) %>%
  mutate(mz=as.numeric(mz)) %>%
  filter(polarity=="pos") %>% 
  arrange(mz, rt)

stan_MS1_data <- pblapply(ms_file_paths, function(file_path){
  raw_EIC <- grabSingleFileData(file_path) %>%
    as.data.table()
  v <- lapply(unique(all_stans$mz), function(mz_i){
    raw_EIC[mz%between%pmppm(as.numeric(mz_i), 5)]
  }) %>% do.call(what = rbind) %>%
    cbind(fileid=basename(file_path))
}) %>% do.call(what = rbind)

saveRDS(stan_MS1_data, file = "stan_data.rds")
