
###############
##
##  Script to rename long name .bam (eg 5905_25_wz3909_TGACTTCG_S35.bam) files from GCF core to wzNUMBER.bam files (wz3909.bam). 
##  Deploy script from the directory that your bams are in: setwd("/DATA/your.name/your_bams/")
##  requires library magicfor. 
##  Tesa Severson, 18 May 2020
##
##############

# clean up any leftover R stuff laying around
rm(list = ls())

# make a new directory to put your renamed files in. 
subDir <- "rename"
dir.create(file.path(subDir))

# get list of .bam files in folder (long names)
dn = list.files(".", pattern="\\.bam$", full.names = FALSE)

install.packages("magicfor")
library(magicfor)               # library to take print function and make a file
magic_for(print, silent = FALSE) # call magic_for()
# rename your files with only wz number
newfiles <- list()
for(i in 1:length(dn)){
  comp = dn[i]
  x <- unlist(strsplit(dn[i],"_"))#split name
  y <- grep("wz", x) # get only 'wz' value out of split
  wzname <- x[[y]] # find for 'wz' name using value from split
  newfiles <- paste("rename/",wzname,".bam", sep="") # make new names with the folder 'rename/' first
  file.copy(from=comp, to =newfiles) # save your renamed files into the rename folder
  print(paste(comp," saved to ", newfiles))
}
tt <-magic_result_as_dataframe()
write.table(tt, file="rename_log.txt")
