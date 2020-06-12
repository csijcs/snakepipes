
###############
##
##  Script to rename long name .bam (eg 5905_25_wz3909_TGACTTCG_S35.bam) files from GCF core to wzNUMBER.bam files (wz3909.bam).
##  Deploy script from the directory that your bams are in: setwd("/DATA/your.name/your_bams/")
##  Tesa Severson, 18 May 2020
##
##############

# clean up any leftover R stuff laying around
rm(list = ls())

# make a new directory to put your renamed files in.
subDir <- "rename"
dir.create(file.path(subDir))

# get list of .bam files in folder (long names)
dn = list.files(".", pattern="\\.bam*", full.names = FALSE)

# rename your files with only wz number
newfiles <- list()
tmp <- list()
for(i in 1:length(dn)){
  comp = dn[i]
  x <- unlist(strsplit(dn[i],"_"))#split name
  y <- grep("wz", x,ignore.case=TRUE) # get only 'wz' value out of split
  z <- grep("bam", x)
  wzname <- x[[y]] # find for 'wz' name using value from split
  ext_tmp <- x[[z]]
  ext <- sub("^[^.]*.", "", ext_tmp)
  newfiles <- paste(subDir,"/",wzname,".",ext, sep="") # make new names with the folder 'rename/' first
  file.copy(from=comp, to =newfiles) # save your renamed files into the rename folder
  tmp[i]<-print(paste(comp," saved to ", newfiles))
}
tt <-as.data.frame(do.call(rbind,tmp))
write.table(tt, file="rename_log.txt",row.names = F,col.names = F)
