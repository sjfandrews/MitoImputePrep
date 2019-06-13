library(tidyverse)
library(readxl)
library(HiMC); data(nodes)
library(dplyr)

'%!in%' = function(x,y)!('%in%'(x,y))

printTime = function(end.time) {
  if (end.time[3] < 60) {
    out.time = end.time[3]
    return(paste0('Process completed in ',round(out.time,2),' seconds'))
  }
  if (end.time[3] >= 60 & end.time[3] < 3600) {
    out.time = end.time[3] / 60
    return(paste0('Process completed in ',round(out.time,2),' minutes'))
  }
  if (end.time[3] >= 3600 & end.time[3] < 86400) {
    out.time = (end.time[3] / 60) / 60
    return(paste0('Process completed in ',round(out.time,2),' hours'))
  }
  if (end.time[3] >= 86400) {
    out.time = ((end.time[3] / 60) / 60) / 24
    return(paste0('Process completed in ',round(out.time,2),' days'))
  }
}

#
start.time = proc.time() # Start the timer!

message("USAGE:")
message("ARGUMENT 1:  PED FILE")
message("ARGUMENT 2:  MAP FILE")
message("ARGUMENT 3:  OUT FILE")

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

ped.file = args[1] # Input PLINK PEDIGREE file
map.file = args[2] # Input PLINK MAP file
out.file = args[3] # Output .csv file

## READ IN .ped AND .map FILES FOR TRUTH SET 
SNP_data = generate_snp_data(map_file = map.file,
                              ped_file = ped.file)
haplogroups = HiMC::getClassifications(SNP_data)

write.csv(haplogroups, out.file, quote = F, row.names = F)

timer = TRUE
if (timer == TRUE) {
  print(printTime(proc.time() - start.time))
}