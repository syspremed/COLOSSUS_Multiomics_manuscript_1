################################################################################
# Figure Panel : 5C
# Description  : Boxplot and t-test of MV5 scores vs perineural invasion (PNI) in COLOSSUS St2/3
#                MSS RAS-mt CRC (n = 98).
# Input data   : 20260502_COLOSSUS_new_MV5_perineural_invasion.txt
# Output       : boxplot (screen) and t-test result (console)
# R packages   : (base R)
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : MV5 is multiplied by -1 so that higher values correspond to the PNI direction.
################################################################################

# perineural vs MV5 

data <- read.delim("20260502_COLOSSUS_new_MV5_perineural_invasion.txt", stringsAsFactors = FALSE)
dim(data)

boxplot(-1*data$MV5~data$perineural_invasion)
t.test(-1*data$MV5~data$perineural_invasion)


