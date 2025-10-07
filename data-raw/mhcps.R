# MHCPS processing script
# Original data source:
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.0006-341X.2002.00064.x&file=BIOM_64_sm_010423.txt
# Also in data-raw/mhcps_raw.txt
# In data-raw/mhcps.txt some corrections have been made:
# Line 86 had two merged values:
# 68.378.3 9999 1 has been changed to 68.3 78.3 9999 1
# Lines 321 and 555 appear to have two lines merged in each.
# These have been separated.
mhcps <- read.delim("data-raw/mhcps.txt", sep = " ")
# Deal with right-censoring
mhcps$Vi[mhcps$Vi == 9999] <- Inf
# Values that appear to be incorrect
mhcps$Ui[mhcps$Ui == 7755] <- 77.55
mhcps$Ui[mhcps$Ui == 873] <- 87.3
mhcps$Ui[mhcps$Ui == 1.55] <- 71.55
mhcps$Ui[mhcps$Ui == 3.3] <- 83.3
mhcps$Ti[mhcps$Ti == 689] <- 68.9
# Save
save(mhcps, file = "data/mhcps.rda")
