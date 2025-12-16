####################################
####                            ####
#### 4) SENTINEL SET EVALUATION ####
####    appendix                ####
####                            ####
####################################

#### append an ordered set of facilities to a fIDsubsets object 
## enables evaluation of a sentinel set in a non-native parameter context
## NB if time-shifting beware use of same code to represent different pre/post merger entities 

### run section-by-section to ENTER values

library(tidyverse)
library(readr)

setwd()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> fIDsubsets object to which the set will be appended

## ENTER file_id1 for the network inherited by the fIDsubsets object 
network_file_id1 <- "ALL-fy2000-wdw182-cw03"

# check exists
stopifnot (file.exists(paste0(paste("fIDsubsets", network_file_id1, sep = "_"), ".rds")))
stopifnot (file.exists(paste0(paste("hnreallife", network_file_id1, sep = "_"), ".csv")))

## read in fIDsubsets, view the existing entries
fIDsubsets <- read_rds(paste0(paste("fIDsubsets", network_file_id1, sep = "_"), ".rds"))
sapply(fIDsubsets, "[[", "txt")
sapply(fIDsubsets, "[[", "fIDvector")


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> (ordered) set of facilities to append to the above
#> create addset = "fID" vector for the ordered set of facilities to add
#> "fID" values must be those used in the HospiNet which generated the fIDsubset
#> THIS MAY REQUIRE MANUAL EYEBALL / ENTRY to ensure real life validity

## ENTER addset = the vector of fID values (in order)
# example here from an optimised priority list with given file_id3
addset_file_id3 <- "DUMMY-fy2020-wdw182-cw03_p00005-sp10-i100_meanBHG10-excSpl-noCap-noFrc"
stopifnot (file.exists(paste0(paste("sentinelcsv", addset_file_id3, sep = "_"), ".csv")))
addset <- read_csv(paste0(paste("sentinelcsv", addset_file_id3, sep = "_"), ".csv")) %>% pull(fID)
# or manual entry values to addset 

## enter freetext values to be used to identify this subset
# Fulltxt is eyeball-friendly description NB suitable for list element name (cf "Specialist")
# Txt is shorthand used in file-naming convention and Step3 parameterisation (cf "Spl")
addset_fulltxt <- "Benchmark"
addset_txt <- "Ben"

## utility function to confirm if a priority list and fIDsubset share matching fID real life information
fn_shared_fID <- function(prioritylist_id3, fIDsubsets_id1) {
  prioritylist_id1 = str_split_i(prioritylist_id3, pattern = "_", i = 1)
  stopifnot (file.exists(paste0(paste("hnreallife", prioritylist_id1, sep = "_"), ".csv")))
  stopifnot (file.exists(paste0(paste("hnreallife", fIDsubsets_id1, sep = "_"), ".csv")))
  if (all.equal(read_csv(paste0(paste("hnreallife", prioritylist_id1, sep = "_"), ".csv")) %>% select(-hnfile, -hnrow),
                read_csv(paste0(paste("hnreallife", fIDsubsets_id1, sep = "_"), ".csv")) %>% select(-hnfile, -hnrow))) {
    print("Matching real life information for fID values")} else {print("MANUAL CHECK REQUIRED")}
  }
fn_shared_fID(prioritylist_id3 = addset_file_id3, fIDsubsets_id1 = network_file_id1)

    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> edit fIDsubset

## check and append 
if (length(setdiff(addset, read_csv(paste0(paste("hnreallife", network_file_id1, sep = "_"), ".csv"))$fID)) > 0) {stop("addset is not a subset of the network facilities")}
if (addset_fulltxt %in% names(fIDsubsets)) {stop(paste(addset_fulltxt, "is already used as a name"))}
if (addset_fulltxt %in% sapply(fIDsubsets, "[[", "txt")) {stop(paste(addset_fulltxt, "is already used as a shorthand"))}
fIDsubsets[[addset_fulltxt]] <- list(txt = addset_txt, fIDvector = addset)
# view
sapply(fIDsubsets, "[[", "txt")
fIDsubsets[[addset_fulltxt]]$fIDvector


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> save appended fIDsubset

write_rds(fIDsubsets,
          file = paste0(paste("fIDsubsets", file_id1, sep = "_"), ".rds"))

