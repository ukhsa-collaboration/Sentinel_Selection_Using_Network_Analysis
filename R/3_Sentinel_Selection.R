###############################
####                       ####
#### 3) SENTINEL SELECTION ####
####                       ####
###############################


#### reads in previous step output 
#### creates priority list via greedy algorithm, optimised vs specified protocol, within exclude/cap/force restrictions on potential sentinels
#### creates and saves sentinelcsv[csv], criteriaCribsheet[rds]  
#### NB 'times' in the sentinelcsv file are NOT the run-by-run times used for set evaluation

library(tidyverse)
library(readr)
library(tictoc)


setwd()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> SET PARAMETERS for step_3 : Sentinel Selections
#> also creates file prefix for output

parameters3 <- list(
  inherit_spread = list(txt = file_id2),   # use file_id2 for the spread simulation
  protocol = list(smp_code = "BHG10"),     # optimise on this protocol's detection time, use smp_code (NA for 'sequence all cases') 
  fID_exclude = list(fIDset = "Spl"),      # subset to exclude from priority list, use shorthand from fIDsubsets for this network (NA for no excluded set)
  fID_cap = list(fIDset = NA,              # subset to restrict use in the priority list, use shorthand from fIDsubsets for this network (NA for no capped set)
                 value = NA),              # max number of restricted set permitted in the priority list (ignored if fIDset = NA)
  fID_force = list(fIDset = NA)            # subset to force inclusion at top of priority list, use shorthand from fIDsubsets for this network (NA for no forced set)   
  )

# parameter checking and file-name suffix generation 
parameters3$inherit_network <- list(txt = str_split_i(parameters3$inherit_spread, "_", i = 1))
stopifnot (sapply(FUN = file.exists, str_replace(c("hnmatrix_NETWORKid.rds", "hnreallife_NETWORKid.csv", "fIDsubsets_NETWORKid.rds"),
                                                 pattern = "NETWORKid", replacement = parameters3$inherit_network$txt)))
stopifnot (sapply(FUN = file.exists, str_replace(c("sampleAffectedTimes_SPREADid.csv", "samplingCodes_SPREADid.csv"), 
                                                 pattern = "SPREADid", replacement =  parameters3$inherit_spread$txt)))
parameters3$protocol$obs <- ifelse(is.na(parameters3$protocol$smp_code), "atimes", paste0("dTimes_", parameters3$protocol$smp_code))
parameters3$protocol$txt <- paste0("mean", ifelse(is.na(parameters3$protocol$smp_code), "atimes", parameters3$protocol$smp_code))
if (!parameters3$protocol$smp_code %in% read_csv(paste0("samplingCodes_", parameters3$inherit_spread, ".csv"))$smp_code) {stop("STOP this protocol is not in the spread data")}
if (!all(sapply(parameters3[str_starts(names(parameters3), "fID")], "[[", "fIDset") %in% c(NA, sapply(read_rds(paste0(paste("fIDsubsets", parameters3$inherit_network, sep = "_"), ".rds")), "[[", "txt")))) {
  stop("STOP cannot find subset definitions wrt inherited network - to append see example in step4")}
if (!is.na(parameters3$fID_cap$fIDset)) {
  if (!(sign(parameters3$fID_cap$value) %in% c(1))) {stop("STOP require a number for the capped subset")}}
parameters3$fID_exclude$txt <- ifelse(is.na(parameters3$fID_exclude$fIDset), "noExc", paste0("exc", parameters3$fID_exclude$fIDset))
parameters3$fID_cap$txt <- ifelse(is.na(parameters3$fID_cap$fIDset), "noCap", paste0("cap", parameters3$fID_cap$fIDset, parameters3$fID_cap$value))
parameters3$fID_force$txt <- ifelse(is.na(parameters3$fID_force$fIDset), "noFrc", paste0("frc", parameters3$fID_force$fIDset))

file_id3 <- paste(parameters3$inherit_spread, 
                  paste(sapply(parameters3[which(!str_starts(names(parameters3), "inherit"))], "[[", "txt"), collapse = "-"),
                  sep = "_")
file_id3


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> READ IN FROM PREVIOUS STEP(S)

hnmatrix <- read_rds(paste0(paste("hnmatrix", parameters3$inherit_network, sep = "_"), ".rds"))
hnreallife <- read_csv(paste0(paste("hnreallife", parameters3$inherit_network, sep = "_"), ".csv"))
fIDsubsets <- read_rds(paste0(paste("fIDsubsets", parameters3$inherit_network, sep = "_"), ".rds"))
samplingCodes <- read_csv(paste0(paste("samplingCodes", parameters3$inherit_spread, sep = "_"), ".csv"))
sampleAffectedTimes <- read_csv(paste0(paste("sampleAffectedTimes", parameters3$inherit_spread, sep = "_"), ".csv"), 
                                col_types = cols(truncated = col_number()))

print(paste(parameters3$inherit_network$txt, "has", nrow(hnmatrix), "fIDs"))
if (!(parameters3$protocol$obs %in% names(sampleAffectedTimes))) {stop("STOP check spread data variable names")}
print(paste(parameters3$inherit_spread$txt, "has", 
      length(which(sampleAffectedTimes[[parameters3$protocol$obs]] < 0)),
      "/", nrow(sampleAffectedTimes), "truncations before complete", 
      paste0("(", str_remove(parameters3$protocol$txt, "mean"), ")")))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> PREPARE criteria subset files used by greedy algorithm
#> also expected length of priority list
#> uses information from parameters3

## vectors of fID to exclude/cap/force in priority list
# donotchoose (exclude)
# maxNchoose (only Ncap may be used)
# forcechoose (include first before greedy algorithm kicks in for remaining priorities)
donotchoose <- "STOP"
if (parameters3$fID_exclude$txt == "noExc") {donotchoose <- NULL} else {
  if (parameters3$fID_exclude$fIDset %in% sapply(fIDsubsets, "[[", "txt")) {
    donotchoose <- fIDsubsets[[which(sapply(fIDsubsets, "[[", "txt") %in% parameters3$fID_exclude$fIDset)]]$fIDvector
    }
  }
donotchoose <- unique(donotchoose)
maxNchoose <- "STOP"
Ncap <- "STOP"
if (parameters3$fID_cap$txt == "noCap") {
  maxNchoose <- NULL
  Ncap <- nrow(hnmatrix) +1
  } else {
    if (parameters3$fID_cap$fIDset %in% sapply(fIDsubsets, "[[", "txt")) {
      maxNchoose <- fIDsubsets[[which(sapply(fIDsubsets, "[[", "txt") %in% parameters3$fID_cap$fIDset)]]$fIDvector
      Ncap <- parameters3$fID_cap$value
    }
  }
maxNchoose <- unique(maxNchoose)
forcechoose <- "STOP"
if (parameters3$fID_force$txt == "noFrc") {forcechoose <- NULL} else {
  if (parameters3$fID_force$fIDset %in% sapply(fIDsubsets, "[[", "txt")) {
    forcechoose <- fIDsubsets[[which(sapply(fIDsubsets, "[[", "txt") %in% parameters3$fID_force$fIDset)]]$fIDvector
  }
}
forcechoose <- unique(forcechoose)                   

## check validity of these criteria and in combination : stop if impossible or undefined
# vs fID in the hnmatrix only
if (any(donotchoose == "STOP", maxNchoose == "STOP", Ncap == "STOP", forcechoose == "STOP")) {stop("STOP undefined criteria shorthand")}
if (!is_empty(
  intersect(intersect(rownames(hnmatrix), donotchoose), intersect(rownames(hnmatrix), forcechoose))
)) {stop("STOP overlap in forced and excluded facilities")}
if (n_distinct(
  intersect(intersect(rownames(hnmatrix), forcechoose), intersect(rownames(hnmatrix), maxNchoose))) 
  >= Ncap) {stop("STOP forced facilities include too many capped facilities")}

## corresponding vectors with hnrow replacing fID
donotchoose_index <- which(rownames(hnmatrix) %in% donotchoose)
maxNchoose_index <- which(rownames(hnmatrix) %in% maxNchoose)
if (length(forcechoose) == 0) {forcechoose_index <- which(rownames(hnmatrix) %in% forcechoose)} else {
  forcechoose_index <- rep(NA_integer_, length(forcechoose))
  forcechoose_index <- sapply(1:length(forcechoose_index),
                              FUN = function(i) {forcechoose_index[i] = which(rownames(hnmatrix) %in% forcechoose[i]) })
} # function ensures forcechoose_index matches order of forcechoose
# doublecheck hnrow <-> fID, as used in hnreallife
stopifnot(setequal(hnreallife$fID[hnreallife$hnrow %in% donotchoose_index], donotchoose),
          setequal(hnreallife$fID[hnreallife$hnrow %in% maxNchoose_index], maxNchoose),
          setequal(hnreallife$fID[hnreallife$hnrow %in% forcechoose_index], forcechoose))

## save in single list
sentinel_criteria <- list(set_exclude = data.frame("fID" = donotchoose, "hnrow" = donotchoose_index),
                          set_cap = data.frame("fID" = maxNchoose, "hnrow" = maxNchoose_index),
                          Ncap = Ncap,
                          set_force = data.frame("fID" = forcechoose, "hnrow" = forcechoose_index))

## expected priority list length
n_sentinels  = 
  n_distinct(setdiff(rownames(hnmatrix), sentinel_criteria$set_exclude$fID)) -
  ifelse(is.null(sentinel_criteria$set_cap$fID), 0,
         n_distinct(intersect(rownames(hnmatrix), sentinel_criteria$set_cap$fID)) - unique(sentinel_criteria$Ncap))
print(paste(n_sentinels, "fIDs expected on priority list"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> PREPARE spread times used by greedy algorithm
#> source/target summary value for nominated measure
#> uses information from parameters3

## tidied copy of sampleAffectedTimes
# obsTimes : observation times to be optimised
subsampleAffectedTimes2 <- sampleAffectedTimes
subsampleAffectedTimes2$obsTimes <- subsampleAffectedTimes2[[parameters3$protocol$obs]]
# replace incomplete obsTimes with truncated value
obs_cutoff <- suppressWarnings(max(subsampleAffectedTimes2$truncated, na.rm = TRUE))
if (any(subsampleAffectedTimes2$obsTimes < 0)) {
  if (!is.finite(obs_cutoff)) {stop("STOP mismatch between trunction flag and incomplete runs")}
  subsampleAffectedTimes2$obsTimes[subsampleAffectedTimes2$obsTimes < 0] <- obs_cutoff
 }

## Dataframe with source/target/cmbT
# calculate the mean observation for each source/target pair
ComboObsTimes <- subsampleAffectedTimes2 %>% group_by(target, source) %>%
  summarise(cmbT = mean(obsTimes), .groups = "drop") %>% as.data.frame()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  SELECTION ALGORITHM FUNCTIONS

  
## Calculate expected detection time given a sentinel set
# sentinelSet = facilities to use as sentinels [vector: hnrow value]
# obsTable = dataframe with source/target/cmbT [one row per source/target pair, hnrow values]
# function returns: detection time [single value]
calcTdetect <- function(sentinelSet, obsTable) {
  mean((obsTable[obsTable$target %in% sentinelSet, ] %>% 
          group_by(source) %>%
          summarise(firstDetect = min(cmbT), .groups = "drop") %>%
          arrange(firstDetect) %>% as.data.frame())$firstDetect)
  } #end calcTdetect

## Find next sentinel to add to priority list - single iteration using GREEDY ALGORITHM
# calls calcTdetect passing observationsTable
# observationsTable = dataframe with source/target/cmbT, optimises to minimise cmbT [one row per source/target pair, hnrow values]
# includedSentinels = priority list so far
# candidateRestrictions = facilities to exclude/cap/force [list with dataframes set_exclude, set_cap, set_force each with hnrow variable, Ncap]  
# function returns: list of results []
findNextSentinel <- function(observationsTable = ComboObsTimes,
                             includedSentinels, 
                             candidateRestrictions = sentinel_criteria){
  
  unavailableFacilities = candidateRestrictions$set_exclude$hnrow
  cappedFacilities = candidateRestrictions$set_cap$hnrow
  cappedN = Ncap
  forceFacilities = candidateRestrictions$set_force$hnrow
  if (setequal(includedSentinels, forceFacilities)) {print(paste(length(includedSentinels), "forced sentinels added"))} # onscreen commentary

  # reset
  forcedChoice <- FALSE
  tiedOptimal <- FALSE
  capsIncluded <- n_distinct(intersect(includedSentinels, cappedFacilities))
  allHospitals <- 1:max(observationsTable$target)
  
  # next sentinel to be chosen from these facilities
  # exclude facilities already on list or in fID_exclude
  # exclude fID_cap facilities if cap value already reached
  # if fID_force incomplete, take next fID_force as only candidate
  if (length(includedSentinels) == 0) {print(paste(length(intersect(allHospitals, unavailableFacilities)), "facilities excluded"))} # onscreen commentary
  trySentinels <- intersect(allHospitals[!allHospitals %in% includedSentinels],
                            allHospitals[!allHospitals %in% unavailableFacilities]) 
  if (length(intersect(includedSentinels, cappedFacilities))>(cappedN-1)) {
    trySentinels <- setdiff(trySentinels, cappedFacilities)
    if (tail(includedSentinels, n = 1) %in% cappedFacilities) {print(paste("list #", length(includedSentinels), "capped reached", paste0("(", cappedN, ")")))} # onscreen commentary
    } 
  candidateCount <- length(trySentinels)
  #print(paste("trySentinels remaining length = ", candidateCount, "list so far", length(includedSentinels)))
  if (length(setdiff(forceFacilities, includedSentinels)) > 0) {
    #print (paste("force loop, adding", setdiff(forceFacilities, includedSentinels)[1])) # onsrceen commentary
    forcedChoice <- TRUE
    trySentinels <- setdiff(forceFacilities, includedSentinels)[1]
  }
  
  # calculate times for all trySentinels and select facility(s) with minimum time
  newTDetects <- unlist(lapply(trySentinels, function(h){
    includedSentinelsTemp <- c(includedSentinels, h)
    calcTdetect(sentinelSet = includedSentinelsTemp, obsTable = observationsTable)
    }))
  nextSentinel <- which(newTDetects == min(newTDetects))
  # break ties
  if(length(nextSentinel) > 1){
    nextSentinel<-sample(nextSentinel, 1)
    print(paste("list #", 1+length(includedSentinels), "randomly selected from equal best options")) # commentary
    tiedOptimal <- TRUE
  }
  
  # return list with values/flags for selected facility
  list(
    tDetect = newTDetects[nextSentinel],       # time
    sentinel = trySentinels[nextSentinel],     # facility number
    tied = tiedOptimal,                        # logical flag for resolved tie
    forced = forcedChoice,                     # logical flag for forced inclusion
    capping = trySentinels[nextSentinel] %in% cappedFacilities,
                                               # capped set included
    capinc = capsIncluded + (trySentinels[nextSentinel] %in% cappedFacilities)
  )
} #end findNextSentinel


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  RUN THE SENTINEL SELECTION
#>  use default values in findNextSentinel ie ComboObsList and sentinel_criteria

## set up empty list
includedSentinelList <- c()
resultsDataFrame <- NULL
## list interation
tic()
for (it in 1:(n_sentinels)) {
  resSentinel <- findNextSentinel(includedSentinels = includedSentinelList)
  resultsDataFrame <- bind_rows(resultsDataFrame, resSentinel)
  includedSentinelList <- c(includedSentinelList, resSentinel$sentinel)
}
toc()
print(paste(sum(resultsDataFrame$tied), "ties broken ( highest tied rank was", min(which(resultsDataFrame$tied)), ")")) # onscreen commentary
print(paste("Priority list length", nrow (resultsDataFrame)))

## results in real life
resultsDataFrame <- left_join(resultsDataFrame, hnreallife, by = c("sentinel" = "hnrow")) %>% 
  mutate(sentinelrank = row_number()) %>% relocate(sentinelrank) %>% 
  rename(hnrow = sentinel) %>% relocate(hnrow, .before = "fID")
head(resultsDataFrame)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> SAVE in working directory using file-naming convention

write_csv(resultsDataFrame,
          file = paste0(paste("sentinelcsv", file_id3, sep = "_"), ".csv"))
write_rds(sentinel_criteria,
          file = paste0(paste("criteriaCribsheet", file_id3, sep = "_"), ".rds"))

