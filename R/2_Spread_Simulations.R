###############################
####                       ####
#### 2) SPREAD SIMULATIONS ####
####                       ####
###############################

#### reads in previous step output 
#### applies overall case incidence factor
#### simulates novel strain spread
#### for each simulation: captures affected times and detection times for nominated sampling protocol(s)
#### creates and saves sampleAffectedTimes[csv], samplingCodes[rds]

library(tidyverse)
library(readr)
library(janitor)
library(tictoc)

setwd()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> SET PARAMETERS for step_2 : Spread Simulations
#> also creates file prefix for output

parameters2 <- list(
  inherit_network = list(txt = file_id1),                 # use file_id1 for the (weighted) network
  transmission_parameter = list(value = 0.00001),         # file_id takes value after decimal point
  protocol_list = list(value = 10),                       # monthly sequenced specimens (max) per facility (NB "sequence all cases" included as a default)
  incidence_multiplier = list(value = 1)                  # file_id takes nearest %age
  )

# additional information for spread simulation
max_timestep <- 20*365          # end incomplete simulation at this timestep
sims_per_seed <- 100            # simulations per seed location (currently uses same value for all locations)


# parameter checking and file-name suffix generation 
stopifnot (sapply(FUN = file.exists, str_replace(c("hnmatrix_NETWORKid.rds", "hnreallife_NETWORKid.csv", "hndenominatorcases_NETWORKid.csv"), 
                                                 pattern = "NETWORKid", replacement = parameters2$inherit_network$txt) ))
parameters2$transmission_parameter$txt <- paste0("p", str_remove(format(parameters2$transmission_parameter$value, sci = 10), "0."))
parameters2$incidence_multiplier$txt <- paste0("i", str_pad(round_half_up(parameters2$incidence_multiplier$value*100), width = 3, side = "left", pad =  "0"))
parameters2$protocol_list$txt = paste0("sp", str_pad(parameters2$protocol_list$value, width = 2, side = "left", pad = "0"))
file_id2 <- paste(parameters2$inherit_network, 
                  paste(sapply(parameters2[-1], "[[", "txt"), collapse = "-"),
                  sep = "_")
file_id2


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  READ IN FROM PREVIOUS STEP(S)

hnmatrix <- read_rds(paste0(paste("hnmatrix", parameters2$inherit_network$txt, sep = "_"), ".rds"))
hnreallife <- read_csv(paste0(paste("hnreallife", parameters2$inherit_network$txt, sep = "_"), ".csv"))
hndenominatorcases <- read_csv(paste0(paste("hndenominatorcases", parameters2$inherit_network$txt, sep = "_"), ".csv"))

print(paste(parameters2$inherit_network$txt, "has", nrow(hnmatrix), "fIDs"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  PREPARE files used in the singleRun function
#>  uses information from parameters2

## calculate mean monthly total cases for protocol denominator - including incidence factor
# vector in same order as hnrow
Nmonth_mean <- hndenominatorcases %>%
  mutate(annualcases_iadjusted = fID_meanannualtotalcases*parameters2$incidence_multiplier$value,
         monthlycases_iadjusted = pmax(1, janitor::round_half_up(annualcases_iadjusted/12))) %>%
  arrange(hnrow) %>% pull(monthlycases_iadjusted)

## sampling_codes table
# variables used in simulations: smp_fn (detection function), smp_spec (specimens per month), smp_name (detection time variable name)
# other required variables, used elsewhere: smp_code
parameters2$protocol_list$sampling_codes <- list(
  list(smp_fn = NA, smp_spec = NA, smp_describe = "use arrival times"),
  list(smp_fn = "BHG", smp_spec = parameters2$protocol_list$value, smp_describe = "use BHG function, prev ~ p(S->A), cases ~ Nmonth_mean")
) %>% bind_rows() %>% 
  mutate(smp_code = ifelse(is.na(smp_fn), NA, paste0(smp_fn, str_pad(smp_spec, width = 2, side = "left", pad = "0"))),
         smp_name = ifelse(is.na(smp_fn), "atimes", paste0("dTimes_", smp_code)))
print(paste(parameters2$protocol_list$txt, "refers to simultaneous application of these sampling protocols:"))
parameters2$protocol_list$sampling_codes



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  DETECTION FUNCTIONS 

## Daily detection function for BHG protocol
# v_MONTHspecimens = specimens per month [vector: each fID] 
# v_MONTHcases = sampling denominator (all monthly cases) [vector: each fID]
# v_prev = novel strain proportion of all cases [vector: each fID]
# function returns : probabilty of detection that day [vector: each fID]
fn_BHG <- function(v_MONTHspecimens, v_MONTHcases, v_prev) {
  v_MONTHnovelcases = ceiling(v_MONTHcases*v_prev)
  v_MONTHmissNovel = dhyper(x = 0,
                            m = v_MONTHnovelcases,
                            n = v_MONTHcases - v_MONTHnovelcases,
                            k = pmin(v_MONTHspecimens, v_MONTHcases))
  v_DAYfindnovel = 1-(v_MONTHmissNovel^(1/30))
  return(v_DAYfindnovel)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  SIMULATION FUNCTIONS 

## Single spread simulation
# calls fn_BHG 
# cmat = weighted adjacency matrix
# startHosp = seed location [hnrow value]
# runNum = this simulation (one of repeated simulations with this seed)
# stopnumber = stops an incomplete run at this timestep
# ptx = tranmission parameter
# sampling_codes = protocol parameter table [dataframe with smp_fn, smp_code, smp_spec, smp_name]
# cases = sampling denominator (all monthly cases) [vector/each fID]
# function returns: dataframe with source/runNum/results(atimes, dTimes)/max_timestep [one row per target location]
singleRun <- function(cmat, startHosp, runNum, 
                      stopnumber = max_timestep,
                      ptx = parameters2$transmission_parameter$value, 
                      sampling_codes = parameters2$protocol_list$sampling_codes, 
                      cases = Nmonth_mean){

  nHosp <- nrow(cmat)
  
  # set-up vectors of facility status: affected hospital <=> onward transmission source
  affectedHosps <- rep(0, nHosp)   # starts as 0, once infected -> 1
  affectedTimes <- rep(-1, nHosp)  # this vector is referenced for onward transmission
  detectedHosps_alt <- lapply(1:(nrow(sampling_codes[!is.na(sampling_codes$smp_fn),])), function(i) {rep(0, nHosp)})
  names(detectedHosps_alt) <- sampling_codes$smp_name[!is.na(sampling_codes$smp_fn)] %>% str_replace("dTimes", "dHosps")
  detectedTimes_alt <- lapply(1:(nrow(sampling_codes[!is.na(sampling_codes$smp_fn),])), function(i) {rep(-1, nHosp)})
  names(detectedTimes_alt) <- sampling_codes$smp_name[!is.na(sampling_codes$smp_fn)]
  
  # first timestep: introduce seed
  it <- 0
  
  # update affected status (startHosp only)
  affectedHosps[startHosp] <- 1
  affectedTimes[startHosp] <- it 
  
  # update detection status (startHosp only)
  # cycle through protocols (with check no unsupported detection functions)
  probChange <- 1 - exp(-ptx*cmat[startHosp, startHosp])
  randomNums_detect <- runif(1)
  for (i_smp in sampling_codes$smp_code[!is.na(sampling_codes$smp_fn)]) {
    i_update <- FALSE
    thisSMP_spec = sampling_codes$smp_spec[sampling_codes$smp_code %in% i_smp]
    # checking no supported functions
    if (sampling_codes$smp_fn[sampling_codes$smp_code %in% i_smp] == "BHG") {
      thisSMP_prob <- probChange
      probDetection <- fn_BHG(v_MONTHspecimens = thisSMP_spec, v_MONTHcases = cases[startHosp], v_prev = thisSMP_prob)
      i_update <- TRUE}
    # end of options
    if (i_update == FALSE) {stop(paste0("STOP", i_smp, "has no detection function"))}
    if (randomNums_detect <= probDetection) {detectedTimes_alt[[paste("dTimes", i_smp, sep = "_")]][startHosp] <- it}
    if (randomNums_detect <= probDetection) {detectedHosps_alt[[paste("dHosps", i_smp, sep = "_")]][startHosp] <- 1}
    }

  # iterate timestep while exist unaffected/undetected facilities (and have not reached the timestep limit)
  while ((it < stopnumber)&(min(c(sum(affectedHosps), sapply(detectedHosps_alt, sum))) < nHosp)){
    it<-it+1

    exportedPats <- base::colSums(cmat*affectedHosps)
    probChange<- 1-exp(-ptx*exportedPats)
    
    # update affected status (unaffected locations only)
    runTheseHospitals <- (affectedTimes == -1)
    randomNums <- rep(1, nHosp)
    randomNums[runTheseHospitals] <- runif(sum(runTheseHospitals))
    affectedTimes[(runTheseHospitals)&(randomNums<probChange)] <- it
    affectedHosps[randomNums<probChange] <- 1    

    # update detection (affected locations only)
    # cycle through protocols
    randomNums_detect <- rep(1, nHosp)
    randomNums_detect <- runif(nHosp)
    for (i in (names(detectedTimes_alt))) { 
    testTheseHospitals <- ((affectedHosps == 1) & (detectedTimes_alt[[i]] == -1))
    thisSMP_spec <- sampling_codes$smp_spec[sampling_codes$smp_name %in% i]
    thisSMP_prob <- probChange
    probDetection <- fn_BHG(v_MONTHspecimens = thisSMP_spec, v_MONTHcases = cases, v_prev = thisSMP_prob)
    detectedTimes_alt[[i]][(testTheseHospitals) & (randomNums_detect <= probDetection)] <- it
    detectedHosps_alt[[str_replace(i, "dTimes", "dHosps")]][(testTheseHospitals) & (randomNums_detect <= probDetection)] <- 1
    } #end i
  } #end while

  return(bind_cols(list(runNum = runNum, source = startHosp,
                        target = 1:length(affectedTimes),
                        atimes = affectedTimes,
                        dt = detectedTimes_alt)))
  } #end singleRun


## RUN SIMULATIONS FOR ALL SEED LOCATIONS  
# calls singleRun using default values for : stopnumber, ptx, sampling_codes, cases 
# contactMatrix = weighted adjacency matrix
# numRuns = number of simulations per seed location ('source')
# startHosps = NULL (cycle through all possible seed locations) or a single location [hnrow value]
# function returns: data frame of results [line per source/target pair per run]
runEntireSim <- function(contactMatrix, numRuns = sims_per_seed, startHosps = NULL){
  # set index case facility
  if (is.null(startHosps)) {
    starts <- (1:nrow(contactMatrix)) 
    } else {
      if (!all(is.numeric(startHosps),
               !length(startHosps) > 1,
               (startHosps %in% 1:nrow(contactMatrix)))) stop("startHosps should be NULL or single hnrow number")
      starts <- startHosps
  }
  # iterate through seed locations and run simulations
  affectedTimes <- bind_rows(
    lapply(starts, function(h) {
      print(paste("Running seed location:", h))
      bind_rows(lapply(1:numRuns, 
                       function(x) singleRun(cmat = contactMatrix, startHosp = h, runNum = x))) })
    ) %>% mutate(truncated = ifelse(if_any(contains("imes"), ~. <0), 
                 max_timestep,
                 NA))
  return(affectedTimes) 
} #end runEntireSim


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  RUN THE SIMULATIONS 

# reminder of values
max_timestep
sims_per_seed

# For timing purposes, run the simulation on a single starting location 
system.time(
  testResults <- runEntireSim(contactMatrix = hnmatrix, startHosps = 2, 
                              numRuns = 100)
  )
summary(testResults)
testResults %>% group_by(source, runNum) %>% summarise(lastAffectedTime = max(atimes)) %>% summary()
print(paste("Truncations before facilities affected:", 
            nrow(testResults %>% filter(atimes<0)), "/", nrow(testResults)))
print(paste("Truncations before detected using", 
            str_remove(colnames(testResults)[str_starts(colnames(testResults), "dTimes")], "dTimes_"), ":",
            sapply(colnames(testResults)[str_starts(colnames(testResults), "dTimes")],
                   function(i) {length(which(testResults %>% pull(i) < 0))}),
            "/", nrow(testResults)))


## Run all starting hospitals
tic()
system.time(
  sampleAffectedTimes <- runEntireSim(contactMatrix = hnmatrix, numRuns = sims_per_seed, startHosps = NULL)
)
toc() 

## onscreen summary
summary(sampleAffectedTimes)
sampleAffectedTimes %>% group_by(source, runNum) %>% summarise(lastAffectedTime = max(atimes)) %>% summary()
ggplot(sampleAffectedTimes %>% group_by(source, runNum) %>% summarise(lastAffectedTime = max(atimes))) + 
         geom_boxplot(aes(y = lastAffectedTime)) + ylim(c(0, max_timestep)) + labs(title = file_id2) + coord_flip()
paste("Truncations before affected:", 
      nrow(sampleAffectedTimes %>% filter(atimes<0)), "/", nrow(sampleAffectedTimes))
paste("Truncations before detection using", 
      str_remove(colnames(sampleAffectedTimes)[str_starts(colnames(sampleAffectedTimes), "dTimes")], "dTimes_"), ":",
      sapply(colnames(sampleAffectedTimes)[str_starts(colnames(sampleAffectedTimes), "dTimes")],
             function(i) {length(which(sampleAffectedTimes %>% pull(i) < 0))}),
      "/", nrow(sampleAffectedTimes))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> SAVE in working directory using file-naming convention
 
write_csv(sampleAffectedTimes,
          file = paste0(paste("sampleAffectedTimes", file_id2, sep = "_"), ".csv"))
write_csv(parameters2$protocol_list$sampling_codes,
          file = paste0(paste("samplingCodes", file_id2, sep = "_"), ".csv"))
read_csv(paste0(paste("sampleAffectedTimes", file_id2, sep = "_"), ".csv"))
read_csv(paste0(paste("samplingCodes", file_id2, sep = "_"), ".csv"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> OPTIONAL summaries and visualisations
#> diagnosing extent of any truncations

## truncations
sATsimtruncated <- sampleAffectedTimes %>% mutate(across(contains("imes"), ~.x < 0, .names = "TRUNC{.col}")) %>%
  mutate(sim = paste(source, runNum, sep = "_")) %>% group_by(sim) %>%
  summarise(source = unique(source), n_targets = n(), across(starts_with("TRUNC", ignore.case = FALSE), sum)) %>%
  mutate(across(starts_with("TRUNC"), ~.x>0, .names = "any{.col}"))
# fID truncated - counts per simulation
sATsimtruncated %>% select(sim, starts_with("TRUNC", ignore.case = FALSE)) %>%
  pivot_longer(!sim, names_to = "measure", values_to = "nTruncations") %>%
  tabyl(measure, nTruncations) %>% adorn_title()
# sims with any truncations - before affected, before detected
sATsimtruncated %>% select(sim, starts_with("anyTRUNC")) %>%
  pivot_longer(!sim, names_to = "measure", values_to = "hasTruncations") %>%
  group_by(measure) %>% summarise(n_sims = n(), hasTruncations = sum(hasTruncations))

