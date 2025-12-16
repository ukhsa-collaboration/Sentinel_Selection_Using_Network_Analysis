####################################
####                            ####
#### 4) SENTINEL SET EVALUATION ####
####                            ####
####################################

#### read in sentinelcsv (priority list)
#### n-member sentinel set is the top n from the priority list

#### NB detection performance estimation uses priority list's "native scenario/protocol" i.e. from inherited network/spread data 
## To evaluate performance of a priority list or n-member set in a non-native parameter context:
# if necessary, revisit step2 (or step1 and step2) to provide spread data for the parameter context
# append evaluation set to fIDsubset for the inherited network of the non-native context (code in 4appendix)
# run step3 using the "force" criteria with that appended subset
# => the subset (in order) is top of resultant priority list which inherits the required context, which can be evaluated here

#### evaluations available (for single priority list)
# 1) create and save sentinelregionPlot[png], sentineltypePlot[png] : region and type profile of n-member sentinel sets
# 2) create and save sentinelsetsTimes[csv] and sentinelsetsTimesPlot[png] : mean detection time estimates of n-member sentinel sets
# 3) create and save sentinelsetsSpecimens[csv] and sentinelsetsSpecimentsPlot[png] : sentinel set sampling burden estimates (from fitted normal)
# 4) create and save sentinelsetsAffected[csv] and sentinelsetsAffectedPlot[png] : number of facilities affected at detection 
# 5) create and save randomsetsTimes[csv] and randomsetsTimesPlot[png] : random sets mean detection time estimates


library(tidyverse)
library(readr)
library(janitor)
library(tictoc)

setwd()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> SPECIFY PRIORITY LIST FOR EVALUATION
#> suffix for the sentinelcsv file to evaluate (after "sentinelcsv_")
#> check exists 

file_id3 <- "DUMMY-fy2020-wdw182-cw03_p00001-sp10-i100_meanBHG10-excSpl-noCap-noFrc"
stopifnot(file.exists(paste0(paste("sentinelcsv", file_id3, sep = "_"), ".csv")))

## nominate a sample set size to view ON SCREEN specific details
# saved files have details for all set sizes
nMember <- 4

## for evaluation (5)
# number of random sets to generate
nRandom <- 100


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> READ IN SENTINELCSV
#> check all other intermediate files exist
#> set shared values used in plotting

priorityList <- read_csv(paste0(paste("sentinelcsv", file_id3, sep = "_"), ".csv"))

# extract corresponding intermediate stage id and check files exist
file_id1 <- str_split_i(file_id3, pattern = "_", i = 1)
file_id2 <- paste(file_id1, str_split_i(file_id3, pattern = "_", i = 2), sep = "_")
stopifnot (sapply(FUN = file.exists, str_replace(c("hnmatrix_NETWORKid.rds", "hnreallife_NETWORKid.csv", "hndenominatorcases_NETWORKid.csv", "hncasesdraws_NETWORKid.csv", "fIDsubsets_NETWORKid.rds"), 
                                                 pattern = "NETWORKid", replacement = file_id1) ))
stopifnot (sapply(FUN = file.exists, str_replace(c("sampleAffectedTimes_SPREADid.csv", "samplingCodes_SPREADid.csv"), 
                                                 pattern = "SPREADid", replacement =  file_id2)))
stopifnot (sapply(FUN = file.exists, str_replace(c("criteriaCribsheet_SENTINELid.rds"), 
                                                 pattern = "SENTINELid", replacement =  file_id3)))

# check nMember vs priority list length
stopifnot (nMember <= nrow(priorityList))

# shared plot information (scaled to title width)
plotwidthdefault = 1 + nchar(file_id3)*0.22
plotheightdefault = 10


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> DETECTION TIME ESTIMATION functions

# estimates detection time using a single set of simulations (from all seed locations)
# oT : oT is dataframe of set of simulations (unique source/target pairs) {runNum, source, target, obsT}
# iSL: priority list in order (hnrow values for facilities)
# function output: vector of times per n-member sentinel sets (n is 1 to whole priority list)
fn_TSingle <- function(oT, iSL) {
  sapply(lapply(1:length(iSL), function(n){iSL[1:n]}),
         function(nSentinelSet) { 
           mean((oT[oT$target %in% nSentinelSet, ] %>% 
                   group_by(source) %>% 
                   summarise(firstDetect = min(obsT), .groups = "drop"))$firstDetect)
         })} #end fn_TSingle

## estimates detection time for sentinel sets taken from ordered vector of facilities [hnvalue as used by sAT] 
# calls fn_TSingle
# sAT : spread simulation output file starting "sampleAffectedTimes" 
# obs : use this protocol's detection time, use smp_code (NA or 'atimes' for 'sequence all cases')
# facility_indices : vector of facilities [hnrow values] from same hnmatrix inherited by sAT
# function output : longformat dataframe with times for first detection by nSentinels-size sentinel set {Sentinels, runN, runDT}
fn_TEstimateAll <- function(sAT, obs, facility_indices) {
  # check sAT has obs and facilities_indices
  if (is.na(obs)) {obs <- "atimes"}
  if (!all(facility_indices %in% sAT$source)) {stop(paste("Not all facility indices are targets in sAT"))}
  if (!any(str_detect(names(sAT), pattern = obs))) {stop(paste(obs, "not in sAT:", paste(names(sAT), collapse = ", ")))}
  obstimes = sAT %>% select(runNum, source, target, contains(obs)) %>% rename("obsT" = last_col())
  sapply(1:max(obstimes$runNum), 
         function(r) {fn_TSingle(oT = obstimes[obstimes$runNum == r, ], 
                                   iSL = facility_indices)}) %>%
    as.data.frame(make.names = TRUE) %>% mutate(nSentinels = row_number()) %>%
    pivot_longer(!nSentinels, names_to = "runN", values_to = "runDT")
} #end fn_TEstimateAll 


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> 1) REGION AND TYPE PROFILE                                               ####

## read in inherited facility information
hnreallife <- read_csv(paste0(paste("hnreallife", file_id1, sep = "_"), ".csv"))

## calculate and plot numbers by strata
# for all n-member sentinel sets (top n on prioritylist)
profilebyregion <- bind_cols(priorityList %>% select(sentinelrank),
                              sapply(unique(hnreallife$fID_region), function(i) {cumsum(priorityList$fID_region == i)}))
profilebytype <- bind_cols(priorityList %>% select(sentinelrank),
                            sapply(unique(hnreallife$fID_type), function(i) {cumsum(priorityList$fID_type == i)}))
profilebyregionPlot <- ggplot(data = profilebyregion %>% pivot_longer(!sentinelrank, names_to = "region", values_to = "n")) +
  geom_line(aes(x = sentinelrank, y = n, colour = region)) + labs(title = file_id3, x = "sentinel set size", y = "facilties included") +
  theme(plot.title = element_text(size = 11))
profilebyregionPlot
profilebytypePlot <- ggplot(data = profilebytype %>% pivot_longer(!sentinelrank, names_to = "type", values_to = "n")) +
  geom_line(aes(x = sentinelrank, y = n, colour = type)) + labs(title = file_id3, x = "sentinel set size", y = "facilties included") +
  theme(plot.title = element_text(size = 11))
profilebytypePlot

## VIEW summary for nominated set size
glimpse(profilebyregion %>% filter(sentinelrank == nMember))  
glimpse(profilebytype %>% filter(sentinelrank == nMember))  

## save plots in working directory using file-naming convention
ggsave(profilebyregionPlot + labs(subtitle = "profilebyregionPlot"),
       filename = paste0(paste("profilebyregionPlot", file_id3, sep = "_"), ".png"),
       units = "cm", width = plotwidthdefault, height = plotheightdefault)
ggsave(profilebytypePlot + labs(subtitle = "profilebytypePlot"),
       filename = paste0(paste("profilebytypePlot", file_id3, sep = "_"), ".png"),
       units = "cm", width = plotwidthdefault, height = plotheightdefault)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> 2) MEAN DETECTION TIME ESTIMATES sentinel sets                           ####

## read inherited spread data and protocol
sampleAffectedTimes <- read_csv(paste0(paste("sampleAffectedTimes", file_id2, sep = "_"), ".csv"), 
                                col_types = cols(truncated = col_number()))
obs <- str_split_i(str_split_i(file_id3, pattern = "_", i = 3), pattern = "-", i = 1) %>% str_remove("mean")

## estimate detection times for all n-member sentinel sets (top n on prioritylist)
# summary stats for all n-member sentinel sets (top n on prioritylist)

sentineltimes <- fn_TEstimateAll(sAT = sampleAffectedTimes, obs = obs, facility_indices = priorityList$hnrow) %>%
  group_by(nSentinels) %>%
  summarise(T_mn = mean(runDT), T_sd = sd(runDT)) %>%
  mutate(nsim = nrow(distinct(sampleAffectedTimes, runNum, source))) %>% 
  mutate(file_id3 = file_id3)

sentineltimesPlot <- ggplot(data = sentineltimes) +
  geom_point(aes(x = nSentinels, y = T_mn), size = 1) +
  geom_errorbar(aes(x = nSentinels, ymin = T_mn - T_sd, ymax = T_mn + T_sd)) +
  labs(title = file_id3, x = "sentinel set size", y = "mean detection time (+/- sd)") +
  theme(plot.title = element_text(size = 11))
sentineltimesPlot

## VIEW summary for nominated set size
glimpse(sentineltimes %>% filter(nSentinels == nMember))

## save summary dataframe and plot in working directory using file-naming convention
write_csv(sentineltimes,
          file = paste0(paste("sentinelsetsTimes", file_id3, sep = "_"), ".csv"))
ggsave(sentineltimesPlot + labs(subtitle = "sentinelsetsTimesPlot"),
       filename = paste0(paste("sentinelsetsTimesPlot", file_id3, sep = "_"), ".png"),
       units = "cm", width = plotwidthdefault, height = plotheightdefault)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> 3) SAMPLING BURDEN ESTIMATES sentinel sets                               ####

## inherited parameters
incidence_multiplier <- str_split_i(str_split_i(file_id3, pattern = "_", i = 2), pattern = "-", i = 3) %>% str_remove("i") %>% as.numeric() * 0.01
obs <- str_split_i(str_split_i(file_id3, pattern = "_", i = 3), pattern = "-", i = 1) %>% str_remove("mean")

hndenominatorcases <- read_csv(paste0(paste("hndenominatorcases", file_id1, sep = "_"), ".csv"))
samplingCodes <- read_csv(paste0(paste("samplingCodes", file_id2, sep = "_"), ".csv"))
hncasesdraws <- read_csv(paste0(paste("hncasesdraws", file_id1, sep = "_"), ".csv"))

## edit hncasesdraws 
# add incidence-adjusted mean annual cases and monthly specimens
hncasesdraws_iSpec <- left_join(
  priorityList %>% select(sentinelrank, fID) %>% mutate(file_id3 = file_id3),
  left_join(
    hndenominatorcases %>% mutate(Nannual_cases = fID_meanannualtotalcases*incidence_multiplier) %>% select(fID, hnrow, Nannual_cases) %>%
      mutate(Nmonth_specimens = ifelse(obs == "atimes", Inf, samplingCodes$smp_spec[which(samplingCodes$smp_code %in% obs)])),
    hncasesdraws) %>%
    mutate(draw_AnnualSpec = pmin(draw_fitted*Nannual_cases, 12*Nmonth_specimens))
  ) %>% select(-draw_fitted) %>% pivot_wider(names_from = id, values_from = draw_AnnualSpec)

## estimate annual specimen burden for all n-member sentinel sets (top n on prioritylist)
# summary stats for all n-member sentinel sets (top n on prioritylist)
sentinelspecimens <- lapply(1:nrow(hncasesdraws_iSpec), 
                            function(i) {colSums(hncasesdraws_iSpec %>% filter(sentinelrank <= i) %>% select(starts_with("draw_")))}) %>%
  bind_rows() %>%
  mutate(nSentinels = row_number(), 
         ndraw = max(as.numeric(str_remove(hncasesdraws$id, "draw_")))) %>% rowwise() %>%
  mutate(file_id3 = file_id3, permonth = unique(hncasesdraws_iSpec$Nmonth_specimens),
         S_mn = mean(c_across(starts_with("draw_"))), S_sd = sd(c_across(starts_with("draw_")))) %>% ungroup() %>%
  relocate(starts_with("draw_"), .after = "S_sd") 
sentinelspecimensPlot <-
  ggplot(data = sentinelspecimens %>% select(nSentinels, starts_with("draw_")) %>% 
           pivot_longer(!nSentinels, names_to = "draw", values_to = "specimens")) + 
  geom_violin(aes(x = nSentinels, y = specimens, group = nSentinels)) +
  labs(title = file_id3, x = "sentinel set size", y = paste("annual specimens", "(", unique(sentinelspecimens$ndraw), "random draws)")) +
  theme(plot.title = element_text(size = 11))
sentinelspecimensPlot 

## VIEW summary for nominated set size
summary(sentinelspecimens %>% filter(nSentinels == nMember) %>% select(nSentinels, starts_with("draw_")) %>%
          pivot_longer(!nSentinels,  names_to = "draw", values_to = "specimens") %>% select(-draw))
ggplot(data = sentinelspecimens %>% filter(nSentinels == nMember) %>% select(nSentinels, starts_with("draw_")) %>% 
         pivot_longer(!nSentinels, names_to = "draw", values_to = "specimens")) + 
  geom_histogram(aes(specimens, after_stat(density))) +
  ggtitle(paste0(nMember, "-member sentinel set, annual specimens"))

## save summary dataframe and plot in working directory using file-naming convention
write_csv(sentinelspecimens,
          file = paste0(paste("sentinelsetsSpecimens", file_id3, sep = "_"), ".csv"))
ggsave(sentinelspecimensPlot + labs(subtitle = "sentinelsetsSpecimensPlot"),
       filename = paste0(paste("sentinelsetsSpecimensPlot", file_id3, sep = "_"), ".png"),
       units = "cm", width = plotwidthdefault, height = plotheightdefault)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> 4) NUMBER OF AFFECTED FACILITIES AT DETECTION sentinel sets              ####

## read inherited spread data and protocol
sampleAffectedTimes <- read_csv(paste0(paste("sampleAffectedTimes", file_id2, sep = "_"), ".csv"), 
                                col_types = cols(truncated = col_number()))
obs <- str_split_i(str_split_i(file_id3, pattern = "_", i = 3), pattern = "-", i = 1) %>% str_remove("mean")

## estimate affected facilities at detection times for all n-member sentinel sets (top n on prioritylist)
# summary stats for all n-member sentinel sets (top n on prioritylist)
fn_AffectedSingle <- function(sAT = sampleAffectedTimes, pList = priorityList, dT = obs, nSet) {
  sAT_nMember <-
    left_join(
      left_join(sAT %>% select(runNum, source, target, atimes),
                sAT %>% select(runNum, source, target, contains(dT)) %>% rename("obsT" = last_col())),
      pList %>% filter(sentinelrank <= nSet) %>%
        mutate(target = hnrow, targetIsSentinel = TRUE) %>% select(target, targetIsSentinel),
      by = "target") %>%
    mutate(id = paste0("s", source, "r", runNum))
  sAT_nMember <- left_join(
    sAT_nMember,
    sAT_nMember %>% filter(targetIsSentinel) %>% group_by(id) %>% summarise(dTimes = min(obsT)) %>% ungroup())
  sAT_nMember$affected <- (sAT_nMember$atimes <= sAT_nMember$dTimes)
  sAT_nMember <- 
    sAT_nMember %>% group_by(id) %>% arrange(id, atimes) %>%
    filter(affected) %>% add_tally(name = "nAffected") %>%
    distinct(source, id, dTimes, nAffected) %>% ungroup() %>%
    mutate(nSentinels = nSet, nSims = n_distinct(id))
  return(sAT_nMember)
}
sentinelaffected <- lapply(1:nrow(priorityList), function(i) {
  fn_AffectedSingle(nSet = i) %>% summarise(nSentinels = unique(nSentinels), A_mn = mean(nAffected), A_sd = sd(nAffected))
  }) %>% bind_rows() %>% 
  mutate(nsim = nrow(distinct(sampleAffectedTimes, runNum, source))) %>% 
  mutate(file_id3 = file_id3)
sentinelaffectedPlot <- 
  ggplot(data = sentinelaffected) + 
  geom_point(aes(x = nSentinels, y = A_mn), size = 1) +
  geom_errorbar(aes(x = nSentinels, ymin = A_mn - A_sd, ymax = A_mn + A_sd)) +
  labs(title = file_id3, x = "sentinel set size", y = "mean affected facilties (+/- sd)") +
  theme(plot.title = element_text(size = 11))
sentinelaffectedPlot

## VIEW summary for nominated set size
fn_AffectedSingle(nSet = nMember) %>% tabyl(nAffected) %>% adorn_totals() 
ggplot(data = fn_AffectedSingle(nSet = nMember)) +
  geom_histogram(aes(nAffected, after_stat(density)), binwidth = 1) +
  ggtitle(paste0(nMember, "-member sentinel set, affected facilities"))

## save summary dataframe and plot in working directory using file-naming convention
write_csv(sentinelaffected,
          file = paste0(paste("sentinelsetsAffected", file_id3, sep = "_"), ".csv"))
ggsave(sentinelaffectedPlot + labs(subtitle = "sentinelsetsAffectedPlot"),
       filename = paste0(paste("sentinelsetsAffectedPlot", file_id3, sep = "_"), ".png"),
       units = "cm", width = plotwidthdefault, height = plotheightdefault)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> 5) MEAN DETECTION TIME ESTIMATES random sets                             ####
#> uses same restrictions as in priorityList

## reminder: number of random sets
nRandom

## read inherited spread data, protocol and criteria
sampleAffectedTimes <- read_csv(paste0(paste("sampleAffectedTimes", file_id2, sep = "_"), ".csv"), 
                                col_types = cols(truncated = col_number()))
obs <- str_split_i(str_split_i(file_id3, pattern = "_", i = 3), pattern = "-", i = 1) %>% str_remove("mean")
criteria <- read_rds(paste0(paste("criteriaCribsheet", file_id3, sep = "_"), ".rds"))

## generate random sets
random_indices <- lapply(1:nRandom,  function(x) {
  f = sample(unique(sampleAffectedTimes$target), replace = FALSE) 
  f = setdiff(f, criteria$set_exclude$hnrow)
  f = c(criteria$set_force$hnrow, setdiff(f, criteria$set_force$hnrow)) 
  if (length(criteria$set_cap$hnrow) > 0) {
    f = setdiff(f, 
                f[which(f %in% criteria$set_cap$hnrow)[(criteria$Ncap+1):n_distinct(criteria$set_cap$hnrow)]]
    )}
  return(f)}) 

# For timing purposes, run for a single random set 
system.time(
  testRandom <- fn_TEstimateAll(sAT = sampleAffectedTimes, obs = obs, facility_indices = random_indices[[1]])
  )

## estimate detection times for all n-member sentinel sets (top n on prioritylist)
# summary stats across all random sets for all n-member sentinel sets (top n on prioritylist)
randomtimes <- lapply(1:nRandom, function(x) {
  print(paste("running random set", x))
  fn_TEstimateAll(sAT = sampleAffectedTimes, obs = obs, facility_indices = random_indices[[x]])
}) %>% 
  bind_rows() %>%
  group_by(nSentinels) %>%
  summarise(T_mn = mean(runDT), T_sd = sd(runDT)) %>%
  mutate(nsim = nrow(distinct(sampleAffectedTimes, runNum, source)), rdmsets = nRandom) %>%
  mutate(file_id3 = file_id3)

randomtimesPlot <- ggplot(data = randomtimes) +
  geom_point(aes(x = nSentinels, y = T_mn), size = 1) +
  geom_ribbon(aes(x = nSentinels, ymin = T_mn - T_sd, ymax = T_mn + T_sd), alpha = 0.2) +
  labs(title = file_id3, x = paste("sentinel set size", paste0("(", unique(randomtimes$rdmsets), " random sets)")), y = "mean detection time (+/- sd)") +
  theme(plot.title = element_text(size = 11))
randomtimesPlot

## VIEW summary for nominated set size
glimpse(randomtimes %>% filter(nSentinels == nMember))

## save summary dataframe and plot in working directory using file-naming convention
write_csv(randomtimes,
          file = paste0(paste("randomsetsTimes", file_id3, sep = "_"), ".csv"))
ggsave(randomtimesPlot + labs(subtitle = "randomsetsTimesPlot"),
       filename = paste0(paste("randomsetsTimesPlot", file_id3, sep = "_"), ".png"),
       units = "cm", width = plotwidthdefault, height = plotheightdefault)

