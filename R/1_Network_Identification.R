##################################################
####                                          ####
#### 1) READ IN DATA / NETWORK IDENTIFICATION ####
####                                          ####
##################################################

#### reads in data (including HospitalNetwork object containing network adjacency matrix) 
#### checks some formatting
#### applies carriage weighting to network movements
#### prepares lookup files used in later stages
#### creates and saves hnmatrix[rds], hnreallife[csv], hndenominatorcases[csv], hncasesmodel[csv], fIDsubsets[rds]

library(tidyverse)
library(readr)
library(janitor)
library(igraph)

setwd()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> SET PARAMETERS for step_1 : Network Identification
#> also creates file prefix for output
#> *** eyeball check consistency of patient/year/window vs data read-in below ***

parameters1 <- list(
  "patients_definition" = list("txt" = "DUMMY"),                # text used in file_id
  "financial_year" = list("fy_enddate" = ymd("2020-03-31")),    # file_id uses year of data period end    
  "window_length" = list("value" = 182),          
  "carriage_weighting" = list("value" = 0.03)                   # (e.g. 0.03 for 3% carriage) file_id takes nearest %age
  )

# additional information for specimen burden model
draws_from_distribution <- 100

# parameter checking and file-name suffix generation 
parameters1$financial_year$txt <- paste0("fy", year(parameters1$financial_year$fy_enddate))
parameters1$window_length$txt <- paste0("wdw", parameters1$window_length$value)
parameters1$carriage_weighting$txt <- paste0("cw", str_pad(round_half_up(parameters1$carriage_weighting$value*100), width = 2, side = "left", pad =  "0"))
file_id1 <- paste(sapply(parameters1, "[[", "txt"), collapse = "-")
file_id1

## crib sheet for spells to be included in the network per financial year
fn_spells_cribsheet <- function(financialyearend, window) {
  fy_start <- ymd(paste(financialyearend-1, "April", "01", sep = "-"))
  fy_end <- ymd(paste(financialyearend, "March", "31", sep = "-"))
  admission_interval <- interval(fy_start, fy_end)
  discharge_interval <- interval(fy_start - days(window), fy_end)
  cat(paste("REMINDER: Parameters refer to HospitalNetwork period = ", paste0("fy", year(fy_end)), 
            " with inter-spell window", window, "days", "\n", 
            "so uses", "\n", "spells which have admission",
            format(int_start(admission_interval), "%d/%b/%Y"), "-", format(int_end(admission_interval), "%d/%b/%Y"), "\n",
            "and spells which have discharge", 
            format(int_start(discharge_interval), "%d/%b/%Y"), "-", format(int_end(discharge_interval), "%d/%b/%Y")
            ))
  }
fn_spells_cribsheet(year(parameters1$financial_year$fy_enddate), parameters1$window_length$value)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  DATA TO READ IN
#>  ensure all fID codes refer to same facility throughout
#>  *** EYEBALL CHECK HospiNet matches patient/year/window in parameters1 ***

## HospiNet object
  # pure patient movements i.e. unadjusted for carriage
  # capture its specifications in hn_details for use in file-naming convention
  # extract just the 'matrix' : adjacency matrix (row/column names are fID)
## facility information 
  # csv file including columns: fID/fID_fullname/fID_type/fID_region
  # information as correct for "fIDyear"
  # facility codes ("fID") as used in HospiNet object generation
## annual cases (protocol denominators e.g. multi-year average)
  # csv file including columns: fID/fID_meanannualtotalcases
  # facility codes ("fID") as used in HospiNet object generation for "fIDyear"
  # csv file including columns: mean, sd
## fitted normal (mean = 1) for annual cases distribution vs denominator values
  # csv file (single row) including columns: mean, sd
  # parameters from a normal which has been fitted across facilities used in HospiNet object generation

hn_file <- list("fIDyear" = 2020,
                "hn_file_location" = getwd(),
                "hn_file_name" = "SMALLDUMMY_patients182days")
fIDinfo_file <- list("fIDyear" = 2020,
                     "fIDinfo_file_location" = getwd(),
                     "fIDinfo_file_name" = "SMALLDUMMY_fIDinfo")
fIDannualcases_file <- list("fIDyear" = 2020,
                            "fIDannualcases_file_location" = getwd(),
                            "fIDannualcases_file_name" = "SMALLDUMMY_annualcases")
fIDfittednormal_file <- list("fIDyear" = 2020,
                             "fIDfittednormal_file_location" = getwd(),
                             "fIDfittednormal_file_name" = "SMALLDUMMY_fittednormal")

if (n_distinct(c(hn_file$hn_details$fIDyear,fIDinfo_file$fIDyear, fIDannualcases_file$fIDyear,  fIDfittednormal_file$fIDyear)) != 1) {stop("STOP check files relating to fID from the HospiNet dataperiod")}
stopifnot(file.exists(paste0(hn_file$hn_file_location, "/", hn_file$hn_file_name, ".rds"),
                      paste0(fIDinfo_file$fIDinfo_file_location, "/", fIDinfo_file$fIDinfo_file_name, ".csv"),
                      paste0(fIDannualcases_file$fIDannualcases_file_location, "/", fIDannualcases_file$fIDannualcases_file_name, ".csv"),
                      paste0(fIDfittednormal_file$fIDfittednormal_file_location, "/", fIDfittednormal_file$fIDfittednormal_file_name, ".csv")))

# read in data
hn <- read_rds(paste0(hn_file$hn_file_location, "/", hn_file$hn_file_name, ".rds"))$matrix
fIDinfo <- read_csv(paste0(fIDinfo_file$fIDinfo_file_location, "/", fIDinfo_file$fIDinfo_file_name, ".csv"))
fIDannualcases <- read_csv(paste0(fIDannualcases_file$fIDannualcases_file_location, "/", fIDannualcases_file$fIDannualcases_file_name, ".csv"))
fIDfittednormal <- read_csv(paste0(fIDfittednormal_file$fIDfittednormal_file_location, "/", fIDfittednormal_file$fIDfittednormal_file_name, ".csv"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  SOME FORMAT CHECKING/CLEANING AND APPLY CARRIAGE WEIGHTING

## hn (adjacency matrix)
# check for NA in adjacency matrix and replace with zero
# apply carriage weighting [here rather than in flexible sentinel analysis]
hn_tidy <- hn
if (sum(is.na(hn)) > 0) {
  print(paste("Replacing", sum(is.na(hn)), "NA values")) 
  hn_tidy[is.na(hn)] <- 0}
hn_tidy <- hn_tidy * parameters1$carriage_weighting$value
dim(hn_tidy)

## fIDinfo 
# check for column names, unique entry by fID, includes all hn's fID
# tidy to match hn$fID
# add columns to reference hn
tmp_requiredcolumns <- c("fID", "fID_fullname", "fID_type", "fID_region")
if (!all(all(tmp_requiredcolumns %in% names(fIDinfo)),
         all.equal(unique(fIDinfo$fID), fIDinfo$fID),
         length(setdiff(rownames(hn), fIDinfo$fID)) == 0)) {stop("STOP fIDinfo incorrect format / missing fID vs hn")}
if (length(setdiff(fIDinfo$fID, rownames(hn))) > 0) {
  print(paste("Removing", length(setdiff(fIDinfo$fID, rownames(hn))), "fIDinfo fID not in hn"))}
if (length(setdiff(names(fIDinfo), tmp_requiredcolumns)) > 0) {
  print(paste("Removing unwanted columns:", setdiff(names(fIDinfo), tmp_requiredcolumns)))}
fIDinfo_tidy <- left_join(data.frame(fID = rownames(hn)) %>% 
                            mutate(hnrow = row_number(), hnfile = file_id1, hnyear =parameters1$financial_year$txt),
                          fIDinfo %>% select(any_of(tmp_requiredcolumns)))
dim(fIDinfo_tidy)
glimpse(fIDinfo_tidy)

## fIDannualcases 
# check for column names, unique entry by fID, includes all hn's fID, no NA values
# tidy to match hn$fID
# add columns to reference hn
tmp_requiredcolumns <- c("fID", "fID_meanannualtotalcases")
if (!all(all(tmp_requiredcolumns %in% names(fIDannualcases)),
        all.equal(unique(fIDannualcases$fID), fIDannualcases$fID),
        length(setdiff(rownames(hn), fIDannualcases$fID)) == 0,
        all(!is.na(fIDannualcases$fID_meanannualtotalcases)))) {stop("STOP fIDannualcases incorrect format / missing fID vs hn")}
if (length(setdiff(fIDannualcases$fID, rownames(hn))) > 0) {
  print(paste("Removing", length(setdiff(fIDannualcases$fID, rownames(hn))), "fIDannualcases fID not in hn"))}
if (length(setdiff(names(fIDannualcases), tmp_requiredcolumns)) > 0) {
  print(paste("Removing unwanted columns:", setdiff(names(fIDannualcases), tmp_requiredcolumns)))}
fIDannualcases_tidy <- left_join(data.frame(fID = rownames(hn)) %>%
                                   mutate(hnrow = row_number(), hnfile = file_id1, hnyear =parameters1$financial_year$txt), 
                                 fIDannualcases %>% select(any_of(tmp_requiredcolumns)))
dim(fIDannualcases_tidy)
glimpse(fIDannualcases_tidy)

## fIDfittednormal
# check column names, unique value, mean = 1
tmp_requiredcolumns <- c("mean", "sd")
if (!all(all(tmp_requiredcolumns %in% names(fIDfittednormal)),
         nrow(fIDfittednormal) == 1,
         fIDfittednormal$mean == 1)) {stop("STOP fIDfittednormal incorrect format /  mean is not 1")}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> PREPARE LOOKUP for use in defining candidate sentinels
#> Can be appended with other definitions in later steps e.g. Benchmark sentinel set
 
## Extract setsets of Specialist and Teaching facilities
# Txt are used in file-naming convention as shorthand for the fID subset
# fIDset is vector of fID from hnreallife
fIDsubsets <- list(
  "Specialist" = list(txt = "Spl",
                      fIDvector = fIDinfo_tidy %>% filter(fID_type == "ACUTE - SPECIALIST") %>% pull(fID)),
  "Teaching" = list(txt = "Tch",
                    fIDvector = fIDinfo_tidy %>% filter(fID_type == "ACUTE - TEACHING" ) %>% pull(fID))
  )


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> PREPARE LOOKUP for use in specimen burden estimation
#> draws from the fitted normal

## draw from fitted normal (mean = 1)
# row per FID each with draws_from_distribution random values
hncasesdraws <- 
  full_join(by = "fID",
            data.frame(fID = rownames(hn)) %>% 
              mutate(hnrow = row_number(), "fit_mean" = fIDfittednormal$mean, "fit_sd" = fIDfittednormal$sd),
            lapply(rownames(hn), function(i) {
              rnd_figure = rnorm(n = draws_from_distribution, mean = 1, sd = fIDfittednormal$sd)
              bind_cols("fID" = i, "draw_fitted" = rnd_figure) %>% mutate(id = paste("draw", row_number(), sep = "_")) 
              }) %>% bind_rows())


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>  IS THE NETWORK A SINGLE STRONGLY-CONNECTED COMPONENT ?
#>  If true, then a strain emerging at any facility can spread to all facilities: procede with analysis
#>  If untrue, consider piecewise analysis (code not given here) 

hngraph <- igraph::graph_from_adjacency_matrix(adjmatrix = hn_tidy,
                                               mode = "directed", weighted = TRUE)
if (igraph::count_components(hngraph, mode = "strong") != 1) {
  stop(paste0("**** STOP ****  this network requires alternative treatment as has",
              igraph::count_components(hngraph, mode = "strong"), "strongly connected components"))}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> SAVE in working directory using file-naming convention
 

write_rds(hn_tidy,
          file = paste0(paste("hnmatrix", file_id1, sep = "_"), ".rds"))
write_csv(fIDinfo_tidy,
          file = paste0(paste("hnreallife", file_id1, sep = "_"), ".csv"))
write_csv(fIDannualcases_tidy,
          file = paste0(paste("hndenominatorcases", file_id1, sep = "_"), ".csv"))
write_csv(hncasesdraws,
          file = paste0(paste("hncasesdraws", file_id1, sep = "_"), ".csv"))
write_rds(fIDsubsets,
          file = paste0(paste("fIDsubsets", file_id1, sep = "_"), ".rds"))

read_rds(paste0(paste("hnmatrix", file_id1, sep = "_"), ".rds"))
read_csv(paste0(paste("hnreallife", file_id1, sep = "_"), ".csv"))
read_csv(paste0(paste("hndenominatorcases", file_id1, sep = "_"), ".csv"))
read_csv(paste0(paste("hncasesdraws", file_id1, sep = "_"), ".csv"))
read_rds(paste0(paste("fIDsubsets", file_id1, sep = "_"), ".rds"))




