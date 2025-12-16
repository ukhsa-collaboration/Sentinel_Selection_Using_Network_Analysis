###################################
####                           ####
#### 0) DUMMY DATA             ####
####                           ####
###################################

#### Provides example set of dummy data files

#### Uses fake patient database from HospitalNetwork package
# clusters version used, to provide a proxy for facility's region
# other heterogeneity in facilities and network connections (type, size etc) not specifiable
# absolute dates not specifiable
#### fn_spells_cribsheet (in 1_Network_Identification) is provided to inform date specification for a real world patient database

library(tidyverse)
library(HospitalNetwork)
library(janitor)
library(igraph)

setwd()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> LINE LIST OF STAYS
## each record has "sID" = subject ID, "fID" = facility ID, "Adate" =  admission date, "Ddate" = discharge date

## WARNING: large scale line list generation is time-consuming (order of hours for patient counts in the millions)
stay_linelist = HospitalNetwork::create_fake_subjectDB_clustered(
  n_subjects = 125000, 
  n_facilities = 27,
  avg_n_stays = 3,
  n_clusters = 3)
summary(stay_linelist)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> create HOSPINET object

# check requirements
stay_linelist <- HospitalNetwork::checkBase(stay_linelist)
# create HospiNet object
DUMMYpatients182days <- hospinet_from_subject_database(base = stay_linelist, window_threshold = 182, noloops = FALSE)
# check single strongly connected component STOP if not
if (igraph::count_components(DUMMYpatients182days$igraph, mode = c("strong")) != 1) {
  stop("REDO LINE LIST as require a single strongly connected component")}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> create random "REAL LIFE" information

DUMMY_fIDinfo <- data.frame("fID" = rownames(DUMMYpatients182days$matrix)) %>%
  mutate("fID_fullname" = paste("Hospital", fID)) %>% 
  left_join(as.data.frame(DUMMYpatients182days$cluster_fast_greedy) %>% 
              rename("fID" = node) %>% mutate(fID_region = paste("Region", LETTERS[cluster_fast_greedy])) %>%
              select(fID, fID_region))
DUMMY_fIDinfo$fID_type <- sapply(
  runif(nrow(DUMMYpatients182days$matrix), min = 0, max = 1),
  function(x) {case_when(x < 0.1 ~ "ACUTE - SPECIALIST",        # 10% specialist
                         x < (0.1 + 0.4) ~ "ACUTE - TEACHING",  # 40% teaching
                         TRUE ~ "ACUTE - OTHER")})              # 50% other
tabyl(DUMMY_fIDinfo, fID_type, fID_region) %>% adorn_totals(c("row", "col"))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> create random mean CASES per facility and nominate value for fitted normal sd

DUMMY_annualcases <-
  data.frame("fID" = rownames(DUMMYpatients182days$matrix))
DUMMY_annualcases$fID_meanannualtotalcases <-
  pmax(1, colSums(DUMMYpatients182days$matrix)*rnorm(n = nrow(DUMMYpatients182days$matrix), mean = 0.001, sd = 0.0005))
ggplot(DUMMY_annualcases, aes(fID_meanannualtotalcases)) + geom_histogram()

DUMMY_fittednormal <- data.frame("mean" = 1, "sd" = 0.143)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> save dummy data files


write_rds(DUMMYpatients182days, file = paste0("DUMMYpatients182days", ".rds"))
write_csv(DUMMY_fIDinfo, file = paste0("DUMMY_fIDinfo", ".csv"))
write_csv(DUMMY_annualcases, file = paste0("DUMMY_annualcases", ".csv"))
write_csv(DUMMY_fittednormal, file = paste0("DUMMY_fittednormal", ".csv"))
