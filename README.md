# Exploiting network analysis to create a novel sentinel surveillance system for efficient, rapid detection of emerging Clostridioides difficile strains in England

The following R scripts are provided:   
`1_Network_Identification.R`    
`2_Spread_Simulations.R`   
`3_Sentinel_Selection.R`    
`4_Sentinel_Set_Evaluation.R`    
`4appendix_Sentinel_Set_Evaluation.R`    
These correspond to the four steps in Figure 1 of the manuscript. Each uses output from the lower numbered scripts. 

Also provided is `0_Dummy_Data.R`.      
This generates dummy data files.
NB the facilities and network connections output by `HospitalNetwork::stay_linelist` are somewhat homogeneous (severely limiting the ability to optimise sentinel choice) so we recommend using primarily to assist with code familiarisation.  

The following R packages are called within these scripts:
`tidyverse`, `readr`, `janitor`, `igraph`, `tictoc`.  
The manuscript analysis used R version 4.4.1, on RStudio (version 2024.04.2), with `tidyverse 2.0.0`, `readr 2.1.5`, `janitor 2.2.0`, `igraph 2.0.3`, `tictoc 1.2.1` and a HospiNet object generated using `HospitalNetwork 0.9.3`.
