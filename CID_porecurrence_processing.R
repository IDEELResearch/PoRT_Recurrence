################################################
# Frequent Po recurrence in coastal Tanzania
# Author: Kelly Carey-Ewend and Guozheng Yang
# Date: 7 February 2025
# Modified: 3 July 2025
################################################


################################################
#######DATA LOADING AND CLEANING################
################################################


#set working directory
setwd("~/Documents/PoRT")

#load needed packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(ggsurvfit)
  library(ggplot2)
  library(icenReg)
  library(janitor)
  library(MCC)
  library(patchwork)
  library(RColorBrewer)
  library(bshazard)
  library(utile.visuals)
  library(gginnards)
  library(forcats)
})

#Determine colorblind-friendly color palettes
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n=12,"Paired")
brewer.pal(n=12,"Paired")

#load in dataset
port_full <- read.csv("CID_Porecurrence_dataset_deidentified.csv", header = TRUE)

#remove empty observations
port_missing <- port_full$screening_id[is.na(port_full$qpcr_reps_pos_screening_po_port)]
port <- port_full[!is.na(port_full$qpcr_reps_pos_screening_po_port),]

#label observations as being from Port or parallel PfTZ cohort
port <- port %>%
  mutate(cohort = ifelse(grepl("PfTZ", screening_id, fixed=TRUE), "pftz", "port"))

#evaluate redundant variables (to be cut, recoded, or merged). port_cut removes variables for easier parsing (not to be used in final analyses)
port_cut <- port[,-c(2:48)]

#determine whether any pf wk2 observations are only shown to be sampled in a single variable
for (i in 1:nrow(port_cut)){
  port_cut$wk2_conc[i] <- 0
  if (!is.na(port_cut$qpcr_reps_pos_wk2[i])){port_cut$wk2_conc[i] <- port_cut$wk2_conc[i] + 1}
  if (!is.na(port_cut$qpcr_reps_pos_wk2fue_ven[i])){port_cut$wk2_conc[i] <- port_cut$wk2_conc[i] + 1}
  if (!is.na(port_cut$qpcr_reps_pos_wk2_port[i])){port_cut$wk2_conc[i] <- port_cut$wk2_conc[i] + 1}
  if (!is.na(port_cut$qpcr_reps_pos_wk2fue_ven_port[i])){port_cut$wk2_conc[i] <- port_cut$wk2_conc[i] + 1}
}

#create list of redudant or unreliable variables to be cut, remove from port dataset
drop <- c("enrollment_age", "enrollment_gender", 
          "qpcr_pfdens_screening", "qpcr_pfdens_screening_port", 
          "qpcr_pfdens_enroll", "qpcr_pfdens_enroll_port",
          "qpcr_pfdens_enroll_cap", "qpcr_pfdens_enroll_cap_port",
          "qpcr_pfdens_wk2", "qpcr_pfdens_wk2fue_cap", "qpcr_reps_pos_wk2fue_cap", "qpcr_pfdens_wk2fue_ven", 
          "qpcr_pfdens_wk2_port", "qpcr_pfdens_wk2fue_cap_port", "qpcr_reps_pos_wk2fue_cap_port", "qpcr_pfdens_wk2fue_ven_port", 
          "qpcr_pfdens_wk4", "qpcr_pfdens_wk4_port", "qpcr_pfdens_wk4fue_cap", "qpcr_reps_pos_wk4fue_cap", 
          "qpcr_pfdens_wk4fue_ven", "qpcr_pfdens_wk4fue_cap_port", "qpcr_pfdens_wk4fue_ven_port", 
          "qpcr_pfdens_wk6_port", "qpcr_pfdens_wk8_port", 
          "qpcr_pfdens_wk10_port", "qpcr_pfdens_wk12_port",
          "qpcr_pfdens_wk14_port", "qpcr_pfdens_wk16_port",
          "qpcr_pfdens_wk18_port", "qpcr_pfdens_wk20_port",
          "qpcr_pfdens_wk22_port", "qpcr_pfdens_wk24_port"
          )

port <- port[ , !(names(port) %in% drop)]

#drop individuals who were not sampled at enrollment (or followed longitudinally)
port <- port[port$port_sampleyn_d0 == 1,]

##################################################################
###create merged, cleaned variables for processing and analyses###
##################################################################

###Po positivity and study missingness

#po_screening: po-positivity at screening (missing, neg {0}, pos {1})
port <- port %>%
  mutate(po_screening = ifelse(
    is.na(qpcr_reps_pos_screening_po) & is.na(qpcr_reps_pos_screening_po_port),
    NA, 
  #NA values cannot be compared to 0, so each evaluation must check whether is NA or a negative value so it will ignore NAs
  #observations with all NA are already counted as missing in the outer if-else
    ifelse(
      (is.na(qpcr_reps_pos_screening_po) | qpcr_reps_pos_screening_po == 0) 
      & 
      (is.na(qpcr_reps_pos_screening_po_port) | qpcr_reps_pos_screening_po_port == 0),
      0, 1)))

#missing_screening: no sampling at screening (sampled {0}, missing {1})
port <- port %>%
  mutate(missing_screening = ifelse(is.na(po_screening), 1, 0))
         
#po_enroll: po-positivity at enrollment (missing, neg {0}, pos {1})
port <- port %>%
  mutate(po_enroll = ifelse(
    is.na(qpcr_po_reps_pos_enroll) & is.na(qpcr_po_reps_pos_enroll_port), 
    NA,
    ifelse(
      (is.na(qpcr_po_reps_pos_enroll) | qpcr_po_reps_pos_enroll == 0) 
      & 
      (is.na(qpcr_po_reps_pos_enroll_port) | qpcr_po_reps_pos_enroll_port == 0),
      0, 1)))

#missing_enroll: no sampling at enrollment (sampled {0}, missing {1})
port <- port %>%
  mutate(missing_enroll = ifelse(is.na(po_enroll), 1, 0))

#po_wk2: po-positivity at wk2 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk2 = ifelse(
    is.na(qpcr_reps_pos_wk2_po) & is.na(qpcr_reps_pos_wk2_po_port) & is.na(qpcr_reps_pos_wk2fue_po) & is.na(qpcr_reps_pos_wk2fue_po_port), 
    NA,
    ifelse(
      (is.na(qpcr_reps_pos_wk2_po) | qpcr_reps_pos_wk2_po == 0) 
      & 
        (is.na(qpcr_reps_pos_wk2_po_port) | qpcr_reps_pos_wk2_po_port == 0)
      & 
        (is.na(qpcr_reps_pos_wk2fue_po) | qpcr_reps_pos_wk2fue_po == 0)
      &
        (is.na(qpcr_reps_pos_wk2fue_po_port) | qpcr_reps_pos_wk2fue_po_port == 0),
      0, 1)))
         
#missing_wk2: no sampling at wk2 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk2 = ifelse(is.na(po_wk2), 1, 0))

#po_wk4: po-positivity at wk4 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk4 = ifelse(
    is.na(qpcr_reps_pos_wk4_po) & is.na(qpcr_reps_pos_wk4_po_port) & is.na(qpcr_reps_pos_wk4fue_ven_po) & is.na(qpcr_reps_pos_wk4fue_ven_po_port), 
    NA,
    ifelse(
      (is.na(qpcr_reps_pos_wk4_po) | qpcr_reps_pos_wk4_po == 0) 
      & 
        (is.na(qpcr_reps_pos_wk4_po_port) | qpcr_reps_pos_wk4_po_port == 0)
      & 
        (is.na(qpcr_reps_pos_wk4fue_ven_po) | qpcr_reps_pos_wk4fue_ven_po == 0)
      &
        (is.na(qpcr_reps_pos_wk4fue_ven_po_port) | qpcr_reps_pos_wk4fue_ven_po_port == 0),
      0, 1)))

#missing_wk4: no sampling at wk4 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk4 = ifelse(is.na(po_wk4), 1, 0))

#po_wk6: po-positivity at wk6 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk6 = ifelse(
    is.na(qpcr_reps_pos_wk6_po_port), NA,
    ifelse(qpcr_reps_pos_wk6_po_port == 0, 0, 1)))

#missing_wk6: no sampling at wk6 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk6 = ifelse(is.na(po_wk6), 1, 0))

#po_wk8: po-positivity at wk8 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk8 = ifelse(
    is.na(qpcr_reps_pos_wk8_po_port), NA,
    ifelse(qpcr_reps_pos_wk8_po_port == 0, 0, 1)))

#missing_wk8: no sampling at wk8 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk8 = ifelse(is.na(po_wk8), 1, 0))

#po_wk10: po-positivity at wk10 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk10 = ifelse(
    is.na(qpcr_reps_pos_wk10_po_port), NA,
    ifelse(qpcr_reps_pos_wk10_po_port == 0, 0, 1)))

#missing_wk10: no sampling at wk10 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk10 = ifelse(is.na(po_wk10), 1, 0))

#po_wk12: po-positivity at wk12 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk12 = ifelse(
    is.na(qpcr_reps_pos_wk12_po_port), NA,
    ifelse(qpcr_reps_pos_wk12_po_port == 0, 0, 1)))

#missing_wk12: no sampling at wk12 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk12 = ifelse(is.na(po_wk12), 1, 0))

#po_wk14: po-positivity at wk14 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk14 = ifelse(
    is.na(qpcr_reps_pos_wk14_po_port), NA,
    ifelse(qpcr_reps_pos_wk14_po_port == 0, 0, 1)))

#missing_wk14: no sampling at wk14 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk14 = ifelse(is.na(po_wk14), 1, 0))

#po_wk16: po-positivity at wk16 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk16 = ifelse(
    is.na(qpcr_reps_pos_wk16_po_port), NA,
    ifelse(qpcr_reps_pos_wk16_po_port == 0, 0, 1)))

#missing_wk16: no sampling at wk16 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk16 = ifelse(is.na(po_wk16), 1, 0))

#po_wk18: po-positivity at wk18 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk18 = ifelse(
    is.na(qpcr_reps_pos_wk18_po_port), NA,
    ifelse(qpcr_reps_pos_wk18_po_port == 0, 0, 1)))

#missing_wk18: no sampling at wk18 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk18 = ifelse(is.na(po_wk18), 1, 0))

#po_wk20: po-positivity at wk20 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk20 = ifelse(
    is.na(qpcr_reps_pos_wk20_po_port), NA,
    ifelse(qpcr_reps_pos_wk20_po_port == 0, 0, 1)))

#missing_wk20: no sampling at wk20 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk20 = ifelse(is.na(po_wk20), 1, 0))

#po_wk22: po-positivity at wk22 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk22 = ifelse(
    is.na(qpcr_reps_pos_wk22_po_port), NA,
    ifelse(qpcr_reps_pos_wk22_po_port == 0, 0, 1)))

#missing_wk22: no sampling at wk22 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk22 = ifelse(is.na(po_wk22), 1, 0))

#po_wk24: po-positivity at wk24 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(po_wk24 = ifelse(
    is.na(qpcr_reps_pos_wk24_po_port), NA,
    ifelse(qpcr_reps_pos_wk24_po_port == 0, 0, 1)))

#missing_wk24: no sampling at wk24 (sampled {0}, missing {1}) 
port <- port %>%
  mutate(missing_wk24 = ifelse(is.na(po_wk24), 1, 0))


#####Pf-positivity

#pf_screening: pf-positivity at screening (missing, neg {0}, pos {1})
port <- port %>%
  mutate(pf_screening = ifelse(
    is.na(qpcr_reps_pos_screening) & is.na(qpcr_reps_pos_screening_port),
    NA, 
    ifelse(
      (is.na(qpcr_reps_pos_screening) | qpcr_reps_pos_screening == 0) 
      & 
        (is.na(qpcr_reps_pos_screening_port) | qpcr_reps_pos_screening_port == 0),
      0, 1)))

#pf_enroll: pf-positivity at enrollment (missing, neg {0}, pos {1})
port <- port %>% 
  #pcr
  mutate(pf_enroll = ifelse(
    is.na(qpcr_reps_pos_enroll) & is.na(qpcr_reps_pos_enroll_port) & is.na(qpcr_reps_pos_enroll_cap) & is.na(qpcr_reps_pos_enroll_cap_port), 
    NA,
    ifelse(
      (is.na(qpcr_reps_pos_enroll) | qpcr_reps_pos_enroll == 0) 
      & 
        (is.na(qpcr_reps_pos_enroll_port) | qpcr_reps_pos_enroll_port == 0)
      & 
        (is.na(qpcr_reps_pos_enroll_cap) | qpcr_reps_pos_enroll_cap == 0)
      &
        (is.na(qpcr_reps_pos_enroll_cap_port) | qpcr_reps_pos_enroll_cap_port == 0),
      0, 1)))

#pf_wk2: pf-positivity at wk2 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk2 = ifelse(
    is.na(qpcr_reps_pos_wk2) & is.na(qpcr_reps_pos_wk2fue_ven) & is.na(qpcr_reps_pos_wk2_port) & is.na(qpcr_reps_pos_wk2fue_ven_port), 
    NA,
    ifelse(
      (is.na(qpcr_reps_pos_wk2) | qpcr_reps_pos_wk2 == 0) 
      & 
        (is.na(qpcr_reps_pos_wk2fue_ven) | qpcr_reps_pos_wk2fue_ven == 0)
      & 
        (is.na(qpcr_reps_pos_wk2_port) | qpcr_reps_pos_wk2_port == 0)
      &
        (is.na(qpcr_reps_pos_wk2fue_ven_port) | qpcr_reps_pos_wk2fue_ven_port == 0),
      0, 1)))

#pf_wk4: pf-positivity at wk4 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk4 = ifelse(
    is.na(qpcr_reps_pos_wk4) & is.na(qpcr_reps_pos_wk4_port) & is.na(qpcr_reps_pos_wk4fue_ven) & is.na(qpcr_reps_pos_wk4fue_cap_port) & is.na(qpcr_reps_pos_wk4fue_ven_port), 
    NA,
    ifelse(
      (is.na(qpcr_reps_pos_wk4) | qpcr_reps_pos_wk4 == 0) 
      & 
        (is.na(qpcr_reps_pos_wk4_port) | qpcr_reps_pos_wk4_port == 0)
      & 
        (is.na(qpcr_reps_pos_wk4fue_ven) | qpcr_reps_pos_wk4fue_ven == 0)
      &
        (is.na(qpcr_reps_pos_wk4fue_cap_port) | qpcr_reps_pos_wk4fue_cap_port == 0)
      &
        (is.na(qpcr_reps_pos_wk4fue_ven_port) | qpcr_reps_pos_wk4fue_ven_port == 0),
      0, 1)))

#pf_wk6: po-positivity at wk6 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk6 = ifelse(
    is.na(qpcr_reps_pos_wk6_port), NA,
    ifelse(qpcr_reps_pos_wk6_port == 0, 0, 1)))

#pf_wk8: po-positivity at wk8 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk8 = ifelse(
    is.na(qpcr_reps_pos_wk8_port), NA,
    ifelse(qpcr_reps_pos_wk8_port == 0, 0, 1)))

#pf_wk10: po-positivity at wk10 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk10 = ifelse(
    is.na(qpcr_reps_pos_wk10_port), NA,
    ifelse(qpcr_reps_pos_wk10_port == 0, 0, 1)))

#pf_wk12: po-positivity at wk12 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk12 = ifelse(
    is.na(qpcr_reps_pos_wk12_port), NA,
    ifelse(qpcr_reps_pos_wk12_port == 0, 0, 1)))

#pf_wk14: po-positivity at wk14 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk14 = ifelse(
    is.na(qpcr_reps_pos_wk14_port), NA,
    ifelse(qpcr_reps_pos_wk14_port == 0, 0, 1)))

#pf_wk16: po-positivity at wk16 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk16 = ifelse(
    is.na(qpcr_reps_pos_wk16_port), NA,
    ifelse(qpcr_reps_pos_wk16_port == 0, 0, 1)))

#pf_wk18: po-positivity at wk18 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk18 = ifelse(
    is.na(qpcr_reps_pos_wk18_port), NA,
    ifelse(qpcr_reps_pos_wk18_port == 0, 0, 1)))

#pf_wk20: po-positivity at wk20 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk20 = ifelse(
    is.na(qpcr_reps_pos_wk20_port), NA,
    ifelse(qpcr_reps_pos_wk20_port == 0, 0, 1)))

#pf_wk22: po-positivity at wk22 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk22 = ifelse(
    is.na(qpcr_reps_pos_wk22_port), NA,
    ifelse(qpcr_reps_pos_wk22_port == 0, 0, 1)))

#pf_wk24: po-positivity at wk24 (missing, neg {0}, pos {1})
port <- port %>% 
  mutate(pf_wk24 = ifelse(
    is.na(qpcr_reps_pos_wk24_port), NA,
    ifelse(qpcr_reps_pos_wk24_port == 0, 0, 1)))

#confirm missingness by Po is concordant with missingess by Pf and by treatment
print(port[c(port$missing_wk2 == 1 & (!is.na(port$pf_wk2) | !is.na(port$port_sample_treatyn_w2))),c("screening_id", "missing_wk2", "pf_wk2", "port_sample_treatyn_w2")])
print(port[c(port$missing_wk2 == 0 & (is.na(port$pf_wk2) | is.na(port$port_sample_treatyn_w2))),c("screening_id", "missing_wk2", "pf_wk2", "port_sample_treatyn_w2")])
print(port[c(port$missing_wk4 == 1 & (!is.na(port$pf_wk4) | !is.na(port$port_sample_treatyn_w4))),c("screening_id", "missing_wk4", "pf_wk4", "port_sample_treatyn_w4")])
print(port[c(port$missing_wk4 == 0 & (is.na(port$pf_wk4) | is.na(port$port_sample_treatyn_w4))),c("screening_id", "missing_wk4", "pf_wk4", "port_sample_treatyn_w4")])
print(port[c(port$missing_wk6 == 1 & (!is.na(port$pf_wk6) | !is.na(port$port_sample_treatyn_w6))),c("screening_id", "missing_wk6", "pf_wk6", "port_sample_treatyn_w6")])
print(port[c(port$missing_wk6 == 0 & (is.na(port$pf_wk6) | is.na(port$port_sample_treatyn_w6))),c("screening_id", "missing_wk6", "pf_wk6", "port_sample_treatyn_w6")])
print(port[c(port$missing_wk8 == 1 & (!is.na(port$pf_wk8) | !is.na(port$port_sample_treatyn_w8))),c("screening_id", "missing_wk8", "pf_wk8", "port_sample_treatyn_w8")])
print(port[c(port$missing_wk8 == 0 & (is.na(port$pf_wk8) | is.na(port$port_sample_treatyn_w8))),c("screening_id", "missing_wk8", "pf_wk8", "port_sample_treatyn_w8")])
print(port[c(port$missing_wk10 == 1 & (!is.na(port$pf_wk10) | !is.na(port$port_sample_treatyn_w10))),c("screening_id", "missing_wk10", "pf_wk10", "port_sample_treatyn_w10")])
print(port[c(port$missing_wk10 == 0 & (is.na(port$pf_wk10) | is.na(port$port_sample_treatyn_w10))),c("screening_id", "missing_wk10", "pf_wk10", "port_sample_treatyn_w10")])
print(port[c(port$missing_wk12 == 1 & (!is.na(port$pf_wk12) | !is.na(port$port_sample_treatyn_w12))),c("screening_id", "missing_wk12", "pf_wk12", "port_sample_treatyn_w12")])
print(port[c(port$missing_wk12 == 0 & (is.na(port$pf_wk12) | is.na(port$port_sample_treatyn_w12))),c("screening_id", "missing_wk12", "pf_wk12", "port_sample_treatyn_w12")])
print(port[c(port$missing_wk14 == 1 & (!is.na(port$pf_wk14) | !is.na(port$port_sample_treatyn_w14))),c("screening_id", "missing_wk14", "pf_wk14", "port_sample_treatyn_w14")])
print(port[c(port$missing_wk14 == 0 & (is.na(port$pf_wk14) | is.na(port$port_sample_treatyn_w14))),c("screening_id", "missing_wk14", "pf_wk14", "port_sample_treatyn_w14")])
print(port[c(port$missing_wk16 == 1 & (!is.na(port$pf_wk16) | !is.na(port$port_sample_treatyn_w16))),c("screening_id", "missing_wk16", "pf_wk16", "port_sample_treatyn_w16")])
print(port[c(port$missing_wk16 == 0 & (is.na(port$pf_wk16) | is.na(port$port_sample_treatyn_w16))),c("screening_id", "missing_wk16", "pf_wk16", "port_sample_treatyn_w16")])
print(port[c(port$missing_wk18 == 1 & (!is.na(port$pf_wk18) | !is.na(port$port_sample_treatyn_w18))),c("screening_id", "missing_wk18", "pf_wk18", "port_sample_treatyn_w18")])
print(port[c(port$missing_wk18 == 0 & (is.na(port$pf_wk18) | is.na(port$port_sample_treatyn_w18))),c("screening_id", "missing_wk18", "pf_wk18", "port_sample_treatyn_w18")])
print(port[c(port$missing_wk20 == 1 & (!is.na(port$pf_wk20) | !is.na(port$port_sample_treatyn_w20))),c("screening_id", "missing_wk20", "pf_wk20", "port_sample_treatyn_w20")])
print(port[c(port$missing_wk20 == 0 & (is.na(port$pf_wk20) | is.na(port$port_sample_treatyn_w20))),c("screening_id", "missing_wk20", "pf_wk20", "port_sample_treatyn_w20")])
print(port[c(port$missing_wk22 == 1 & (!is.na(port$pf_wk22) | !is.na(port$port_sample_treatyn_w22))),c("screening_id", "missing_wk22", "pf_wk22", "port_sample_treatyn_w22")])
print(port[c(port$missing_wk22 == 0 & (is.na(port$pf_wk22) | is.na(port$port_sample_treatyn_w22))),c("screening_id", "missing_wk22", "pf_wk22", "port_sample_treatyn_w22")])
print(port[c(port$missing_wk24 == 1 & (!is.na(port$pf_wk24) | !is.na(port$port_sample_treatyn_w24))),c("screening_id", "missing_wk24", "pf_wk24", "port_sample_treatyn_w24")])
print(port[c(port$missing_wk24 == 0 & (is.na(port$pf_wk24) | is.na(port$port_sample_treatyn_w24))),c("screening_id", "missing_wk24", "pf_wk24", "port_sample_treatyn_w24")])


#merge final variables together into new dataset (dropping other redundant variables)
full_drops <- c(
  'qpcr_reps_pos_screening_po', 
  'qpcr_reps_pos_screening_po_port', 
  'qpcr_po_reps_pos_enroll', 
  'qpcr_po_reps_pos_enroll_port', 
  'qpcr_reps_pos_wk2_po', 
  'qpcr_reps_pos_wk2_po_port', 
  'qpcr_reps_pos_wk2fue_po', 
  'qpcr_reps_pos_wk2fue_po_port', 
  'qpcr_reps_pos_wk4_po', 
  'qpcr_reps_pos_wk4_po_port', 
  'qpcr_reps_pos_wk4fue_ven_po', 
  'qpcr_reps_pos_wk4fue_ven_po_port', 
  'qpcr_reps_pos_wk6_po_port', 
  'qpcr_reps_pos_wk8_po_port', 
  'qpcr_reps_pos_wk10_po_port', 
  'qpcr_reps_pos_wk12_po_port', 
  'qpcr_reps_pos_wk14_po_port', 
  'qpcr_reps_pos_wk16_po_port', 
  'qpcr_reps_pos_wk18_po_port', 
  'qpcr_reps_pos_wk20_po_port', 
  'qpcr_reps_pos_wk22_po_port', 
  'qpcr_reps_pos_wk24_po_port',
  'qpcr_pfdens_screening', 
  'qpcr_reps_pos_screening', 
  'qpcr_pfdens_screening_port', 
  'qpcr_reps_pos_screening_port', 
  'qpcr_pfdens_enroll', 
  'qpcr_reps_pos_enroll', 
  'qpcr_pfdens_enroll_port', 
  'qpcr_reps_pos_enroll_port', 
  'qpcr_pfdens_enroll_cap', 
  'qpcr_reps_pos_enroll_cap', 
  'qpcr_pfdens_enroll_cap_port', 
  'qpcr_reps_pos_enroll_cap_port', 
  'qpcr_pfdens_wk2', 
  'qpcr_reps_pos_wk2', 
  'qpcr_pfdens_wk2fue_cap', 
  'qpcr_reps_pos_wk2fue_cap', 
  'qpcr_pfdens_wk2fue_ven', 
  'qpcr_reps_pos_wk2fue_ven', 
  'qpcr_pfdens_wk2_port', 
  'qpcr_reps_pos_wk2_port', 
  'qpcr_pfdens_wk2fue_cap_port', 
  'qpcr_reps_pos_wk2fue_cap_port', 
  'qpcr_pfdens_wk2fue_ven_port', 
  'qpcr_reps_pos_wk2fue_ven_port', 
  'qpcr_pfdens_wk4', 
  'qpcr_reps_pos_wk4', 
  'qpcr_pfdens_wk4_port', 
  'qpcr_reps_pos_wk4_port', 
  'qpcr_pfdens_wk4fue_cap', 
  'qpcr_reps_pos_wk4fue_cap', 
  'qpcr_pfdens_wk4fue_ven', 
  'qpcr_reps_pos_wk4fue_ven', 
  'qpcr_pfdens_wk4fue_cap_port', 
  'qpcr_reps_pos_wk4fue_cap_port', 
  'qpcr_pfdens_wk4fue_ven_port', 
  'qpcr_reps_pos_wk4fue_ven_port', 
  'qpcr_pfdens_wk6_port', 
  'qpcr_reps_pos_wk6_port', 
  'qpcr_pfdens_wk8_port', 
  'qpcr_reps_pos_wk8_port', 
  'qpcr_pfdens_wk10_port', 
  'qpcr_reps_pos_wk10_port', 
  'qpcr_pfdens_wk12_port', 
  'qpcr_reps_pos_wk12_port', 
  'qpcr_pfdens_wk14_port', 
  'qpcr_reps_pos_wk14_port', 
  'qpcr_pfdens_wk16_port', 
  'qpcr_reps_pos_wk16_port', 
  'qpcr_pfdens_wk18_port', 
  'qpcr_reps_pos_wk18_port', 
  'qpcr_pfdens_wk20_port', 
  'qpcr_reps_pos_wk20_port', 
  'qpcr_pfdens_wk22_port', 
  'qpcr_reps_pos_wk22_port', 
  'qpcr_pfdens_wk24_port', 
  'qpcr_reps_pos_wk24_port'
)


port <- port[ , !(names(port) %in% full_drops)]


#determine individuals who dropped out (missed po visits at one point and all following timepoints) or were administratively censored
port <- port %>%
#first, create variable representing time of first visit missed in complete LTFU (no following visits)
  mutate(t_ltfu = ifelse(
    #nested if-else statements to determine individuals who were continuously missing after a certain point
    #first, LTFU at wk2
    missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
    & missing_wk16 == 1 & missing_wk14 == 1 & missing_wk12 == 1 & missing_wk10 == 1
    & missing_wk8 == 1 & missing_wk6 == 1 & missing_wk4 == 1  & missing_wk2 == 1,
    14, 
    #LTFU at wk4 
    ifelse(
      missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
      & missing_wk16 == 1 & missing_wk14 == 1 & missing_wk12 == 1 & missing_wk10 == 1
      & missing_wk8 == 1 & missing_wk6 == 1 & missing_wk4 == 1,
      28,
      #LTFU at wk6 
      ifelse(
        missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
        & missing_wk16 == 1 & missing_wk14 == 1 & missing_wk12 == 1 & missing_wk10 == 1
        & missing_wk8 == 1 & missing_wk6 == 1,
        42,
        #LTFU at wk8 
        ifelse(
          missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
          & missing_wk16 == 1 & missing_wk14 == 1 & missing_wk12 == 1 & missing_wk10 == 1
          & missing_wk8 == 1,
          56,
          #LTFU at wk10 
          ifelse(
            missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
            & missing_wk16 == 1 & missing_wk14 == 1 & missing_wk12 == 1 & missing_wk10 == 1,
            70,
            #LTFU at wk12 
            ifelse(
              missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
              & missing_wk16 == 1 & missing_wk14 == 1 & missing_wk12 == 1,
              84,
              #LTFU at wk14 
              ifelse(
                missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
                & missing_wk16 == 1 & missing_wk14 == 1,
                98,
                #LTFU at wk16 
                ifelse(
                  missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1
                  & missing_wk16 == 1,
                  112,
                  #LTFU at wk18 
                  ifelse(
                    missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1 & missing_wk18 == 1,
                    126,
                    #LTFU at wk20
                    ifelse(
                      missing_wk24 == 1 & missing_wk22 == 1 & missing_wk20 == 1,
                      140,
                      #LTFU at wk22
                      ifelse(
                        missing_wk24 == 1 & missing_wk22 == 1,
                        154,
                        #LTFU at wk24
                        ifelse(
                          missing_wk24 == 1,
                          168, NA
                      )
                     )
                    )
                   )
                  )
                 )
                )
               )
              )
             )
            )
           )
          )

################################################
####### DETERMINE LOSS TO FOLLOW UP ############
################################################

#create t_visit to show timing of last visit before LTFU

port <- port %>%
  mutate(t_visit = t_ltfu - 14)

#include only individuals screened through week 4
port <- port[(port$t_visit >= 28 | is.na(port$t_visit)),]

#create list of days corresponding to week2 - week24
days <- c(14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 168)

#initialize time to 1st po_recurrence
port <- port %>%
  mutate(t_po1st = NA)

#write csv of all mqtz ids for all port observations
write.csv(port[,c("screening_id", "qpcr_mqtzid_port")],"port_mqtz_ids.csv", row.names = TRUE)

################################################
####### PLOTTING HELPER FUNCTIONS ##############
################################################

#add helper functions for survminer survival curve plotting package
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

# add method to grid.draw (needed to ggsave survplot objects)
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}


################################################
###### TIME TO 1ST PO POSITIVE TEST ############
################################################

#create variable for first po recurrence (not accounting for negativity between)
for (i in 1:nrow(port)) {
  print(port$screening_id[i])
  #create vector of po status from wk2 to wk24
  screen <- c(port$po_wk2[i], port$po_wk4[i], port$po_wk6[i], port$po_wk8[i],
              port$po_wk10[i], port$po_wk12[i], port$po_wk14[i], port$po_wk16[i],
              port$po_wk18[i], port$po_wk20[i], port$po_wk22[i], port$po_wk24[i])
  for (j in 1:length(screen)) {
    if (screen[j] == 1 & !is.na(screen[j])) {
      port$t_po1st[i] <- days[j]
      break
    }
  }
}

#merge t_visit and t_po1st to create t, delta, and outcome
# t should represent first visit with Po, last visit prior to LTFU, or final visit at wk24 with no Po (whichever comes first)
#delta should be 1 if individual was LTFU at t and 0 for all else
#outcome po1st should be 1 if individual was po-positive at t and 0 for all else
#individuals who make it to the end of the study with no Po should have t=168, delta = 0, and po1st = 0
port <- port %>%
  mutate(t = ifelse(
    #if person was censored (t_visit exists) and had a po positive (t_po1st exists) and time to positive was before censoring
    !is.na(t_po1st) & !is.na(t_visit) & t_po1st <= t_visit 
    | 
    #or if the person was never censored and had a po positive
    !is.na(t_po1st) & is.na(t_visit), 
    #return time to first po positive
    t_po1st, 
    #otherwise
    ifelse(
      #if individual was censored, report time of last visit. otherwise, give end of study
    !is.na(t_visit), t_visit, 168
  ))) %>%
  #if individual is missing final LTFU visit, or if final LTFU visit is after first po positive, drop = 0. Otherwise, drop = 1
  mutate(drop = ifelse(is.na(t_visit) | !is.na(t_visit) & !is.na(t_po1st) & t_po1st < t_visit, 0, 1)) %>%
  #if had po positive, delta = 1
  mutate(delta = ifelse(!is.na(t_po1st), 1, 0))


print(port[,c("t_po1st", "t_visit", "t", "delta", "drop")])

#determine age cutpoint with maximal survival difference
surv_cutpoint(port, time = "t", event = "delta", c("enrollment_age_port"), minprop = 0.1, progressbar = TRUE)

#encode new age binary variable with cutpoint of 33 determined based on dataset
port <- port %>%
  mutate(age_cat_emp = ifelse(enrollment_age_port >= 31, 1, 0))

#encode new age binary variable with a priori cutoff at 18yo
port <- port %>%
  mutate(adult = ifelse(enrollment_age_port >= 18, 1, 0))

#subset to port observations only
port_only <- port[(port$cohort == "port"),]

#basic Kaplan-meier survival plot without accounting for censoring
km <- survfit(Surv(t, delta) ~ 1, data = port_only)
plot(km)
summary(km)
total_surv <- ggsurvplot(km, data = port_only, risk.table = TRUE, color = "red2", cumcensor = TRUE, fun = "event")
total_surv
ggsave("figures/surv_total_1strecur_nointerval.png", total_surv, height = 6, width = 6, dpi = 600)

#divide by initial variables of interest: sex, age, pf-positivity, treatment at start

#sex
km_sex <- survfit(Surv(t, delta) ~ enrollment_gender_port, data = port_only)
summary(km_sex)
plot(km_sex)
sex_surv <- ggsurvplot(km_sex, data = port_only, risk.table = TRUE, cumcensor = TRUE, pval = TRUE, conf.int = TRUE,
                     fun = "event",
                       palette = c("skyblue2","pink1"), legend.labs = c("Male","Female"))
sex_surv
ggsave("figures/surv_sex_1strecur_nointerval.png", sex_surv, height = 6, width = 6, dpi = 600)

#Cox proportional hazards model by sex
coxph_sex <- coxph(Surv(t, delta) ~ enrollment_gender_port, data = port_only)
summary(coxph_sex)

#pf-infection at enrollment
km_pf <- survfit(Surv(t, delta) ~ pf_enroll, data = port_only)
summary(km_pf)
plot(km_pf)
pf_surv <- ggsurvplot(km_pf, data = port_only, risk.table = TRUE, cumcensor = TRUE, pval = TRUE, conf.int = TRUE,
              fun = "event",
              palette = c("skyblue2","red2"), legend.labs = c("Pf-","Pf+"))
pf_surv
ggsave("figures/surv_pfenroll_1strecur_nointerval.png", pf_surv, height = 6, width = 6, dpi = 600)

#Cox proportional hazards model by pf-infection at enrollment
coxph_pf <- coxph(Surv(t, delta) ~ pf_enroll, data = port_only)
summary(coxph_pf)

#ACT treatment at enrollment
km_treat <- survfit(Surv(t, delta) ~ port_sample_treatyn_do, data = port_only)
summary(km_treat)
plot(km_treat)
treat_surv <- ggsurvplot(km_treat, data = port_only, risk.table = TRUE, cumcensor = TRUE, pval = TRUE, conf.int = TRUE,
                         fun = "event",
                         palette = c("skyblue2","red2"), legend.labs = c("untreated","treated"))
treat_surv
ggsave("figures/surv_treatenroll_1strecur_nointerval.png", treat_surv, height = 6, width = 6, dpi = 600)

#Cox proportional hazards model by treatment at enrollment
coxph_treat <- coxph(Surv(t, delta) ~ port_sample_treatyn_do, data = port_only)
summary(coxph_treat)


#plot survival between individuals over and under 31yo
km_age <- survfit(Surv(t, delta) ~ age_cat_emp, data = port_only)
summary(km_age)
age_surv <- ggsurvplot(km_age, data = port_only, risk.table = TRUE, cumcensor = TRUE, pval = TRUE, conf.int = TRUE,
                       fun = "event",
                       palette = c("skyblue2","red2"), legend.labs = c("under 31yo","31yo+"))
age_surv
ggsave("figures/surv_age_1strecur_nointerval.png", age_surv, height = 6, width = 6, dpi = 600)

#Cox PH model by age cutoff
coxph_age <- coxph(Surv(t, delta) ~ age_cat_emp, data = port_only)
summary(coxph_age)

#plot survival between children and adults
km_adult <- survfit(Surv(t, delta) ~ adult, data = port_only)
summary(km_adult)
adult_surv <- ggsurvplot(km_adult, data = port_only, risk.table = TRUE, cumcensor = TRUE, pval = TRUE, conf.int = TRUE,
                         fun = "event",
                         palette = c("skyblue2","red2"), legend.labs = c("child","adult"))
adult_surv
ggsave("figures/surv_adult_1strecur_nointerval.png", adult_surv, height = 6, width = 6, dpi = 600)

#Cox PH model between children and adults
coxph_adult <- coxph(Surv(t, delta) ~ adult, data = port_only)
summary(coxph_adult)


######################################################
##### INTERVAL CENSORING, REMOVE IMMORTAL TIME #######
######################################################

#interval censoring (Surv() type = "interval2"), for analysis that accounts for biweekly screening and missed visits prior to positivity
###intervals, reflecting known event-free time, defined as:
#(NA, time2] for left-censoring, none in this study
#[time, NA) for right censoring, like LTFU and admin censoring at end of study. Time reflects final timepoint with no event
# [time, time2] where time == time2 for exact event time (none in study)
# (time, time2] for interval censoring, with known event by time2 and final prior negative test at time

#approach: use two lists of variables (po_wkXX) and times (days per variable)
#iterate through the variables, changing time and time 2 based on missing and positive values of po_wkXX and draw days from time list

#create list of days corresponding to enrollment, then week2 - week24
days <- c(0, 14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 168)

#create empty time time variables in port
port <- port %>%
  mutate(time = NA) %>%
  mutate(time2 = NA) %>%
  mutate(time_enroll = NA) %>%
  mutate(time2_enroll = NA) %>%
  mutate(t_1stponeg = NA) %>%
  mutate(time_pf = NA) %>%
  mutate(time2_pf = NA) %>%
  mutate(time_enroll_pf = NA) %>%
  mutate(time2_enroll_pf = NA) %>%
  mutate(t_1stpfneg = NA)

#iterate through all individuals in port in order, adding intervals for po and pf positivity following their first negative
for (i in 1:nrow(port)) {
  #create vector of po status from day0 - wk24
  #if from port, day 0 is enrollment (those positive at screening and negative at enrollment are at risk for recurrence starting at day 0)
  #enrollment missing for PfTZ, so screening will be the day 0 timepoint
  if (port$cohort[i] == "port") {
    screen <- c(port$po_enroll[i], port$po_wk2[i], port$po_wk4[i], port$po_wk6[i], port$po_wk8[i],
              port$po_wk10[i], port$po_wk12[i], port$po_wk14[i], port$po_wk16[i],
              port$po_wk18[i], port$po_wk20[i], port$po_wk22[i], port$po_wk24[i])
  } else {
    screen <- c(port$po_screening[i], port$po_wk2[i], port$po_wk4[i], port$po_wk6[i], port$po_wk8[i],
              port$po_wk10[i], port$po_wk12[i], port$po_wk14[i], port$po_wk16[i],
              port$po_wk18[i], port$po_wk20[i], port$po_wk22[i], port$po_wk24[i])
    }
  #set flag variable as 0, which flips the first time the individual screens negative (and is at risk for recurrence)
  neg_flag = 0
  #for each item in the list of screening results. j can be used to index a screening result and the corresponding day count in the days list
  for (j in 1:length(screen)) {
    #the first time they are negative, flag turns to 1
    if (neg_flag == 0 & !is.na(screen[j]) & screen[j] == 0) {
      neg_flag <- 1
      port$t_1stponeg[i] <- days[j]}
    if (neg_flag == 1) {
        #if a given day is po-negative, update, time to day[j] to reflect most recent negative test
        if (!is.na(screen[j]) & screen[j] == 0) {
          port$time_enroll[i] <- days[j]
          port$time[i] <- days[j] - port$t_1stponeg[i]}
        #if a given day is po-positive, update time2 to day[j] to reflect first recurrence, 
        #then break for loop
        if (!is.na(screen[j]) & screen[j] == 1) {
          port$time2_enroll[i] <- days[j]
          port$time2[i] <- days[j] - port$t_1stponeg[i]
          break
        }
      }
    }
  #create vector of pf status from enrollment to wk2 - wk24
  if (port$cohort[i] == "port") {
  screen_pf <- c(port$pf_enroll[i], port$pf_wk2[i], port$pf_wk4[i], port$pf_wk6[i], port$pf_wk8[i],
              port$pf_wk10[i], port$pf_wk12[i], port$pf_wk14[i], port$pf_wk16[i],
              port$pf_wk18[i], port$pf_wk20[i], port$pf_wk22[i], port$pf_wk24[i])
  } else {
    screen_pf <- c(port$pf_screening[i], port$pf_wk2[i], port$pf_wk4[i], port$pf_wk6[i], port$pf_wk8[i],
                   port$pf_wk10[i], port$pf_wk12[i], port$pf_wk14[i], port$pf_wk16[i],
                   port$pf_wk18[i], port$pf_wk20[i], port$pf_wk22[i], port$pf_wk24[i])
  }
  #set flag variable as 0, which flips the first time the individual screens negative (and is at risk for recurrence)
  neg_flag_pf = 0
  for (j in 1:length(screen_pf)) {
    #the first time they are negative, flag turns to 1
    if (neg_flag_pf == 0 & !is.na(screen_pf[j]) & screen_pf[j] == 0) {
      neg_flag_pf <- 1
      port$t_1stpfneg[i] <- days[j]}
    if (neg_flag_pf == 1) {
      #if a given day is po-negative, update, time to day[j] to reflect most recent negative test
      if (!is.na(screen_pf[j]) & screen_pf[j] == 0) {
        port$time_enroll_pf[i] <- days[j]
        port$time_pf[i] <- days[j] - port$t_1stpfneg[i]}
      #if a given day is po-positive, update time2 to day[j] to reflect first recurrence, 
      #then break for loop
      if (!is.na(screen_pf[j]) & screen_pf[j] == 1) {
        port$time2_enroll_pf[i] <- days[j]
        port$time2_pf[i] <- days[j] - port$t_1stpfneg[i]
        break
      }
    }
  }
}

#examine distribution of missingness without immediately-following recurrence
print(port[,c("screening_id", "missing_wk2", "missing_wk4", "missing_wk6", "missing_wk8", "missing_wk10", "missing_wk12", "missing_wk14", "missing_wk16", "missing_wk18", "missing_wk20", "missing_wk22", "missing_wk24")])

#Determine number of individuals who were positive at screening but negative at enrollment
screen_pos <- length(port$screening_id[port$po_screening == 1 & !is.na(port$po_screening)])
snpos_enrollneg <- length(port$screening_id[port$po_screening == 1 & !is.na(port$po_screening) & port$po_enroll == 0 & !is.na(port$po_enroll)])
snpos_enrollneg/screen_pos

#determine duration of baseline po infection
sort(port$t_1stponeg)
quantile(port$t_1stponeg, probs = seq(0, 1, 0.25), type = 4)

#determine Ct among individuals with persistent and single-visit baseline infections
quantile(port$qpcr_ct_screening_po_port[port$t_1stponeg > 0 & port$cohort == "port"], probs = seq(0, 1, 0.25), type = 4)
quantile(port$qpcr_ct_screening_po_port[port$t_1stponeg == 0& port$cohort == "port"], probs = seq(0, 1, 0.25), type = 4)

######################################################
########### CONVENTIONAL SURVIVAL CODING #############
######################################################

#create conventional survival coding for coxph investigation
port <- port %>%
  mutate(timesurv = case_when(!is.na(time2) ~ time2, is.na(time2) ~ time)) %>%
  mutate(event = case_when(!is.na(time2) ~ 1, is.na(time2) ~ 0)) %>%
  mutate(drop = case_when(!is.na(time2) ~ 0, is.na(time2) ~ 1)) %>%
  mutate(weeksurv = timesurv/7) %>%
  mutate(timesurv_pf = case_when(!is.na(time2_pf) ~ time2_pf, is.na(time2_pf) ~ time_pf)) %>%
  mutate(event_pf = case_when(!is.na(time2_pf) ~ 1, is.na(time2_pf) ~ 0)) %>%
  mutate(drop = case_when(!is.na(time2_pf) ~ 0, is.na(time2_pf) ~ 1)) %>%
  mutate(weeksurv_pf = timesurv_pf/7) %>%
  mutate(timesurv_enroll = case_when(!is.na(time2_enroll) ~ time2_enroll, is.na(time2_enroll) ~ time_enroll)) %>%
  mutate(event_enroll = case_when(!is.na(time2_enroll) ~ 1, is.na(time2_enroll) ~ 0))

#encode time from enrollment for plotting with delayed entry
#time on the x-axis will be absolute
#delayed entry coded using counting process structure

#determine first day of screening
print(port$port_sample_date_d0)
#"2020-10-15"

#encode new time variable to reflect absolute time since study entry
port <- port %>%
  mutate(timesurv1_abs = as.numeric(as.Date(port_sample_date_d0) + t_1stponeg - as.Date("2020-10-15"))) %>%
  mutate(timesurv2_abs = timesurv + timesurv1_abs) %>%
  mutate(weeksurv1_abs = timesurv1_abs/7) %>%
  mutate(weeksurv2_abs = timesurv2_abs/7)

######################################################
################ COVARIATE CODING ####################
######################################################

#calculate median po18S Ct of baseline infection
port <- port %>%
  mutate(po18s_ct = ifelse(qpcr_ct_screening_po_port > 0, qpcr_ct_screening_po_port, qpcr_po_ct_enroll_port))
quantile(port[port$cohort == "port", c("po18s_ct")], seq = c(0,1, by = 0.25))

#calculate median time between screening and enrollment
port <- port %>%
  mutate(screen_to_enroll_d = as.numeric(as.Date(port_sample_date_d0) - as.Date(qpcr_mqtz_date_port)))
quantile(port[port$cohort == "port", c("screen_to_enroll_d")], seq = c(0,1, by = 0.25))

#create variable reflecting pf positivity at screening, enrollment, wk2, or wk4
port <- port %>%
  mutate(pf_early = ifelse((!is.na(pf_screening) & pf_screening == 1) |
                             (!is.na(pf_enroll) & pf_enroll == 1) | 
                             (!is.na(pf_wk2) & pf_wk2 == 1) |
                             (!is.na(pf_wk4) & pf_wk4 == 1), 1, 0)) %>%
  mutate(pf_baseline = ifelse((!is.na(pf_screening) & pf_screening == 1) |
                             (!is.na(pf_enroll) & pf_enroll == 1), 1, 0))
sum(port[port$cohort == "port",c("pf_baseline")])

#create variable reflecting persistent po parasitemia from screening to week 4
port <- port %>%
  mutate(po_cont = ifelse(
    #all combinations of being pf positive early and continously positive throughb wk4
    #no skips allowed, as then a subsequent positive counts as this variable AND an early incident infection
      (!is.na(po_screening) & po_screening == 1 & !is.na(po_enroll) & po_enroll == 1) |
      (!is.na(po_screening) & po_screening == 1 & !is.na(po_enroll) & po_enroll == 1 & !is.na(po_wk2) & po_wk2 == 1) |
      (!is.na(po_screening) & po_screening == 1 & !is.na(po_enroll) & po_enroll == 1 & !is.na(po_wk2) & po_wk2 == 1 & !is.na(po_wk4) & po_wk4 == 1)
    , 1, 0)
  )

#encode age categories
port <- port %>%
  mutate(adult = ifelse(enrollment_age_port >= 18, 1, 0)) %>%
  mutate(over15 = ifelse(enrollment_age_port >= 16, 1, 0))

quantile(port$enrollment_age_port)

#encode season (w 6w lag) of screening
port <- port %>%
  mutate(date_6w = as.Date(port_sample_date_d0) + 42) %>%
  mutate(month_6w = format(as.Date(date_6w, format="%Y-%m-%d"),"%m")) %>%
  mutate(wetdry_6w = ifelse(month_6w %in% c("03", "04", "05", "10", "11" ,"12"), 1, 0)) %>%
  mutate(season_6w = ifelse(month_6w %in% c("03", "04", "05"), "masika", ifelse(month_6w %in% c("10", "11" ,"12"), "vuli", "dry")))


#define early recurrences
###identify and evaluate "early recurrences" likely to represent persistent parasitemia dipping beneath LLD
#individuals with t_1stponeg at enrollment, wk2, or wk4, who are then positive at wk2, wk4, or wk6
#first, examine individuals 1st negative at enrollment (positive at screening)
print(port[port$t_1stponeg == 0, c("po_enroll", "po_wk2", "po_wk4", "po_wk6")])
#6 quick recurrences at week 2, none persistent to week 4
#second, individuals 1st negative at wk2
print(port[port$t_1stponeg == 14, c("po_enroll", "po_wk2", "po_wk4", "po_wk6")])
#1 additional quick recurrence at w4 after being negative 1st at wk2
#finaly, individuals 1st negative at wk4
print(port[port$t_1stponeg == 28, c("po_enroll", "po_wk2", "po_wk4", "po_wk6")])
#no early recurrences

#label individuals with early recurrences
port <- port %>%
  mutate(po_wk2_early = case_when(
    po_enroll == 0 & po_wk2 == 1 ~ 1, 
    TRUE ~ 0)) %>%
  mutate(po_wk4_early = case_when(
    po_enroll == 1 & po_wk2 == 0 & po_wk4 == 1 ~ 1,
    TRUE ~ 0
  ))

####recalculate interval censoring with early recurrences excluded
#change positivity for possible persistent parasitemia
port_noearly <- port %>%
  mutate(po_enroll = case_when(po_wk2_early == 1 ~ 1, TRUE ~ po_enroll)) %>%
  mutate(po_wk2 = case_when(po_wk4_early == 1 ~ 1, TRUE ~ po_wk2)) %>%
  #clear variables
  mutate(time = NA) %>%
  mutate(time2 = NA) %>%
  mutate(t_1stponeg = NA) %>%
  mutate(time_enroll = NA) %>%
  mutate(time2_enroll = NA)

#investigate associations with po recurrence
summary(coxph(Surv(timesurv, event) ~ enrollment_gender_port + over15 + pf_early + po_cont, data = port))

#iterate through all individuals in port in order, adding intervals for po and pf positivity following their first negative
for (i in 1:nrow(port_noearly)) {
  #create vector of po status from day0 - wk24
  #if from port, day 0 is enrollment (those positive at screening and negative at enrollment are at risk for recurrence starting at day 0)
  #enrollment missing for PfTZ, so screening will be the day 0 timepoint
  if (port_noearly$cohort[i] == "port") {
    screen <- c(port_noearly$po_enroll[i], port_noearly$po_wk2[i], port_noearly$po_wk4[i], port_noearly$po_wk6[i], port_noearly$po_wk8[i],
                port_noearly$po_wk10[i], port_noearly$po_wk12[i], port_noearly$po_wk14[i], port_noearly$po_wk16[i],
                port_noearly$po_wk18[i], port_noearly$po_wk20[i], port_noearly$po_wk22[i], port_noearly$po_wk24[i])
  } else {
    screen <- c(port_noearly$po_screening[i], port_noearly$po_wk2[i], port_noearly$po_wk4[i], port_noearly$po_wk6[i], port_noearly$po_wk8[i],
                port_noearly$po_wk10[i], port_noearly$po_wk12[i], port_noearly$po_wk14[i], port_noearly$po_wk16[i],
                port_noearly$po_wk18[i], port_noearly$po_wk20[i], port_noearly$po_wk22[i], port_noearly$po_wk24[i])
  }
  #set flag variable as 0, which flips the first time the individual screens negative (and is at risk for recurrence)
  neg_flag = 0
  #for each item in the list of screening results. j can be used to index a screening result and the corresponding day count in the days list
  for (j in 1:length(screen)) {
    #the first time they are negative, flag turns to 1
    if (neg_flag == 0 & !is.na(screen[j]) & screen[j] == 0) {
      neg_flag <- 1
      port_noearly$t_1stponeg[i] <- days[j]}
    if (neg_flag == 1) {
      #if a given day is po-negative, update, time to day[j] to reflect most recent negative test
      if (!is.na(screen[j]) & screen[j] == 0) {
        port_noearly$time_enroll[i] <- days[j]
        port_noearly$time[i] <- days[j] - port_noearly$t_1stponeg[i]}
      #if a given day is po-positive, update time2 to day[j] to reflect first recurrence, 
      #then break for loop
      if (!is.na(screen[j]) & screen[j] == 1) {
        port_noearly$time2_enroll[i] <- days[j]
        port_noearly$time2[i] <- days[j] - port_noearly$t_1stponeg[i]
        break
      }
    }
  }
}

#recode for dataset without early recurrences
port_noearly <- port_noearly %>%
  mutate(timesurv = case_when(!is.na(time2) ~ time2, is.na(time2) ~ time)) %>%
  mutate(event = case_when(!is.na(time2) ~ 1, is.na(time2) ~ 0)) %>%
  mutate(drop = case_when(!is.na(time2) ~ 0, is.na(time2) ~ 1)) %>%
  mutate(weeksurv = timesurv/7)

#when was next recurrence? after treatment?
#PoTZ-010: never, treated wk6
print(port_noearly[port_noearly$screening_id == "PoTZ-010",])
#PoTZ-019: wk6, then treated at wk8
print(port_noearly[port_noearly$screening_id == "PoTZ-019",])
#PoTZ-026: wk24, treated at wk4
print(port_noearly[port_noearly$screening_id == "PoTZ-026",])
#PoTZ-031: never, never treated
print(port_noearly[port_noearly$screening_id == "PoTZ-031",])
#PoTZ-033: wk24, treated at wk4
print(port_noearly[port_noearly$screening_id == "PoTZ-033",])
#PoTZ-034: wk22, never treated
print(port_noearly[port_noearly$screening_id == "PoTZ-034",])
#PoTZ-069: never, treated at week 2
print(port_noearly[port_noearly$screening_id == "PoTZ-069",])


#add flag variable for early recurrences, then merge datasets
port <- port %>%
  mutate(early_flag = 0)

port_noearly <- port_noearly %>%
  mutate(early_flag = 1)

######################################################
############## PO RECURRENCE PLOTS ###################
######################################################

#subset to port observations only
port_only <- port[(port$cohort == "port"),]

#subset to port observations only
port_noearly <- port_noearly[(port_noearly$cohort == "port"),]

port_earlyandnoearly <- rbind(port_only, port_noearly)

#basic Kaplan-meier survival plot from enrollment with interval censoring
km <- survfit(Surv(time_enroll, time2_enroll, type = "interval2") ~ 1, data = port_only)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port_only, risk.table = TRUE, color = "red2", fun = "event", cumcensor = TRUE)
total_surv
ggsave("figures/surv_total_1strecurenroll_incen.png", total_surv, height = 6, width = 6, dpi = 600)

#KM from study enrollment with
km <- survfit(Surv(weeksurv1_abs, weeksurv2_abs, event) ~ 1, data = port_only)
plot(km)
summary(km)
sum(km$n.event)
abs_surv <- ggsurvplot(km, data = port_only, 
                         risk.table = TRUE, 
                         color = "#1F78B4", 
                         cumcensor = TRUE, 
                         fun = "event",
                         ylab = "Cumulative Po incidence",
                         xlab = "Time (w) since study start",
                         censor.shape = "X",
                         break.x.by = 12,
                         ylim = c(0,1),
                         xlim = c(0, 112),
                         legend.title = "") 

#basic Kaplan-meier survival plot from 1st po neg with interval censoring
km <- survfit(Surv(time, time2, type = "interval2") ~ 1, data = port_only)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port_only, 
                         risk.table = TRUE, 
                         color = "red2", 
                         cumcensor = TRUE, 
                         fun = "event",
                         xlab = "Time (d) since first negative test")
total_surv
ggsave("figures/surv_total_incen.png", total_surv, height = 6, width = 6, dpi = 600)

#basic Kaplan-meier survival plot from 1st po neg with conventional coding
km <- survfit(Surv(weeksurv, event, type = "right") ~ 1, data = port_only)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port_only, 
                         risk.table = TRUE, 
                         color = "#1F78B4", 
                         cumcensor = TRUE, 
                         cumevents = TRUE,
                         fun = "event",
                         ylab = "Cumulative incidence",
                         xlab = "Time (w) since first negative test",
                         censor.shape = "X",
                         break.x.by = 2,
                         ylim = c(0,1),
                         legend.title = "")
total_surv_plot <- total_surv$plot/total_surv$table/total_surv$cumevents/total_surv$ncensor.plot + plot_layout(heights=c(8,2, 2, 2))
total_surv_plot
ggsave("figures/surv_total.png", total_surv, height = 8, width = 6, dpi = 600)

#KM conventional coding to 1st Po recurrence, exclude early recurrences in one plot
km <- survfit(Surv(weeksurv, event, type = "right") ~ early_flag, data = port_earlyandnoearly)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port_earlyandnoearly, 
                         risk.table = TRUE, 
                         palette = c("#1F78B4", "#B15928"), 
                         legend.labs = c("All", "2 neg. tests req."),
                         cumcensor = TRUE, 
                         conf.int = TRUE,
                         fun = "event",
                         ylab = "Cumulative incidence",
                         xlab = "Time (w) since first negative test",
                         censor.shape = "X",
                         break.x.by = 2,
                         ylim = c(0,1),
                         legend.title = "")
total_surv_plot <- total_surv$plot/total_surv$table/total_surv$ncensor.plot + plot_layout(heights=c(8,2, 2))
total_surv_plot
ggsave("figures/surv_total_earlyrecur.png", total_surv_plot, height = 6, width = 7, dpi = 600)

#divide by initial variables of interest: sex, age, pf-positivity, treatment at start

###sex interval-censored
km_sex <- survfit(Surv(time, time2, type = "interval2") ~ enrollment_gender_port, data = port_only)
summary(km_sex)
plot(km_sex)
sex_surv <- ggsurvplot(km_sex, data = port_only, 
                       risk.table = TRUE, 
                       cumcensor = TRUE, 
                       conf.int = TRUE,
                       fun = "event",
                       palette = c("skyblue2","pink1"), 
                       legend.labs = c("Male","Female"),
                       xlab = "Time (d) since first negative test")
sex_surv
ggsave("figures/surv_sex_incen.png", sex_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
sex_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ enrollment_gender_port, port_only, 
                     model = "ph", 
                     weights = NULL,
                     B = c(0,1))
summary(sex_surv_sp)
sex_data <- data.frame(enrollment_gender_port = c(0,1))
rownames(sex_data) <- c('female', 'male')
plot(sex_surv_sp, sex_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ enrollment_gender_port, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)

#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
sex_surv_weibull <- ic_par(Surv(time, time2, type = "interval2") ~ enrollment_gender_port, port_only, model = "ph", dist = "weibull", weights = NULL)
sex_surv_weibull
plot(sex_surv_weibull)

###sex conventional survival coding
km_sex <- survfit(Surv(weeksurv, event, type = "right") ~ enrollment_gender_port, data = port_only)
plot(km_sex)
summary(km_sex)
sum(km_sex$n.event)
sex_surv <- ggsurvplot(km_sex, data = port_only, 
                         risk.table = TRUE, 
                         palette = c("skyblue2","pink1"), 
                         legend.labs = c("Male", "Female"),
                         #cumcensor = TRUE, 
                         conf.int = TRUE,
                         fun = "event",
                         ylab = "Cumulative incidence",
                         xlab = "Time (w) since first negative test",
                         censor.shape = "",
                         break.x.by = 2,
                         ylim = c(0,1),
                         legend.title = "")
sex_surv_plot <- sex_surv$plot/sex_surv$table + plot_layout(heights=c(6,2))
sex_surv_plot
#log-rank test p-value 0.3
ggsave("figures/surv_sex.png", sex_surv_plot, height = 4, width = 5, dpi = 600)

#test risk difference at final timepoint
#subset to final week
sex_24w <- summary(km_sex, times = 24)
#risk difference
rd_sex <- sex_24w$surv[1] - sex_24w$surv[2]
#standard error of RD
diffSE_sex <- sqrt(sex_24w$std.err[1]^2 + sex_24w$std.err[2]^2)
#z-test
zStat_sex <- rd_sex/diffSE_sex
2*pnorm(abs(zStat_sex), lower.tail=FALSE)

coxph(Surv(weeksurv, event, type = "right") ~ enrollment_gender_port, data = port_only)


#pf-positivity at screening, enrollment, wk2, or wk4
km_pf <- survfit(Surv(time, time2, type = "interval2") ~ pf_early, data = port_only)
summary(km_pf)
plot(km_pf)
pf_surv <- ggsurvplot(km_pf, data = port_only, 
                      risk.table = TRUE, 
                      cumcensor = TRUE, 
                      conf.int = TRUE,
                      fun = "event",
                      palette = c("skyblue2","red2"), 
                      legend.labs = c("Pf-","Pf+"),
                      xlab = "Time (d) since first negative test")
pf_surv
ggsave("figures/surv_pfearly_incen.png", pf_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by early pf status
pf_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ pf_early, port_only, 
                     model = "ph", 
                     weights = NULL,
                     B = c(0,1))
summary(pf_surv_sp)
pf_data <- data.frame(pf_early = c(0,1))
rownames(pf_data) <- c('Pf-', 'Pf+')
plot(pf_surv_sp, pf_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ pf_early, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ pf_early, port_only, model = "ph", dist = "weibull", weights = NULL)

#pf-positivity at screening, enrollment, wk2, or wk4, including PfTZ
km_pf <- survfit(Surv(time, time2, type = "interval2") ~ pf_early + cohort, data = port)
summary(km_pf)
plot(km_pf)
pf_surv <- ggsurvplot(km_pf, data = port, 
                      risk.table = TRUE, 
                      cumcensor = TRUE, 
                      conf.int = TRUE,
                      fun = "event",
                      palette = c("skyblue2","green4", "red2"), 
                      legend.labs = c("Pf-", "PfTZ (Pf+)", "Pf+"),
                      xlab = "Time (d) since first negative test")
pf_surv
ggsave("figures/surv_pfearly_portpftz_incen.png", pf_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by early pf status
pf_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ pf_early + cohort, data = port, 
                    model = "ph", 
                    weights = NULL,
                    B = c(0,1))
summary(pf_surv_sp)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ pf_early + cohort, data = port, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ pf_early + cohort, data = port, model = "ph", dist = "weibull", weights = NULL)


#plot survival between individuals over and under 18yo
km_adult <- survfit(Surv(time, time2, type = "interval2") ~ adult, data = port_only)
summary(km_adult)
adult_surv <- ggsurvplot(km_adult, data = port_only, 
                       risk.table = TRUE, 
                       cumcensor = TRUE,
                       conf.int = TRUE,
                       fun = "event",
                       palette = c("skyblue2","red2"), 
                       legend.labs = c("under 18yo","18yo+"),
                       xlab = "Time (d) since first negative test")
adult_surv
ggsave("figures/surv_adult_incen.png", adult_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
adult_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ adult, port_only, 
                    model = "ph", 
                    weights = NULL,
                    B = c(0,1))
summary(adult_surv_sp)
adult_data <- data.frame(adult = c(0,1))
rownames(adult_data) <- c('<18yo', '>=18yo')
plot(adult_surv_sp, adult_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ adult, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ adult, port_only, model = "ph", dist = "weibull", weights = NULL)

#plot survival between individuals over and under 16yo
km_over15 <- survfit(Surv(time, time2, type = "interval2") ~ over15, data = port_only)
summary(km_over15)
over15_surv <- ggsurvplot(km_over15, data = port_only, 
                         risk.table = TRUE, 
                         cumcensor = TRUE, 
                         conf.int = TRUE,
                         fun = "event",
                         palette = c("skyblue2","red2"), 
                         legend.labs = c("under 16yo","16yo+"),
                         xlab = "Time (d) since first negative test")
over15_surv
ggsave("figures/surv_over15_incen.png", over15_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
over15_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ over15, port_only, 
                       model = "ph", 
                       weights = NULL,
                       B = c(0,1))
summary(over15_surv_sp)
over15_data <- data.frame(over15 = c(0,1))
rownames(over15_data) <- c('<16yo', '>=16yo')
plot(over15_surv_sp, over15_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ over15, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ over15, port_only, model = "ph", dist = "weibull", weights = NULL)

#conventional survival coding of over-under 15yo
km_over15 <- survfit(Surv(weeksurv, event) ~ over15, data = port_only)
summary(km_over15)
over15_surv <- ggsurvplot(km_over15, data = port_only, 
                          risk.table = TRUE, 
                          palette = c("#33A02C","#FF7F00"), 
                          legend.labs = c("Under 16yo","16yo+"),
                          #cumcensor = TRUE, 
                          conf.int = TRUE,
                          fun = "event",
                          ylab = "Cumulative incidence",
                          xlab = "Time (w) since first negative test",
                          censor.shape = "",
                          break.x.by = 2,
                          ylim = c(0,1),
                          legend.title = "")
over15_surv_main <- over15_surv$plot/over15_surv$table + plot_layout(heights = c(6,2))
over15_surv_main
ggsave("figures/surv_over15.png", over15_surv_main, height = 4, width = 5, dpi = 600)

#test risk difference at final timepoint
#subset to final week
over15_24w <- summary(km_over15, times = 24)
#risk difference
rd_over15 <- over15_24w$surv[1] - over15_24w$surv[2]
#standard error of RD
diffSE_over15 <- sqrt(over15_24w$std.err[1]^2 + over15_24w$std.err[2]^2)
#z-test
zStat_over15 <- rd_over15/diffSE_over15
2*pnorm(abs(zStat_over15), lower.tail=FALSE)

###by season of screening

#survival anlysis by wet or dry season 6w prior to enrollment date
km_wetdry <- survfit(Surv(time, time2, type = "interval2") ~ wetdry_6w, data = port_only)
summary(km_wetdry)
plot(km_wetdry)
wetdry_surv <- ggsurvplot(km_wetdry, data = port_only, 
                          risk.table = TRUE, 
                          cumcensor = TRUE, 
                          conf.int = TRUE,
                          fun = "event",
                          palette = c("skyblue2","red2"), 
                          legend.labs = c("Dry","Wet"),
                          xlab = "Time (d) since first negative test")
wetdry_surv
ggsave("figures/surv_wetdry6w_incen.png", wetdry_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
wetdry_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ wetdry_6w, port_only, 
                        model = "ph", 
                        weights = NULL,
                        B = c(0,1))
summary(wetdry_surv_sp)
wetdry_data <- data.frame(wetdry_6w = c(0,1))
rownames(wetdry_data) <- c('Dry', 'Wet')
plot(wetdry_surv_sp, wetdry_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ wetdry_6w, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ wetdry_6w, port_only, model = "ph", dist = "weibull", weights = NULL)

#survival analysis to 1st Po recurrence by season
km_season <- survfit(Surv(time, time2, type = "interval2") ~ season_6w, data = port_only)
summary(km_season)
plot(km_season)
season_surv <- ggsurvplot(km_season, data = port_only, 
                          risk.table = TRUE, 
                          cumcensor = TRUE, 
                          conf.int = TRUE,
                          fun = "event",
                          palette = c("skyblue2","green4", "red2"), 
                          legend.labs = c("Dry","Masika", "Vuli"),
                          xlab = "Time (d) since first negative test")
season_surv
ggsave("figures/surv_season6w_incen.png", season_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
season_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ season_6w, port_only, 
                        model = "ph", 
                        weights = NULL,
                        B = c(0,1))
summary(season_surv_sp)
#season_data <- data.frame(season_6w = c("Dry", "Masika", "Vuli"))
#rownames(season_data) <- c("Dry", "Masika", "Vuli")
#plot(season_surv_sp, season_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ season_6w, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ season_6w, port_only, model = "ph", dist = "weibull", weights = NULL)


######################################################
######## PO LINE DIAGRAM  ############################
######################################################

#create event variable (1 if time2 is not missing)
port <- port %>%
  mutate(po_recur1 = ifelse(!is.na(time2),1,0))

#subset to port observations only
port_only <- port[(port$cohort == "port"),]

#create variable for plotting in row
port_only <- port_only %>%
  mutate(id = row_number())

line_plot <- ggplot(aes(x = id), data=port_only) +
  geom_segment(aes(xend = id, yend = t_1stponeg/7), y = 0, linetype = "dotted") +
  #for follow-up, plot to time2_enroll if individual had event and until time_enroll if censored
  geom_segment(aes(xend = id, y = t_1stponeg/7, yend = ifelse(!is.na(time2_enroll), time2_enroll/7, time_enroll/7))) +
  geom_point(aes(y = time2_enroll/7, size = po_recur1), shape = 16, show.legend = FALSE) +
  scale_size_continuous(range = c(0, 2)) +
  scale_y_continuous(breaks = seq(0,24, by = 2)) +
  ylab("Time (w) since enrollment") +
  xlab("Participant") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        panel.grid = element_blank())
line_plot

port_sort <- port_only[order(port_only$t_1stponeg, -port_only$time2_enroll),] %>%
  mutate(id = row_number())


line_plot_wide <- ggplot(aes(x = id), data=port_sort) +
  geom_segment(aes(xend = id, yend = t_1stponeg/7), y = 0, linetype = "dotted") +
  #for follow-up, plot to time2_enroll if individual had event and until time_enroll if censored
  geom_segment(aes(xend = id, y = t_1stponeg/7, yend = ifelse(!is.na(time2_enroll), time2_enroll/7, time_enroll/7))) +
  geom_point(aes(y = time2_enroll/7, size = po_recur1), shape = 16, show.legend = FALSE) +
  scale_size_continuous(range = c(0, 2)) +
  scale_y_continuous(breaks = seq(0,24, by = 2)) +
  ylab("Time (w) since enrollment") +
  xlab("Participant") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        panel.grid = element_blank())
line_plot_wide

ggsave("figures/surv_linediagram.png", line_plot, height = 7, width = 6, dpi = 600)
ggsave("figures/surv_linediagram_wide.png", line_plot_wide, height = 6, width = 9, dpi = 600)

######################################################
######## INTERVAL CENSORED PF RECURRENCE #############
######################################################

#create variable reflecting persistent pf positive at beginning of study period,
#with no skips or negative tests that would be counted as an early incident infection
port <- port %>%
  mutate(pf_cont = ifelse(
    #all combinations of being pf positive early and continously positive throughb wk4
    #no skips allowed, as then a subsequent positive counts as this variable AND an early incident infection
    (!is.na(pf_screening) & pf_screening == 1 ) |
    (!is.na(pf_screening) & pf_screening == 1 & !is.na(pf_enroll) & pf_enroll == 1) |
    (!is.na(pf_screening) & pf_screening == 1 & !is.na(pf_enroll) & pf_enroll == 1 & !is.na(pf_wk2) & pf_wk2 == 1) |
    (!is.na(pf_screening) & pf_screening == 1 & !is.na(pf_enroll) & pf_enroll == 1 & !is.na(pf_wk2) & pf_wk2 == 1 & !is.na(pf_wk4) & pf_wk4 == 1)
    , 1, 0)
    )
 
#subset to port observations only
port_only <- port[(port$cohort == "port"),]

#basic Kaplan-meier survival plot from 1st pf neg with interval censoring
km_pf <- survfit(Surv(time_pf, time2_pf, type = "interval2") ~ 1, data = port_only)
plot(km_pf)
summary(km_pf)
sum(km_pf$n.event)
total_surv_pf <- ggsurvplot(km_pf, data = port_only, 
                         risk.table = TRUE, 
                         color = "red2", 
                         cumcensor = TRUE, 
                         fun = "event",
                         xlab = "Time (d) since first negative test")
total_surv_pf
ggsave("figures/surv_total_pf_incen.png", total_surv_pf, height = 6, width = 6, dpi = 600)

###sex
km_sex_pf <- survfit(Surv(time_pf, time2_pf, type = "interval2") ~ enrollment_gender_port, data = port_only)
summary(km_sex_pf)
plot(km_sex_pf)
sex_surv_pf <- ggsurvplot(km_sex_pf, data = port_only, 
                       risk.table = TRUE, 
                       cumcensor = TRUE, 
                       conf.int = TRUE,
                       fun = "event",
                       palette = c("skyblue2","pink1"), 
                       legend.labs = c("Male","Female"),
                       xlab = "Time (d) since first negative test")
sex_surv_pf
ggsave("figures/surv_sex_pf_incen.png", sex_surv_pf, height = 6, width = 6, dpi = 600)


#semi-parametric model with proportional hazards by sex
sex_surv_pf_sp <- ic_sp(Surv(time_pf, time2_pf, type = "interval2") ~ enrollment_gender_port, port_only, 
                     model = "ph", 
                     weights = NULL,
                     B = c(0,1))
summary(sex_surv_pf_sp)
sex_data <- data.frame(enrollment_gender_port = c(0,1))
rownames(sex_data) <- c('female', 'male')
plot(sex_surv_pf_sp, sex_data)

#fit semi-parametric
diag_baseline(Surv(time_pf, time2_pf, type = "interval2") ~ enrollment_gender_port, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
sex_surv_pf_weibull <- ic_par(Surv(time_pf, time2_pf, type = "interval2") ~ enrollment_gender_port, port_only, model = "ph", dist = "weibull", weights = NULL)
sex_surv_pf_weibull
plot(sex_surv_pf_weibull)

###early pf positivity, no PfTZ
km_pf_pf <- survfit(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont, data = port_only)
summary(km_pf_pf)
plot(km_pf_pf)
pf_surv_pf <- ggsurvplot(km_pf_pf, data = port_only, 
                         risk.table = TRUE, 
                         cumcensor = TRUE, 
                         conf.int = TRUE,
                         fun = "event",
                         palette = c("skyblue2","red3"), 
                         legend.labs = c("Pf-","Pf+"),
                         xlab = "Time (d) since first negative test")
pf_surv_pf
ggsave("figures/surv_pf_pf_incen.png", pf_surv_pf, height = 6, width = 6, dpi = 600)


#semi-parametric model with proportional hazards by pf status
pf_surv_pf_sp <- ic_sp(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont, port_only, 
                       model = "ph", 
                       weights = NULL,
                       B = c(0,1))
summary(pf_surv_pf_sp)
pf_data <- data.frame(pf_early = c(0,1))
rownames(pf_data) <- c("Pf-","Pf+")
plot(pf_surv_pf_sp, pf_data)

#fit semi-parametric
diag_baseline(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
ic_par(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont, port_only, model = "ph", dist = "weibull", weights = NULL)

###early pf positivity, including PfTZ (broader PoRT object)
km_pf_pf <- survfit(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont + cohort, data = port)
summary(km_pf_pf)
plot(km_pf_pf)
pf_surv_pf <- ggsurvplot(km_pf_pf, data = port, 
                          risk.table = TRUE, 
                          cumcensor = TRUE, 
                          conf.int = TRUE,
                          fun = "event",
                          palette = c("skyblue2", "green4", "red3"), 
                          legend.labs = c("Pf-", "PfTZ (Pf+)","Pf+"),
                          xlab = "Time (d) since first negative test")
pf_surv_pf
ggsave("figures/surv_pf_pf_portpftz_incen.png", pf_surv_pf, height = 6, width = 6, dpi = 600)


#semi-parametric model with proportional hazards by pf status, including PfTZ

#first, omit observations where time_pf and time2_pf are missing
port_fullcase <- port[(!is.na(port$time_pf) | !is.na(port$time2_pf)),]

#then, encode semi-parametric
pf_surv_pf_sp <- ic_sp(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont + cohort, port_fullcase, 
                        model = "ph", 
                        weights = NULL,
                        B = c(0,1))
summary(pf_surv_pf_sp)

#fit semi-parametric
diag_baseline(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont + cohort, port_fullcase, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
ic_par(Surv(time_pf, time2_pf, type = "interval2") ~ pf_cont + cohort, port_fullcase, model = "ph", dist = "weibull", weights = NULL)

#plot survival between individuals over and under 16yo
km_over15_pf <- survfit(Surv(time_pf, time2_pf, type = "interval2") ~ over15, data = port_only)
summary(km_over15_pf)
over15_surv_pf <- ggsurvplot(km_over15_pf, data = port_only, 
                          risk.table = TRUE, 
                          cumcensor = TRUE, 
                          conf.int = TRUE,
                          fun = "event",
                          palette = c("skyblue2","red2"), 
                          legend.labs = c("under 16yo","16yo+"),
                          xlab = "Time (d) since first negative test")
over15_surv_pf
ggsave("figures/surv_over15_pf_incen.png", over15_surv_pf, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
over15_surv_pf_sp <- ic_sp(Surv(time_pf, time2_pf, type = "interval2") ~ over15, port_only, 
                        model = "ph", 
                        weights = NULL,
                        B = c(0,1))
summary(over15_surv_pf_sp)
over15_data <- data.frame(over15 = c(0,1))
rownames(over15_data) <- c('<16yo', '>=16yo')
plot(over15_surv_pf_sp, over15_data)

#fit semi-parametric
diag_baseline(Surv(time_pf, time2_pf, type = "interval2") ~ over15, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time_pf, time2_pf, type = "interval2") ~ over15, port_only, model = "ph", dist = "weibull", weights = NULL)

#survival anlysis by wet or dry season 6w prior to enrollment date
km_wetdry_pf <- survfit(Surv(time_pf, time2_pf, type = "interval2") ~ wetdry_6w, data = port_only)
summary(km_wetdry_pf)
plot(km_wetdry_pf)
wetdry_surv_pf <- ggsurvplot(km_wetdry_pf, data = port_only, 
                          risk.table = TRUE, 
                          cumcensor = TRUE, 
                          conf.int = TRUE,
                          fun = "event",
                          palette = c("skyblue2","red2"), 
                          legend.labs = c("Dry","Wet"),
                          xlab = "Time (d) since first negative test")
wetdry_surv_pf
ggsave("figures/surv_wetdry6w_pf_incen.png", wetdry_surv_pf, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
wetdry_surv_pf_sp <- ic_sp(Surv(time_pf, time2_pf, type = "interval2") ~ wetdry_6w, port_only, 
                        model = "ph", 
                        weights = NULL,
                        B = c(0,1))
summary(wetdry_surv_pf_sp)
wetdry_data <- data.frame(wetdry_6w = c(0,1))
rownames(wetdry_data) <- c('Dry', 'Wet')
plot(wetdry_surv_pf_sp, wetdry_data)

#fit semi-parametric
diag_baseline(Surv(time_pf, time2_pf, type = "interval2") ~ wetdry_6w, port_only, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time_pf, time2_pf, type = "interval2") ~ wetdry_6w, port_only, model = "ph", dist = "weibull", weights = NULL)


######################################################
######## TIME TO FIRST PO AND PF RECURRENCE ##########
######################################################

#create flag variable for po or pf identity
port <- port %>%
  mutate(po_flag = 1) %>%
  mutate(pf_flag = 0)

#create new dataset with each observation twice
#carrying time and time2 for first reptition
#time_pf and time2_pf for other
#no fully missing individuals, which are PfTZ individuals who were pf-positive at all screened visits
po_pf_recur <- port[(!is.na(port$time_pf) | !is.na(port$time2_pf)),c("screening_id","time","time2", "po_flag")]
po_pf_recur[(nrow(po_pf_recur)+1):(nrow(po_pf_recur)*2),] <- port[(!is.na(port$time_pf) | !is.na(port$time2_pf)),c("screening_id","time_pf","time2_pf", "pf_flag")]

#plot po and pf recurrence
km_popf <- survfit(Surv(time, time2, type = "interval2") ~ po_flag, data = po_pf_recur)
summary(km_popf)
plot(km_popf)
popf_surv <- ggsurvplot(km_popf, data =  po_pf_recur, 
                             risk.table = TRUE, 
                             cumcensor = TRUE, 
                             conf.int = TRUE,
                             fun = "event",
                             palette = c("red2", "skyblue2"), 
                             legend.labs = c("Pf","Po"),
                             xlab = "Time (d) since first negative test")
popf_surv
ggsave("figures/surv_pftz_popf_incen.png", popf_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
popf_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, 
                           model = "ph", 
                           weights = NULL,
                           B = c(0,1))

summary(popf_surv_sp)
popf_data <- data.frame(po_flag = c(0,1))
rownames(popf_data) <- c("Pf","Po")
plot(popf_surv_sp, popf_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, model = "ph", dist = "exponential", weights = NULL)
ic_par(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, model = "ph", dist = "weibull", weights = NULL)

##time to Po and Pf recurrence without PfTZ

#subset to port observations only
port_only <- port[(port$cohort == "port"),]

#create new dataset with each observation twice
#carrying time and time2 for first reptition
#time_pf and time2_pf for other

po_pf_recur <- port_only[,c("screening_id","time","time2", "po_flag")]
po_pf_recur[(nrow(po_pf_recur)+1):(nrow(po_pf_recur)*2),] <- port_only[,c("screening_id","time_pf","time2_pf", "pf_flag")]

#plot po and pf recurrence
km_popf <- survfit(Surv(time, time2, type = "interval2") ~ po_flag, data = po_pf_recur)
summary(km_popf)
plot(km_popf)
popf_surv <- ggsurvplot(km_popf, data =  po_pf_recur, 
                        risk.table = TRUE, 
                        cumcensor = TRUE, 
                        conf.int = TRUE,
                        fun = "event",
                        palette = c("red2", "skyblue2"), 
                        legend.labs = c("Pf","Po"),
                        xlab = "Time (d) since first negative test")
popf_surv
ggsave("figures/surv_popf_incen.png", popf_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by sex
popf_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, 
                      model = "ph", 
                      weights = NULL,
                      B = c(0,1))

summary(popf_surv_sp)
popf_data <- data.frame(po_flag = c(0,1))
rownames(popf_data) <- c("Pf","Po")
plot(popf_surv_sp, popf_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, model = "ph", dist = "exponential", weights = NULL)
ic_par(Surv(time, time2, type = "interval2") ~ po_flag, po_pf_recur, model = "ph", dist = "weibull", weights = NULL)

######################################################
### LTFU #############################################
######################################################

#LTFU is anyone with time2/time2_pf NA and time/time_pf not 168
#censor is anyone without event observed (including admin censor)
port <- port %>%
  mutate(ltfu = ifelse(is.na(time2) & time < 168, 1, 0)) %>%
  mutate(censor = ifelse (is.na(time2), 1, 0))

###LTFU by demographics

#sex
port %>% tabyl(enrollment_gender_port, ltfu) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages(c("row")) %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")
#not associated
port %>% tabyl(enrollment_gender_port, censor) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")
#not associated

#age
port %>% tabyl(adult, ltfu) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#not associated
port %>% tabyl(adult, censor) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#not associated
port %>% tabyl(over15, ltfu) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#not associated
port %>% tabyl(over15, censor) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#not associated

#early pf
port %>% tabyl(pf_early, ltfu) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#individuals with early pf more likely to be censored
port %>% tabyl(pf_early, censor) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#individuals with early pf somewhat more likely to be censored

#early persistent po
port %>% tabyl(po_cont, ltfu) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#individuals with continuous po positivity at beginning more likely to be censored
#likely because delays study entry. Uncensored individuals are those who had recurrence
port %>% tabyl(po_cont, censor) %>%
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2)
#admin censoring at maximum and admin censoring from persistent parasitemia were driving the earlier difference

######################################################
######## NUMBER OF PO RECURRENCES PER SUBJECT ########
######################################################

####evaluate number of recurrences per patient and timing of positive tests

#create empty dataframe of positive tests
pos <- data.frame(matrix(nrow = 0, ncol = 3))
names(pos) <- c("screening_id", "test_date", "recur_count")

#create days list including screening
days_sn <- c(0, 0, 14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 168)

#switched naming convention here
#port_pftz contains both cohorts
port_pftz <- port
#port only contains port
port <- port[(port$cohort == "port"),]

#iterate through dataset, identify number of recurrences with intervening negative tests
for (i in 1:nrow(port)) {
  print(port$screening_id[i])
  #create vector of po status from enrollment to wk2 - wk24
  screen <- c(port$po_screen[i], port$po_enroll[i], port$po_wk2[i], port$po_wk4[i], port$po_wk6[i], port$po_wk8[i],
              port$po_wk10[i], port$po_wk12[i], port$po_wk14[i], port$po_wk16[i],
              port$po_wk18[i], port$po_wk20[i], port$po_wk22[i], port$po_wk24[i])
  #set flag variable as 0, which flips the when time the individual screens negative (and is at risk for recurrence)
  recur_flag <- 0
  port$recur_count[i] <- 0
  #for each item in the list of screening results. j can be used to index a screening result and the corresponding day count in the days list
  for (j in 1:length(screen)) {
    #if they are positive for the first time (screening or enrollment), record as recurrence 0
    if (port$recur_count[i] == 0  & !is.na(screen[j]) & screen[j] == 1) {
      recur_flag <- 0
      pos[nrow(pos)+1,] <- c(port$screening_id[i], format(as.Date(port$port_sample_date_d0[i]) + days_sn[j],"%Y-%m-%d"), port$recur_count[i])
      port$recur_count[i] <- port$recur_count[i] + 1}
      #add date of test (date of enrollment + corresponding "days" item) and screening_id to positive test data frame    }
    #the if they are negative, flag turns to 1
    if (recur_flag == 0 & !is.na(screen[j]) & screen[j] == 0) {
      recur_flag <- 1}
    if (recur_flag == 1 & !is.na(screen[j]) & screen[j] == 1) {
        recur_flag <- 0
        #add date of test (date of enrollment + corresponding "days" item) and screening_id to positive test data frame
        pos[nrow(pos) + 1,] <- c(port$screening_id[i], format(as.Date(port$port_sample_date_d0[i]) + days_sn[j],"%Y-%m-%d"), port$recur_count[i])
        port$recur_count[i] <- port$recur_count[i] + 1
    }
  }
}

#extract month and month 6w prior to test for plotting
pos <- pos %>%
  mutate(month = as.numeric(format(as.Date(test_date, format = "%Y-%m-%d"), "%m"))) %>%
  mutate(month_6w = as.numeric(format(as.Date(as.Date(test_date) - 42), "%m"))) 

#graph histogram of positive tests by month
month_hist <- ggplot(pos, aes(x=month, fill = recur_count)) + 
  geom_histogram(binwidth = 1, position = "dodge") + 
  scale_x_continuous(name = "Month of test", breaks = seq(1,12,1), limits = c(0.5,12.5)) +
  labs(fill = "Recurrence\nNumber")
month_hist
ggsave("figures/hist_recurbymonth.png", month_hist, height = 6, width = 6, dpi = 600)


month6w_hist <- ggplot(pos, aes(x=month_6w, fill = recur_count)) + 
  geom_histogram(binwidth = 1, position = "dodge") + 
  scale_x_continuous(name = "Month 6w prior to test", breaks = seq(1,12,1), limits = c(0.5,12.5)) +
  labs(fill = "Recurrence\nNumber")
month6w_hist
ggsave("figures/hist_recurbymonth6w.png", month6w_hist, height = 6, width = 6, dpi = 600)

######################################################
######## NUMBER PF RECURRENCES #######################
######################################################

#create empty dataframe of positive tests
pos_pf <- data.frame(matrix(nrow = 0, ncol = 3))
names(pos_pf) <- c("screening_id", "test_date", "recur_count_pf")

#iterate through dataset, identify number of recurrences with intervening negative tests
for (i in 1:nrow(port)) {
  print(port$screening_id[i])
  #create vector of po status from enrollment to wk2 - wk24
  screen <- c(port$pf_screen[i], port$pf_enroll[i], port$pf_wk2[i], port$pf_wk4[i], port$pf_wk6[i], port$pf_wk8[i],
              port$pf_wk10[i], port$pf_wk12[i], port$pf_wk14[i], port$pf_wk16[i],
              port$pf_wk18[i], port$pf_wk20[i], port$pf_wk22[i], port$pf_wk24[i])
  #set flag variable as 0, which flips the when time the individual screens negative (and is at risk for recurrence)
  recur_flag <- 0
  port$recur_count_pf[i] <- 0
  #for each item in the list of screening results. j can be used to index a screening result and the corresponding day count in the days list
  for (j in 1:length(screen)) {
    #if they are positive for the first time (screening or enrollment), record as recurrence 0
    if (port$recur_count_pf[i] == 0  & !is.na(screen[j]) & screen[j] == 1) {
      recur_flag <- 0
      pos_pf[nrow(pos_pf)+1,] <- c(port$screening_id[i], format(as.Date(port$port_sample_date_d0[i]) + days_sn[j],"%Y-%m-%d"), port$recur_count_pf[i])
      port$recur_count_pf[i] <- port$recur_count_pf[i] + 1}
    #add date of test (date of enrollment + corresponding "days" item) and screening_id to positive test data frame    }
    #the if they are negative, flag turns to 1
    if (recur_flag == 0 & !is.na(screen[j]) & screen[j] == 0) {
      recur_flag <- 1}
    if (recur_flag == 1 & !is.na(screen[j]) & screen[j] == 1) {
      recur_flag <- 0
      #add date of test (date of enrollment + corresponding "days" item) and screening_id to positive test data frame
      pos_pf[nrow(pos_pf) + 1,] <- c(port$screening_id[i], format(as.Date(port$port_sample_date_d0[i]) + days_sn[j],"%Y-%m-%d"), port$recur_count_pf[i])
      port$recur_count_pf[i] <- port$recur_count_pf[i] + 1
    }
  }
}

#extract month and month 6w prior to test for plotting
pos_pf <- pos_pf %>%
  mutate(month = as.numeric(format(as.Date(test_date, format = "%Y-%m-%d"), "%m"))) %>%
  mutate(month_6w = as.numeric(format(as.Date(as.Date(test_date) - 42), "%m"))) 

#graph histogram of positive tests by month
month_hist_pf <- ggplot(pos_pf, aes(x=month, fill = recur_count_pf)) + 
  geom_histogram(binwidth = 1, position = "dodge") + 
  scale_x_continuous(name = "Month of test", breaks = seq(1,12,1), limits = c(0.5,12.5)) +
  labs(fill = "Recurrence\nNumber")
month_hist_pf
ggsave("figures/hist_recurbymonth_pf.png", month_hist_pf, height = 6, width = 6, dpi = 600)


month6w_hist_pf <- ggplot(pos_pf, aes(x=month_6w, fill = recur_count_pf)) + 
  geom_histogram(binwidth = 1, position = "dodge") + 
  scale_x_continuous(name = "Month 6w prior to test", breaks = seq(1,12,1), limits = c(0.5,12.5)) +
  labs(fill = "Recurrence\nNumber")
month6w_hist_pf
ggsave("figures/hist_recurbymonth6w_pf.png", month6w_hist_pf, height = 6, width = 6, dpi = 600)



######################################################
############### SECOND PO RECURRENCES  ###############
######################################################

#duplicate port
port_allrecur <- port

#iterate through. For second recurrences, duplicates row and make time and time2 for that second recurrence
for (i in 1:nrow(port_allrecur)) {
  #create vector of po status from enrollment to wk2 - wk24
  screen <- c(port_allrecur$po_enroll[i], port_allrecur$po_wk2[i], port_allrecur$po_wk4[i], port_allrecur$po_wk6[i], port_allrecur$po_wk8[i],
              port_allrecur$po_wk10[i], port_allrecur$po_wk12[i], port_allrecur$po_wk14[i], port_allrecur$po_wk16[i],
              port_allrecur$po_wk18[i], port_allrecur$po_wk20[i], port_allrecur$po_wk22[i], port_allrecur$po_wk24[i])
  neg_flag <- 0
  port_allrecur$recur2_index[i] <- NA
  port_allrecur$count[i] <- "first"
  #for each item in the list of screening results. j can be used to index a screening result and the corresponding day count in the days list
  for (j in 1:length(screen)) {
    #if screening period is beyond first recurrence and second recurrence has NOT already been identified
    if (days[j] > port_allrecur$time2_enroll[i] & !is.na(port_allrecur$time2_enroll[i])  & neg_flag != 2){
      #if screening is negative (but present)
      if (screen[j] == 0 & !is.na(screen[j])){
        port_allrecur$last_neg_test[i] <- days[j]
        #if this is the first negative test after a positive, update first_neg_after_recur variable and create new data entry
        if (neg_flag == 0) {
          #assign first negative after recurrence
          port_allrecur$first_neg_after_recur[i] <- days[j]
          #create new row
          port_allrecur[nrow(port_allrecur) + 1,] <- port_allrecur[i,]
          #label as second recurrence
          port_allrecur$count[nrow(port_allrecur)] <- "second"
          #initialize time variables for second recurrence
          port_allrecur$time[nrow(port_allrecur)] <- NA
          port_allrecur$time_enroll[nrow(port_allrecur)] <- NA
          port_allrecur$time2[nrow(port_allrecur)] <- NA
          port_allrecur$time2_enroll[nrow(port_allrecur)] <- NA
          #create index variable in first recurrence to locate second recurrence row
          port_allrecur$recur2_index[i] <- nrow(port_allrecur)
          #make sure neg_flag reflects negative test
          neg_flag <- 1
        }
        #set time_enroll as last_neg_test recorded for this observation minus date of first recurrence
        port_allrecur$time_enroll[port_allrecur$recur2_index[i]] <- port_allrecur$last_neg_test[i] - port_allrecur$time2_enroll[i]
        #set time as last_neg_test minus the first negative after their first recurrence
        port_allrecur$time[port_allrecur$recur2_index[i]] <- port_allrecur$last_neg_test[i] - port_allrecur$first_neg_after_recur[i]
      }
      #if individual is positive following the first recurrence and has had a negative test
      if (screen[j] == 1 & !is.na(screen[j]) & neg_flag == 1){
        #update time2_enroll of this new row to be the days since last recurrence
        port_allrecur$time2_enroll[port_allrecur$recur2_index[i]] <- days[j] - port_allrecur$time2_enroll[i]
        #update time2 of this new row to be the days since the first negative after last recurrence
        port_allrecur$time2[port_allrecur$recur2_index[i]] <- days[j] - port_allrecur$first_neg_after_recur[i]
        #set negative flag back to 2 so no additional recursion
        neg_flag <- 2
      }
    }
  }
}

#calculate median time to first recurrence
port_allrecur[(port_allrecur$count == "first" & !is.na(port_allrecur$time2)),c("time2")]
median(port_allrecur[(port_allrecur$count == "first" & !is.na(port_allrecur$time2)),c("time2")])
quantile(port_allrecur[(port_allrecur$count == "first" & !is.na(port_allrecur$time2)),c("time2")], probs = seq(0, 1, 0.25), type = 4)

#calculate median time to second recurrence
port_allrecur[(port_allrecur$count == "second" & !is.na(port_allrecur$time2)),c("screening_id", "time2")]
port_allrecur[(port_allrecur$count == "second" & !is.na(port_allrecur$time2)),c("time2")]
median(port_allrecur[(port_allrecur$count == "second" & !is.na(port_allrecur$time2)),c("time2")])
quantile(port_allrecur[(port_allrecur$count == "second" & !is.na(port_allrecur$time2)),c("time2")], probs = seq(0, 1, 0.25), type = 1)

#plot time to recurrences, inlcuding 1st and 2nd
#basic Kaplan-meier survival plot from 1st po neg with interval censoring
km <- survfit(Surv(time, time2, type = "interval2") ~ count, data = port_allrecur)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port_allrecur, 
                         risk.table = TRUE, 
                         conf.int = TRUE,
                         palette = c("skyblue2","red2"), 
                         cumcensor = TRUE,
                         fun = "event",
                         xlab = "Time (d) since first negative test",
                         legend.labs = c("1st Recurrence","2nd Recurrence"))
total_surv
ggsave("figures/surv_total_allrecur_incen.png", total_surv, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by recurrence number
recur_surv_sp <- ic_sp(Surv(time, time2, type = "interval2") ~ count, port_allrecur, 
                        model = "ph", 
                        weights = NULL,
                        B = c(0,1))
summary(recur_surv_sp)

recur_data <- data.frame(count = c("first","second"))
rownames(recur_data) <- c("1st Recurrence","2nd Recurrence")
plot(recur_surv_sp, recur_data)

#fit semi-parametric
diag_baseline(Surv(time, time2, type = "interval2") ~ count, port_allrecur, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time, time2, type = "interval2") ~ count, port_allrecur, model = "ph", dist = "weibull", weights = NULL)



######################################################
######## SECOND PF RECURRENCES #######################
######################################################

#duplicate port
port_allrecur_pf <- port

#iterate through. For second recurrences, duplicates row and make time and time2 for that second recurrence
for (i in 1:nrow(port_allrecur_pf)) {
  #create vector of po status from enrollment to wk2 - wk24
  screen <- c(port_allrecur_pf$pf_enroll[i], port_allrecur_pf$pf_wk2[i], port_allrecur_pf$pf_wk4[i], port_allrecur_pf$pf_wk6[i], port_allrecur_pf$pf_wk8[i],
              port_allrecur_pf$pf_wk10[i], port_allrecur_pf$pf_wk12[i], port_allrecur_pf$pf_wk14[i], port_allrecur_pf$pf_wk16[i],
              port_allrecur_pf$pf_wk18[i], port_allrecur_pf$pf_wk20[i], port_allrecur_pf$pf_wk22[i], port_allrecur_pf$pf_wk24[i])
  neg_flag <- 0
  port_allrecur_pf$recur2_index[i] <- NA
  port_allrecur_pf$count[i] <- "first"
  #for each item in the list of screening results. j can be used to index a screening result and the corresponding day count in the days list
  for (j in 1:length(screen)) {
    #if screening period is beyond first recurrence and second recurrence has NOT already been identified
    if (days[j] > port_allrecur_pf$time2_enroll_pf[i] & !is.na(port_allrecur_pf$time2_enroll_pf[i])  & neg_flag != 2){
      #if screening is negative (but present)
      if (screen[j] == 0 & !is.na(screen[j])){
        port_allrecur_pf$last_neg_test[i] <- days[j]
        #if this is the first negative test after a positive, update first_neg_after_recur variable and create new data entry
        if (neg_flag == 0) {
          #assign first negative after recurrence
          port_allrecur_pf$first_neg_after_recur[i] <- days[j]
          #create new row
          port_allrecur_pf[nrow(port_allrecur_pf) + 1,] <- port_allrecur_pf[i,]
          #label as second recurrence
          port_allrecur_pf$count[nrow(port_allrecur_pf)] <- "second"
          #initialize time variables for second recurrence
          port_allrecur_pf$time_pf[nrow(port_allrecur_pf)] <- NA
          port_allrecur_pf$time_enroll_pf[nrow(port_allrecur_pf)] <- NA
          port_allrecur_pf$time2_pf[nrow(port_allrecur_pf)] <- NA
          port_allrecur_pf$time2_enroll_pf[nrow(port_allrecur_pf)] <- NA
          #create index variable in first recurrence to locate second recurrence row
          port_allrecur_pf$recur2_index[i] <- nrow(port_allrecur_pf)
          #make sure neg_flag reflects negative test
          neg_flag <- 1
        }
        #set time_enroll as last_neg_test recorded for this observation minus date of first recurrence
        port_allrecur_pf$time_enroll_pf[port_allrecur_pf$recur2_index[i]] <- port_allrecur_pf$last_neg_test[i] - port_allrecur_pf$time2_enroll_pf[i]
        print(port_allrecur_pf$screening_id[i])
        print(port_allrecur_pf$time_enroll_pf[port_allrecur_pf$recur2_index[i]])
        #set time as last_neg_test minus the first negative after their first recurrence
        port_allrecur_pf$time_pf[port_allrecur_pf$recur2_index[i]] <- port_allrecur_pf$last_neg_test[i] - port_allrecur_pf$first_neg_after_recur[i]
      }
      #if individual is positive following the first recurrence and has had a negative test
      if (screen[j] == 1 & !is.na(screen[j]) & neg_flag == 1){
        #update time2_enroll of this new row to be the days since last recurrence
        port_allrecur_pf$time2_enroll_pf[port_allrecur_pf$recur2_index[i]] <- days[j] - port_allrecur_pf$time2_enroll_pf[i]
        #update time2 of this new row to be the days since the first negative after last recurrebce
        port_allrecur_pf$time2_pf[port_allrecur_pf$recur2_index[i]] <- days[j] - port_allrecur_pf$first_neg_after_recur[i]
        #set negative flag back to 2 so no additional recursion
        neg_flag <- 2
      }
    }
  }
}

#plot time to recurrences, inlcuding 1st and 2nd
#basic Kaplan-meier survival plot from 1st po neg with interval censoring
km_pf <- survfit(Surv(time_pf, time2_pf, type = "interval2") ~ count, data = port_allrecur_pf)
plot(km_pf)
summary(km_pf)
sum(km_pf$n.event)
total_surv_pf <- ggsurvplot(km_pf, data = port_allrecur_pf, 
                         risk.table = TRUE, 
                         conf.int = TRUE,
                         palette = c("skyblue2","red2"), 
                         cumcensor = TRUE,
                         fun = "event",
                         xlab = "Time (d) since first negative test",
                         legend.labs = c("1st Recurrence","2nd Recurrence"))
total_surv_pf
ggsave("figures/surv_total_allrecur_pf_incen.png", total_surv_pf, height = 6, width = 6, dpi = 600)

#semi-parametric model with proportional hazards by recurrence number
recur_surv_sp_pf <- ic_sp(Surv(time_pf, time2_pf, type = "interval2") ~ count, port_allrecur_pf, 
                       model = "ph", 
                       weights = NULL,
                       B = c(0,1))
summary(recur_surv_sp_pf)

recur_data <- data.frame(count = c("first","second"))
rownames(recur_data) <- c("1st Recurrence","2nd Recurrence")
plot(recur_surv_sp_pf, recur_data)

#fit semi-parametric
diag_baseline(Surv(time_pf, time2_pf, type = "interval2") ~ count, port_allrecur_pf, model = "ph",
              dists = c("exponential", "weibull", 
                        "gamma", "lnorm", "loglogistic"), cols = NULL)
#interval-censored parametric proportional-hazards model with weibull response distribution
#all distributions appear to fit data well for single covariate
ic_par(Surv(time_pf, time2_pf, type = "interval2") ~ count, port_allrecur_pf, model = "ph", dist = "weibull", weights = NULL)


######################################################
############## PFTZ COMPARISON #######################
######################################################

#PfTZ folks sampled in long wet season, so subset to the same for PoRT
port_pftz_comp <- port_pftz %>%
  filter(as.numeric(format(as.Date(port_sample_date_d0), "%m")) >= 3 & as.numeric(format(as.Date(port_sample_date_d0), "%m")) <= 5) %>%
  #truncate data to wk12 (maximum follow-up in PfTZ)
  mutate(timesurv = case_when(timesurv > 84 ~ 84, TRUE ~ timesurv)) %>%
  mutate(weeksurv = timesurv/7) %>%
  mutate(event = case_when(timesurv > 84 ~ 0, TRUE ~ event)) %>%
  mutate(timesurv_pf = case_when(timesurv_pf > 84 ~ 84, TRUE ~ timesurv_pf)) %>%
  mutate(weeksurv_pf = timesurv_pf/7) %>%
  mutate(event_pf = case_when(timesurv_pf > 84 ~ 0, TRUE ~ event_pf))


#plot Po incident cases by being Pf positive at baseline (PfTZ), Po-positive (PoRT where pf_baseline = 0), and both (PoRT where pf_baseline = 1)
km_pftz <- survfit(Surv(weeksurv, event) ~ pf_baseline + cohort, data = port_pftz_comp)
summary(km_pftz)
plot(km_pftz)
pftz_surv <- ggsurvplot(km_pftz, data = port_pftz_comp, 
                      risk.table = TRUE, 
                      cumcensor = TRUE, 
                      conf.int = TRUE,
                      fun = "event",
                      palette = c("#1F78B4", "#E31A1C", "#FDBF6F"), 
                      legend.labs = c("Po+", "Pf+", "Po+ and Pf+"),
                      xlab = "Time (w) since first negative test",
                      legend.title = "",
                      ylim = c(0,1),
                      censor.shape = "X",
                      break.x.by = 2,
                      ylab = "Cumulative Po incidence")
pftz_surv_main <- pftz_surv$plot/pftz_surv$table + plot_layout(heights = c(8,2))
pftz_surv_main
ggsave("figures/surv_po_pftz.png", pftz_surv_main, height = 5, width = 5, dpi = 600)

#plot Pf incident cases by being Pf positive at baseline (PfTZ), Po-positive (PoRT where pf_baseline = 0), and both (PoRT where pf_baseline = 1)
km_pftz_pf <- survfit(Surv(weeksurv_pf, event_pf) ~ pf_baseline + cohort, data = port_pftz_comp)
summary(km_pftz_pf)
plot(km_pftz_pf)
pftz_surv_pf <- ggsurvplot(km_pftz_pf, data = port_pftz_comp, 
                        risk.table = TRUE, 
                        cumcensor = TRUE, 
                        conf.int = TRUE,
                        fun = "event",
                        palette = c("#1F78B4", "#E31A1C", "#FDBF6F"), 
                        legend.labs = c("Po+", "Pf+", "Po+ and Pf+"),
                        xlab = "Time (w) since first negative test",
                        legend.title = "",
                        ylim = c(0,1),
                        censor.shape = "X",
                        break.x.by = 2,
                        ylab = "Cumulative Pf incidence")
pftz_surv_pf_main <- pftz_surv_pf$plot/pftz_surv_pf$table + plot_layout(heights = c(8,2))
pftz_surv_pf_main
ggsave("figures/surv_pf_pftz.png", pftz_surv_pf_main, height = 5, width = 5, dpi = 600)

#test risk difference at final timepoint
#subset to final week
km_pftz_8w <- summary(km_pftz_pf, times = 8)

#risk difference
rd_po_v_pf <- km_pftz_8w$surv[1] - km_pftz_8w$surv[2]
#standard error of RD
diffSE_po_v_pf <- sqrt(km_pftz_8w$std.err[1]^2 + km_pftz_8w$std.err[2]^2)
lCI <- rd_po_v_pf - 1.96 *diffSE_po_v_pf
uCI <- rd_po_v_pf + 1.96 *diffSE_po_v_pf
#z-test
zStat_po_v_pf <- rd_po_v_pf/diffSE_po_v_pf
2*pnorm(abs(zStat_po_v_pf), lower.tail=FALSE)

#risk difference
rd_popf_v_pf <- km_pftz_8w$surv[3] - km_pftz_8w$surv[2]
#standard error of RD
diffSE_popf_v_pf <- sqrt(km_pftz_8w$std.err[3]^2 + km_pftz_8w$std.err[2]^2)
lCI <- rd_popf_v_pf - 1.96 *diffSE_popf_v_pf
uCI <- rd_popf_v_pf + 1.96 *diffSE_popf_v_pf
#z-test
zStat_popf_v_pf <- rd_popf_v_pf/diffSE_popf_v_pf
2*pnorm(abs(zStat_popf_v_pf), lower.tail=FALSE)


######################################################
######## PROCESS PO SPECIES DATA #####################
######################################################


port <- port %>% 
  #poc-positive at screening
  mutate(poc_screening = qpcr_po_speciation_port___2) %>%
  #pow-positive at screening
  mutate(pow_screening = qpcr_po_speciation_port___3) %>%
  #poc-positive at enrollment
  mutate(poc_enroll = qpcr_po_speciation_enroll_port___2) %>%
  #pow-positive at enrollment
  mutate(pow_enroll = qpcr_po_speciation_enroll_port___3) %>%
  #poc-positive at week2
  mutate(poc_wk2 = ifelse(
    is.na(qpcr_po_speciation_wk2_port___2) & is.na(qpcr_po_speciation_wk2fue_port___2), 
    NA,
    ifelse((is.na(qpcr_po_speciation_wk2_port___2) | qpcr_po_speciation_wk2_port___2 == 0) &
            (is.na(qpcr_po_speciation_wk2fue_port___2) | qpcr_po_speciation_wk2fue_port___2 == 0), 
           0, 1)))  %>%
  #pow-positive at week2
  mutate(pow_wk2 = ifelse(
    is.na(qpcr_po_speciation_wk2_port___3) & is.na(qpcr_po_speciation_wk2fue_port___3), 
    NA,
    ifelse((is.na(qpcr_po_speciation_wk2_port___3) | qpcr_po_speciation_wk2_port___3 == 0) &
        (is.na(qpcr_po_speciation_wk2fue_port___3) | qpcr_po_speciation_wk2fue_port___3 == 0), 
        0, 1)))  %>%
  #poc-positive at week4
  mutate(poc_wk4 = ifelse(
    is.na(qpcr_po_speciation_wk4_port___2) & is.na(qpcr_po_speciation_wk4fue_port___2), 
    NA,
    ifelse((is.na(qpcr_po_speciation_wk4_port___2) | qpcr_po_speciation_wk4_port___2 == 0) &
             (is.na(qpcr_po_speciation_wk4fue_port___2) | qpcr_po_speciation_wk4fue_port___2 == 0), 
           0, 1)))  %>%
  #pow-positive at week4
  mutate(pow_wk4 = ifelse(
    is.na(qpcr_po_speciation_wk4_port___3) & is.na(qpcr_po_speciation_wk4fue_port___3), 
    NA,
    ifelse((is.na(qpcr_po_speciation_wk4_port___3) | qpcr_po_speciation_wk4_port___3 == 0) &
        (is.na(qpcr_po_speciation_wk4fue_port___3) | qpcr_po_speciation_wk4fue_port___3 == 0), 
        0, 1)))  %>%
  #poc-positive at week6
  mutate(poc_wk6 = qpcr_po_speciation_wk6_port___2) %>%
  #pow-positive at week6
  mutate(pow_wk6 = qpcr_po_speciation_wk6_port___3) %>%
  #poc-positive at week8
  mutate(poc_wk8 = qpcr_po_speciation_wk8_port___2) %>%
  #pow-positive at week8
  mutate(pow_wk8 = qpcr_po_speciation_wk8_port___3) %>%
  #poc-positive at week10
  mutate(poc_wk10 = qpcr_po_speciation_wk10_port___2) %>%
  #pow-positive at week10
  mutate(pow_wk10 = qpcr_po_speciation_wk10_port___3) %>%
  #poc-positive at week12
  mutate(poc_wk12 = qpcr_po_speciation_wk12_port___2) %>%
  #pow-positive at week10
  mutate(pow_wk12 = qpcr_po_speciation_wk12_port___3) %>%
  #poc-positive at week14
  mutate(poc_wk14 = qpcr_po_speciation_wk14_port___2) %>%
  #pow-positive at week14
  mutate(pow_wk14 = qpcr_po_speciation_wk14_port___3) %>%
  #poc-positive at week16
  mutate(poc_wk16 = qpcr_po_speciation_wk16_port___2) %>%
  #pow-positive at week16
  mutate(pow_wk16 = qpcr_po_speciation_wk16_port___3) %>%
  #poc-positive at week18
  mutate(poc_wk18 = qpcr_po_speciation_wk18_port___2) %>%
  #pow-positive at week18
  mutate(pow_wk18 = qpcr_po_speciation_wk18_port___3) %>%
  #poc-positive at week20
  mutate(poc_wk20 = qpcr_po_speciation_wk20_port___2) %>%
  #pow-positive at week20
  mutate(pow_wk20 = qpcr_po_speciation_wk20_port___3) %>%
  #poc-positive at week22
  mutate(poc_wk22 = qpcr_po_speciation_wk22_port___2) %>%
  #pow-positive at week10
  mutate(pow_wk22 = qpcr_po_speciation_wk22_port___3) %>%
  #poc-positive at week24
  mutate(poc_wk24 = qpcr_po_speciation_wk24_port___2) %>%
  #pow-positive at week10
  mutate(pow_wk24 = qpcr_po_speciation_wk24_port___3)

#drop unnecessary variables
port <- port %>%
  select(-contains(c("qpcr_po_speciation_")))

#create empty time time variables in port
port <- port %>%
  mutate(time_poc = NA) %>%
  mutate(time2_poc = NA) %>%
  mutate(time_enroll_poc = NA) %>%
  mutate(time2_enroll_poc = NA) %>%
  mutate(t_1stpocneg = NA) %>%
  mutate(time_pow = NA) %>%
  mutate(time2_pow = NA) %>%
  mutate(time_enroll_pow = NA) %>%
  mutate(time2_enroll_pow = NA) %>%
  mutate(t_1stpowneg = NA)

#iterate through all individuals in port in order, adding intervals for po and pf positivity following their first negative
for (i in 1:nrow(port)) {
  #create vector of poc status from day0 - wk24
  screen <- c(port$poc_enroll[i], port$poc_wk2[i], port$poc_wk4[i], port$poc_wk6[i], port$poc_wk8[i],
                port$poc_wk10[i], port$poc_wk12[i], port$poc_wk14[i], port$poc_wk16[i],
                port$poc_wk18[i], port$poc_wk20[i], port$poc_wk22[i], port$poc_wk24[i])
  #speciation data doesn't show missingness, so must overwrite from po screening
  screen_po <- c(port$po_enroll[i], port$po_wk2[i], port$po_wk4[i], port$po_wk6[i], port$po_wk8[i],
              port$po_wk10[i], port$po_wk12[i], port$po_wk14[i], port$po_wk16[i],
              port$po_wk18[i], port$po_wk20[i], port$po_wk22[i], port$po_wk24[i])
  #set flag variable as 0, which flips the first time the individual screens negative (and is at risk for recurrence)
  neg_flag = 0
  #for each item in the list of screening results. j can be used to index a screening result and the corresponding day count in the days list
  for (j in 1:length(screen)) {
    #first, overwrite missing tests into Poc screening (from Po screening)
    if (is.na(screen_po[j])) {
      screen[j] <- NA}
    #the first time they are negative, flag turns to 1
    if (neg_flag == 0 & !is.na(screen[j]) & screen[j] == 0) {
      neg_flag <- 1
      port$t_1stpocneg[i] <- days[j]}
    if (neg_flag == 1) {
      #if a given day is po-negative, update, time to day[j] to reflect most recent negative test
      if (!is.na(screen[j]) & screen[j] == 0) {
        port$time_enroll_poc[i] <- days[j]
        port$time_poc[i] <- days[j] - port$t_1stpocneg[i]}
      #if a given day is po-positive, update time2 to day[j] to reflect first recurrence, 
      #then break for loop
      if (!is.na(screen[j]) & screen[j] == 1) {
        port$time2_enroll_poc[i] <- days[j]
        port$time2_poc[i] <- days[j] - port$t_1stpocneg[i]
        break
      }
    }
  }
  #create vector of pow status from enrollment to wk2 - wk24
  screen_pow <- c(port$pow_enroll[i], port$pow_wk2[i], port$pow_wk4[i], port$pow_wk6[i], port$pow_wk8[i],
                   port$pow_wk10[i], port$pow_wk12[i], port$pow_wk14[i], port$pow_wk16[i],
                   port$pow_wk18[i], port$pow_wk20[i], port$pow_wk22[i], port$pow_wk24[i])
  #set flag variable as 0, which flips the first time the individual screens negative (and is at risk for recurrence)
  neg_flag_pow = 0
  for (j in 1:length(screen_pow)) {
    #first, overwrite missing tests into Pow screening (from Po screening)
    if (is.na(screen_po[j])) {
      screen_pow[j] <- NA}
    #the first time they are negative, flag turns to 1
    if (neg_flag_pow == 0 & !is.na(screen_pow[j]) & screen_pow[j] == 0) {
      neg_flag_pow <- 1
      port$t_1stpowneg[i] <- days[j]}
    if (neg_flag_pow == 1) {
      #if a given day is po-negative, update, time to day[j] to reflect most recent negative test
      if (!is.na(screen_pow[j]) & screen_pow[j] == 0) {
        port$time_enroll_pow[i] <- days[j]
        port$time_pow[i] <- days[j] - port$t_1stpowneg[i]}
      #if a given day is po-positive, update time2 to day[j] to reflect first recurrence, 
      #then break for loop
      if (!is.na(screen_pow[j]) & screen_pow[j] == 1) {
        port$time2_enroll_pow[i] <- days[j]
        port$time2_pow[i] <- days[j] - port$t_1stpowneg[i]
        break
      }
    }
  }
}

#evaluate potential relapses
print(port[port$poc_enroll == 1 & !is.na(port$time2_poc),])
print(port[port$pow_enroll == 1 & !is.na(port$time2_pow),])

#conventional survival coding
port <- port %>%
  mutate(timesurv_poc = case_when(!is.na(time2_poc) ~ time2_poc, is.na(time2_poc) ~ time_poc)) %>%
  mutate(weeksurv_poc = timesurv_poc/7) %>%
  mutate(event_poc = case_when(!is.na(time2_poc) ~ 1, is.na(time2_poc) ~ 0)) %>%
  mutate(timesurv_pow = case_when(!is.na(time2_pow) ~ time2_pow, is.na(time2_pow) ~ time_pow)) %>%
  mutate(event_pow = case_when(!is.na(time2_pow) ~ 1, is.na(time2_pow) ~ 0)) %>%
  mutate(weeksurv_pow = timesurv_pow/7)


#evaluate survival for species-specific survival
#basic Kaplan-meier survival plot to Poc incidence from 1st poc neg with interval censoring
km <- survfit(Surv(time_poc, time2_poc, type = "interval2") ~ 1, data = port)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port, 
                         risk.table = TRUE, 
                         color = "red2", 
                         conf.int = TRUE,
                         cumcensor = TRUE, 
                         fun = "event",
                         ylim = c(0,1),
                         xlab = "Time (d) since first negative test")
total_surv
ggsave("figures/surv_total_poc_incen.png", total_surv, height = 6, width = 6, dpi = 600)

#basic Kaplan-meier survival plot to Pow incidence from 1st poc neg with interval censoring
km <- survfit(Surv(time_pow, time2_pow, type = "interval2") ~ 1, data = port)
plot(km)
summary(km)
sum(km$n.event)
total_surv_pow <- ggsurvplot(km, data = port, 
                         risk.table = TRUE, 
                         color = "red2", 
                         conf.int = TRUE,
                         cumcensor = TRUE, 
                         fun = "event",
                         ylim = c(0,1),
                         xlab = "Time (d) since first negative test")
total_surv_pow
ggsave("figures/surv_total_pow_incen.png", total_surv_pow, height = 6, width = 6, dpi = 600)

#copy species-time into same dataset
#create flag variable for po or pf identity
port <- port %>%
  mutate(poc_flag = 1) %>%
  mutate(pow_flag = 0)

#create new dataset with each observation twice
#carrying time and time2 for first repetition
#time_pf and time2_pf for other
#no fully missing individuals, which are PfTZ individuals who were pf-positive at all screened visits
species <- port[,c("screening_id","weeksurv_poc","event_poc", "poc_flag")]
species[(nrow(species)+1):(nrow(species)*2),] <- port[,c("screening_id","weeksurv_pow","event_pow", "pow_flag")]
print(species)

#poc/pow survival
km <- survfit(Surv(weeksurv_poc, event_poc) ~ poc_flag, data = species)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = species, 
                             #risk.table = TRUE, 
                             palette = c("skyblue2", "green2"),
                             legend.labs = c("Pow","Poc"),
                             legend.title = "",
                             #conf.int = TRUE,
                             fun = "event",
                             ylim = c(0,0.25),
                             censor.shape = "",
                             break.x.by = 2,
                             xlab = "Time (d) since first negative test",
                             ylab = "Cumulative incidence")
total_surv
ggsave("figures/surv_total_pocpow.png", total_surv, height = 6, width = 6, dpi = 600)


######################################################
#### LONG DATA FORM FOR TIME-VARYING COVARIATES ######
######################################################


#Reorient to long data and use Surv() for broad analysis by biweekly intervals
#long dataset should have a three rows for each two week window per participant, one showing po positivity or missingness, one showing pf positivity or missingness, one showing treatment or missingness
#time (open) and time2 (closed) should define interval
# event will be 0 =right censored, 1 =event at time, 2 =left censored, 3 =interval censored
#should also include all additional variables for each individual (positivity at each week, sex, age, etc)
#initial survival analysis of po recurrence will only use the po rows
#how will this account for multiple recurrences?

#duplicate observations, one per po observation, by po bi-weekly screening variables
port_po_long_allvisits <- port %>%
  pivot_longer(c(po_enroll, po_wk2, po_wk4, po_wk6, po_wk8, po_wk10, po_wk12, po_wk14, po_wk16, po_wk18, po_wk20, po_wk22, po_wk24,
                 pf_enroll, pf_wk2, pf_wk4, pf_wk6, pf_wk8, pf_wk10, pf_wk12, pf_wk14, pf_wk16, pf_wk18, pf_wk20, pf_wk22, pf_wk24),
    names_to = c(".value", "interval"), names_sep = "_") %>%
  #create time and time2 variables corresponding to each interval
  mutate(in_ = case_when(interval == "enroll" ~ -14, interval == "wk2" ~ 0, interval == "wk4" ~ 14, 
                         interval == "wk6" ~ 28, interval == "wk8" ~ 42, 
                         interval == "wk10" ~ 56, interval == "wk12" ~ 70,
                         interval == "wk14" ~ 84, interval == "wk16" ~ 98,
                         interval == "wk18" ~ 112, interval == "wk20" ~ 126,
                         interval == "wk22" ~ 140, interval == "wk24" ~ 154)) %>%
  mutate(out_ = in_ + 14) %>%
  #add in treatment data, reflecting treatment at start of interval
  mutate(treat_in = case_when(interval == "enroll" ~ NA, interval == "wk2" ~ port_sample_treatyn_do, interval == "wk4" ~ port_sample_treatyn_w2, 
                              interval == "wk6" ~ port_sample_treatyn_w4, interval == "wk8" ~ port_sample_treatyn_w6, 
                              interval == "wk10" ~ port_sample_treatyn_w8, interval == "wk12" ~ port_sample_treatyn_w10,
                              interval == "wk14" ~ port_sample_treatyn_w12, interval == "wk16" ~ port_sample_treatyn_w14,
                              interval == "wk18" ~ port_sample_treatyn_w16, interval == "wk20" ~ port_sample_treatyn_w18,
                              interval == "wk22" ~ port_sample_treatyn_w20, interval == "wk24" ~ port_sample_treatyn_w22)) %>%
  #add in treatment at end of interval
  mutate(treat_out = case_when(interval == "enroll" ~ port_sample_treatyn_do, interval == "wk2" ~ port_sample_treatyn_w2, interval == "wk4" ~ port_sample_treatyn_w4, 
                                                 interval == "wk6" ~ port_sample_treatyn_w6, interval == "wk8" ~ port_sample_treatyn_w8, 
                                                 interval == "wk10" ~ port_sample_treatyn_w10, interval == "wk12" ~ port_sample_treatyn_w12,
                                                 interval == "wk14" ~ port_sample_treatyn_w14, interval == "wk16" ~ port_sample_treatyn_w16,
                                                 interval == "wk18" ~ port_sample_treatyn_w18, interval == "wk20" ~ port_sample_treatyn_w20,
                                                 interval == "wk22" ~ port_sample_treatyn_w22, interval == "wk24" ~ port_sample_treatyn_w24)) %>%
  #add in symptomatic at end of interval
  mutate(symp_out = case_when(
    interval == "enroll" ~ port_sample_sympyn_d0, interval == "wk2" ~ port_sample_sympyn_w2, interval == "wk4" ~ port_sample_sympyn_w4, 
    interval == "wk6" ~ port_sample_sympyn_w6, interval == "wk8" ~ port_sample_sympyn_w8, 
    interval == "wk10" ~ port_sample_sympyn_w10, interval == "wk12" ~ port_sample_sympyn_w12,
    interval == "wk14" ~ port_sample_sympyn_w14, interval == "wk16" ~ port_sample_sympyn_w16,
    interval == "wk18" ~ port_sample_sympyn_w18, interval == "wk20" ~ port_sample_sympyn_w20,
    interval == "wk22" ~ port_sample_sympyn_w22, interval == "wk24" ~ port_sample_sympyn_w24
  )) %>%
  #group observations by subject
  group_by(screening_id) %>%
  #for missed visits, carry forward missing data
  fill(po, .direction = c("down")) %>%
  fill(pf, .direction = c("down")) %>%
  #if enrollment data is still missing, draw from screening
  mutate(po = ifelse(is.na(po) & interval == "enroll", po_screening, po)) %>%
  mutate(pf = ifelse(is.na(pf) & interval == "enroll", pf_screening, pf)) %>%
  #encode lagged pf-positivity, reflecting positive at start of interval (includes forward-carrying)
  mutate(pf_in = ifelse(interval == "enroll", pf_screening, lag(pf, n = 1, NA))) %>%
  #encode pf positive within past two visits
  mutate(pf_4wk = ifelse(pf_in == 1 | (!is.na(lag(pf_in, 1, NA)) & lag(pf_in, 1, NA) == 1), 1, 0)) %>%
  #individuals with missing treatment data were not treated
  mutate(treat_in = ifelse(is.na(treat_in), 0, treat_in)) %>%
  mutate(treat_out = ifelse(is.na(treat_out), 0, treat_out)) %>%
  #encode individuals treated in last month
  mutate(treat_4wk = ifelse(treat_in == 1 | (!is.na(lag(treat_in, 1, NA)) & lag(treat_in, 1, NA) == 1), 1, 0)) %>%
  #encode individuals treated in last 2 months
  mutate(treat_8wk = ifelse(treat_in == 1 | (!is.na(lag(treat_in, 1, NA)) & lag(treat_in, 1, NA) == 1) | (!is.na(lag(treat_in, 2, NA)) & lag(treat_in, 2, NA) == 1) | (!is.na(lag(treat_in, 3, NA)) & lag(treat_in, 3, NA) == 1), 1, 0)) %>%
  #encode individuals treated in last 2 months
  mutate(treat_10wk = ifelse(treat_in == 1 | (!is.na(lag(treat_in, 1, NA)) & lag(treat_in, 1, NA) == 1) | (!is.na(lag(treat_in, 2, NA)) & lag(treat_in, 2, NA) == 1) | (!is.na(lag(treat_in, 3, NA)) & lag(treat_in, 4, NA) == 1) | (!is.na(lag(treat_in, 3, NA)) & lag(treat_in, 4, NA) == 1), 1, 0)) %>%
  #encode reasoning for treat_out
  mutate(treat_out_reason = case_when(
    #blank if not treated at start of interval
    treat_out == 0 ~ "NA",
    #if symptomatic at the timepoint
    treat_out == 1 & symp_out == 1 ~ "symp",
    #pf persistent treated if persistently positive at prior two visits
    treat_out == 1 & pf_in == 1 & lag(pf_in, 1, NA) == 1 ~ "pf-pers",
    #po persistent treated if po-positive at prior two visits
    treat_out == 1 & lag(po, 1, NA) == 1 & lag(po, 2, NA) == 1~ "po-pers",
    #pf treated if pf-positive at previous visit or current visit
    treat_out == 1 & pf_in == 1 ~ "pf",
    treat_out == 1 & lead(pf_in, 1, NA) == 1 ~ "pf",
    #po treated if po-positive at last visit or current visit
    treat_out == 1 & lag(po, 1, NA) == 1 ~ "po",
    treat_out == 1 & po == 1 ~ "po",
    #otherwise, unknown
    treat_out == 1 ~ "unknown")) %>%
  #treat_in_reason is treat_out lagged
  mutate(treat_in_reason = lag(treat_out_reason, 1, NA)) %>%
  #encode treatment for symptomatic, persistent, or single pf
  mutate(treat_pf = case_when(treat_in_reason == "pf" | treat_in_reason == "pf-pers" | treat_in_reason == "symp" ~ 1, TRUE ~ 0)) %>%
  #pf_treatment in last 3 mo
  mutate(treat_pf_3mo = case_when(
    treat_pf == 1 | lag(treat_pf, 1, NA) |
    lag(treat_pf, 2, NA) | lag(treat_pf, 3, NA) |
    lag(treat_pf, 3, NA) | lag(treat_pf, 4, NA) ~ 1, TRUE ~ 0)) %>%
  #check for missed po treatment
  mutate(po_treat_miss = case_when(po == 1 & treat_out == 0 & lead(treat_out, 1, NA) == 0 ~ 1  & lag(treat_out, 1, NA) == 0 & out_ > 28, TRUE ~ 0)) %>%
  #encode LTFU as a drop variable representing the final visit in which they were screened
  mutate(drop = ifelse(!is.na(t_visit) & out_ == t_visit, 1, 0)) %>%
  #remove intervals after t_visit (LTFU) indexed from enrollment (not first negative test)
  filter(out_ <= t_visit | is.na(t_visit)) %>%
  #calculate month 6w prior to out_
  mutate(date_6w_lag = as.Date(as.Date(port_sample_date_d0) + out_ - 42)) %>%
  mutate(month_vary = as.numeric(format(date_6w_lag, "%m"))) %>%
  #determine season 6w prior to each test
  mutate(season_6w = case_when(
    month_vary >= 10 & month_vary <= 12 ~ "shortwet",
    month_vary >= 3 & month_vary <= 5 ~ "longwet",
    TRUE ~ "dry")) %>%
  mutate(long_wet_6w = case_when(
    month_vary >= 3 & month_vary <= 5 ~ 1, 
    TRUE ~ 0)) %>%
  mutate(short_wet_6w = case_when(
    month_vary >= 10 & month_vary <= 12 ~ 1, 
    TRUE ~ 0)) %>%
  #calculate the month 6w prior to interval onset
  mutate(month_onset = as.numeric(format(as.Date(as.Date(port_sample_date_d0) + in_ - 42), "%m"))) %>%
  #determine season 6w prior to interval onset
  mutate(season_6w_onset = case_when(
    month_onset >= 10 & month_onset <= 12 ~ "shortwet",
    month_onset >= 3 & month_onset <= 5 ~ "longwet",
    TRUE ~ "dry")) %>%
  mutate(long_wet_6w_onset = case_when(
    month_onset >= 3 & month_onset <= 5 ~ 1, 
    TRUE ~ 0)) %>%
  mutate(short_wet_6w_onset = case_when(
    month_onset >= 10 & month_onset <= 12 ~ 1, 
    TRUE ~ 0)) %>%
  #calculate month no lag prior to out_
  mutate(month_nolag = as.numeric(format(as.Date(as.Date(port_sample_date_d0) + out_), "%m"))) %>%
  #determine season no lag prior to each test
  mutate(season_nolag = case_when(
    month_nolag >= 10 & month_nolag <= 12 ~ "short wet",
    month_nolag >= 3 & month_nolag <= 5 ~ "long wet",
    TRUE ~ "dry")) %>%
  mutate(long_wet_nolag = case_when(
    month_nolag >= 3 & month_nolag <= 5 ~ 1, 
    TRUE ~ 0)) %>%
  mutate(short_wet_nolag = case_when(
    month_nolag >= 10 & month_nolag <= 12 ~ 1, 
    TRUE ~ 0)) %>%
  #subtract 1st poneg from each interval so they reflect time since becoming at risk
  mutate(in_ = in_ - t_1stponeg) %>%
  mutate(out_ = out_ - t_1stponeg) %>%
  #recode days as weeks
  mutate(in_wk = in_/7) %>%
  mutate(out_wk = out_/7)


#examine treatment reason
print(port_po_long_allvisits[port_po_long_allvisits$treat_out_reason == "unknown",c("screening_id", "interval", "treat_out_reason")])
print(port_po_long_allvisits[port_po_long_allvisits$treat_out_reason == "symp",c("screening_id", "interval", "treat_out_reason", "po", "pf")])
print(port_po_long_allvisits[port_po_long_allvisits$treat_out_reason == "po" | port_po_long_allvisits$treat_out_reason == "po-pers",c("screening_id", "interval")], n = 100)
print(port_po_long_allvisits[port_po_long_allvisits$po_treat_miss == 1,c("screening_id", "interval", "po")], n = 100)


#recode inconsistent use of treatment
port_po_long_allvisits <- port_po_long_allvisits %>%
  mutate(treat_out = case_when(treat_out_reason == "unknown" ~ 0, TRUE ~ treat_out),
         treat_out_reason = case_when(treat_out_reason == "unknown" ~ "NA", TRUE ~ treat_out_reason))

#examine overall breakdown of treatment
port_po_long_allvisits %>% tabyl(treat_out_reason) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#remove rows with negative time
port_po_long <- port_po_long_allvisits %>%
filter(in_ >= 0) %>%
  filter(out_ > 0)

#determine individuals positive for pf during study follow-up
port_pf <- port_po_long %>%
  group_by(screening_id) %>%
  summarise(pf = sum(pf)) %>%
  mutate(pf_study = case_when(pf > 0 ~ 1, pf == 0 ~ 0)) 

port_pf <- subset(port_pf, select = -c(pf))

port_pf %>% tabyl(pf_study) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#add flag for individuals positive for Pf during follow-up
port_po_long <- left_join(port_po_long, port_pf, by = "screening_id" )

#determine individuals treated during study follow-up
port_treat <- port_po_long %>%
  summarise(treat_in = sum(treat_in)) %>%
  mutate(treat_study = case_when(treat_in > 0 ~ 1, treat_in == 0 ~ 0))

port_treat <- subset(port_treat, select = -c(treat_in))

port_treat %>% tabyl(treat_study) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#determine treatment, attach to broader dataset
port <- left_join(port, port_treat, by = "screening_id" )
port_po_long <- left_join(port_po_long, port_treat, by = "screening_id" )

#determine individuals treated for pf
port_treatpf <- port_po_long %>%
  summarise(treat_pf = sum(treat_pf)) %>%
  mutate(treat_pf_study = case_when(treat_pf > 0 ~ 1, treat_pf == 0 ~ 0))

port_treatpf <- subset(port_treatpf, select = -c(treat_pf))

sum(port_treatpf$treat_pf_study)

#attach pf treatment to broader dataset
port_po_long <- left_join(port_po_long, port_treatpf, by = "screening_id" )

#determine individuals symptomatic during study
port_symp <- port_po_long_allvisits %>%
  group_by(screening_id) %>%
  summarise(symp_total = sum(symp_out)) %>%
  mutate(symp_study = case_when(symp_total > 0 ~ 1, symp_total == 0 ~ 0)) 

port_symp %>% tabyl(symp_total) %>%
  #Checked by hand; NAs assume no symptoms
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#determine season at study onset
port_po_long %>% filter(in_ == 0) %>% tabyl(season_6w_onset) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#determine season across all tests
port_po_long %>% tabyl(season_6w, po) %>%
adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

port_po_long %>% filter(season_6w == "dry") %>%
  tabyl(month_vary, po) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#limit to first recurrence
port_po_long_1st <- port_po_long %>%
  #remove individuals after their first po recurrene
  filter(out_ <= time2 | is.na(time2))

######plot recurrence risk alongside treatment in study

#basic Kaplan-meier survival plot from 1st po neg with conventional coding
km <- survfit(Surv(weeksurv, event, type = "right") ~ 1, data = port_only)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port_only, 
                         risk.table = TRUE, 
                         color = "#1F78B4", 
                         cumcensor = TRUE, 
                         fun = "event",
                         ylab = "Cumulative incidence",
                         xlab = "Time (w) since first negative test",
                         censor.shape = "X",
                         break.x.by = 2,
                         ylim = c(0,1),
                         legend.title = "")

#encode treatment by week
#could switch to port_po_long_1st if we only want to look at treatment before the first recurrence
port_treat_sum <- port_po_long_allvisits %>%
  filter(treat_out_reason != "NA") %>%
  group_by(out_wk) %>%
  count(treat_out_reason) %>%
  mutate(reason_factor = as.factor(treat_out_reason)) %>%
  mutate(reason_factor = fct_relevel(reason_factor, 'symp','pf-pers', 'pf', 'po-pers', 'po'))

#cut individuals treated prior to entering risk set
port_treat_sum_d0 <- port_treat_sum %>%
  filter(out_wk >= 0)


#plot with bar for treated
total_surv_addbar <-  total_surv$plot + 
  #add bar showing number treated by week, but divide by 10 for matching of the two axes
  geom_bar(aes(x=out_wk, y = n/6, fill = reason_factor), stat = 'identity', alpha = 0.3, data = port_treat_sum_d0) +
  #add manual y scale with additional axis for treatment, scaled by 10 relative to the incidence
  scale_y_continuous(breaks = seq(0,1,by = 0.25), limits = c(0, 1), sec.axis = sec_axis(~.*6, name = '# treated', breaks = seq(0,10,1))) +
  #manual fill scale
  scale_fill_manual(values = c("pf" = "#E31A1C", "po" = "#6A3D9A", "pf-pers" = "#FF7F00", "po-pers" = "#B15928", "symp" = "#FFFF99"), labels = c("pf" = "Late Pf", "po" = "Late Po", "pf-pers" = "Pers. Pf", "po-pers" = "Pers. Po", "symp" = "Symp. Pf"))
total_surv_addbar <- move_layers(total_surv_addbar, "GeomBar", position = "bottom")
total_surv_addbar
total_surv_plot <- total_surv_addbar/total_surv$table/total_surv$ncensor.plot + plot_layout(heights=c(8,2, 2))
total_surv_plot
ggsave("figures/surv_total_no-treated.png", total_surv_plot, height = 9, width = 7, dpi = 600)

#plot treatment bar graph on its own
treat_bar <-  ggplot(data = port_treat_sum) + 
  #add bar showing number treated by week, but divide by 10 for matching of the two axes
  geom_bar(aes(x=out_wk, y = n, fill = reason_factor), stat = 'identity', alpha = 0.8) +
  #add manual x and y scale
  scale_y_continuous(breaks = seq(0,6,by = 1)) +
  scale_x_continuous(breaks = seq(-8,24,by = 2)) +
  #manual fill scale
  scale_fill_manual(values = c("pf" = "#E31A1C", "po" = "#6A3D9A", "pf-pers" = "#FF7F00", "po-pers" = "#B15928", "symp" = "#FFFF99"), labels = c("pf" = "Late Pf", "po" = "Late Po", "pf-pers" = "Pers. Pf", "po-pers" = "Pers. Po", "symp" = "Symp. Pf")) +
  xlab("Time since first negative test (w)") +
  ylab("# treated") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14),
        legend.position = "top", 
        panel.grid = element_blank())
treat_bar
ggsave("figures/treated_histo.png", treat_bar, height = 5, width = 6, dpi = 600)

#determine treatment prior to 1st recurrence
port_treat1st <- port_po_long_1st %>%
  summarise(treat_in = sum(treat_in)) %>%
  mutate(treat_1st = case_when(treat_in > 0 ~ 1, treat_in == 0 ~ 0))

port_treat1st <- subset(port_treat1st, select = -c(treat_in))

#attach treatment prior to 1st recur to broader dataset
port <- left_join(port, port_treat1st, by = "screening_id" )

port_po_long_1st <- left_join(port_po_long_1st, port_treat1st, by = "screening_id" )

#plot survival by any treatment in study prior to first po occurrence
km_treat1st <- survfit(Surv(weeksurv, event) ~ treat_1st, data = port)
summary(km_treat1st)
treat1st_surv <- ggsurvplot(km_treat1st, data = port, 
                              risk.table = TRUE, 
                              palette = c("#1F78B4","#E31A1C"), 
                              legend.labs = c("Untreated","Treated"),
                              #cumcensor = TRUE, 
                              conf.int = TRUE,
                              fun = "event",
                              ylab = "Cumulative incidence",
                              xlab = "Time (w) since first negative test",
                              censor.shape = "",
                              break.x.by = 2,
                              ylim = c(0,1),
                              legend.title = "")
treat1st_surv_main <- treat1st_surv$plot/treat1st_surv$table + plot_layout(heights = c(8,2))
treat1st_surv_main
ggsave("figures/surv_treat1st.png", treat1st_surv_main, height = 5, width = 5, dpi = 600)

#test risk difference at final timepoint
#subset to final week
treat1st_24w <- summary(km_treat1st, times = 24)
#risk difference
rd_treat1st <- treat1st_24w$surv[1] - treat1st_24w$surv[2]
#standard error of RD
diffSE_treat1st <- sqrt(treat1st_24w$std.err[1]^2 + treat1st_24w$std.err[2]^2)
lci <- rd_treat1st - 1.96*diffSE_treat1st
lci
uci <- rd_treat1st + 1.96*diffSE_treat1st
uci 
#z-test
zStat_treat1st <- rd_treat1st/diffSE_treat1st
2*pnorm(abs(zStat_treat1st), lower.tail=FALSE)

#investigate timing of treatment prior to 1st recurrences
quantile(deframe(port_po_long_1st[port_po_long_1st$treat_in == 1,c("in_wk")]), probs = seq(0, 1, 0.25))

#overall survival
km <- survfit(Surv(in_, out_, po) ~ 1, data = port_po_long_1st)
plot(km)
summary(km)
sum(km$n.event)
total_surv <- ggsurvplot(km, data = port_po_long_1st, 
                         risk.table = TRUE, 
                         palette = "red2", 
                         conf.int = TRUE,
                         cumcensor = TRUE, 
                         fun = "event",
                         ylim = c(0,1),
                         xlab = "Time (d) since first negative test")
total_surv
ggsave("figures/surv_long.png", total_surv, height = 6, width = 6, dpi = 600)

#survival by pf positivity at beginning of interval (or most recent test)
km_pf <- survfit(Surv(in_wk, out_wk, po) ~ pf_in, data = port_po_long_1st)
plot(km_pf)
summary(km_pf)
sum(km_pf$n.event)
pf_surv <- ggsurvplot(km_pf, data = port_po_long_1st, 
                    risk.table = TRUE, 
                    palette = c("#1F78B4","#E31A1C"), 
                    legend.labs = c("Pf-", "Pf+ at last visit"),
                    legend.title = "",
                    conf.int = TRUE,
                    fun = "event",
                    ylim = c(0,1),
                    censor.shape = "",
                    break.x.by = 2,
                    xlab = "Time (w) since first negative test",
                    ylab = "Cumulative incidence")
pf_surv_main <- pf_surv$plot/pf_surv$table + plot_layout(heights = c(8,2))
pf_surv_main
ggsave("figures/surv_long_pfin.png", pf_surv_main, height = 5, width = 5, dpi = 600)

coxph_pf <- coxph(Surv(in_, out_, po) ~ pf_in, data = port_po_long_1st)
summary(coxph_pf)

#survival by pf positivity at last two visits
km_pf <- survfit(Surv(in_wk, out_wk, po) ~ pf_4wk, data = port_po_long_1st)
plot(km_pf)
summary(km_pf)
sum(km_pf$n.event)
pf_surv <- ggsurvplot(km_pf, data = port_po_long_1st, 
                      risk.table = TRUE, 
                      palette = c("#1F78B4","#E31A1C"), 
                      legend.labs = c("Pf-", "Pf+"),
                      legend.title = "",
                      conf.int = TRUE,
                      fun = "event",
                      ylim = c(0,1),
                      censor.shape = "",
                      break.x.by = 2,
                      xlab = "Time (w) since first negative test",
                      ylab = "Cumulative incidence")
pf_surv_main <- pf_surv$plot/pf_surv$table + plot_layout(heights = c(6,2))
pf_surv_main
ggsave("figures/surv_long_pf4wk.png", pf_surv_main, height = 4, width = 5, dpi = 600)

#test risk difference at final timepoint
pf_24w <- summary(km_pf, times = 24)
rd_pf <- pf_24w$surv[1] - pf_24w$surv[2]
diffSE <- sqrt(pf_24w$std.err[1]^2 + pf_24w$std.err[2]^2)
lCI <- pf_24w$surv[1] - pf_24w$surv[2] - 1.96 *diffSE
uCI <- pf_24w$surv[1] - pf_24w$surv[2] + 1.96 *diffSE
zStat <- (pf_24w$surv[1] - pf_24w$surv[2])/diffSE
2*pnorm(abs(zStat), lower.tail=FALSE)

#fit proportional hazards model
coxph_pf <- coxph(Surv(in_wk, out_wk, po) ~ pf_4wk, data = port_po_long_1st)
summary(coxph_pf)

#survival by treatment at beginning of interval (or most recent test)
km_treat <- survfit(Surv(in_wk, out_wk, po) ~ treat_in, data = port_po_long_1st)
plot(km_treat)
summary(km_treat)
sum(km_treat$n.event)
treat_surv <- ggsurvplot(km_treat, data = port_po_long_1st, 
                      risk.table = TRUE, 
                      palette = c("#1F78B4","#E31A1C"), 
                      legend.labs = c("Untreated", "Treated last visit"),
                      legend.title = "",
                      conf.int = TRUE,
                      fun = "event",
                      ylim = c(0,1),
                      censor.shape = "",
                      break.x.by = 2,
                      xlab = "Time (w) since first negative test",
                      ylab = "Cumulative incidence")
treat_surv_main <- treat_surv$plot/treat_surv$table + plot_layout(heights = c(8,2))
treat_surv_main
ggsave("figures/surv_long_treat.png", treat_surv_main, height = 6, width = 6, dpi = 600)

coxph_pf <- coxph(Surv(in_, out_, po) ~ pf_in, data = port_po_long_1st)
summary(coxph_pf)

#individuals treated in last month
km_treat4wk <- survfit(Surv(in_wk, out_wk, po) ~ treat_4wk, data = port_po_long_1st)
plot(km_treat4wk)
summary(km_treat4wk)
sum(km_treat4wk$n.event)
treat4wk_surv <- ggsurvplot(km_treat4wk, data = port_po_long_1st, 
                         risk.table = TRUE, 
                         palette = c("#1F78B4","#E31A1C"), 
                         legend.labs = c("Untreated", "Treated in last mo"),
                         legend.title = "",
                         conf.int = TRUE,
                         fun = "event",
                         ylim = c(0,1),
                         censor.shape = "",
                         break.x.by = 2,
                         xlab = "Time (w) since first negative test",
                         ylab = "Cumulative incidence")
treat4wk_surv_main <- treat4wk_surv$plot/treat4wk_surv$table + plot_layout(heights = c(8,2))
treat4wk_surv_main
ggsave("figures/surv_long_treat4wk.png", treat4wk_surv_main, height = 6, width = 6, dpi = 600)

coxph_treat_4wk <- coxph(Surv(in_, out_, po) ~ treat_4wk, data = port_po_long_1st)
summary(coxph_treat_4wk)

coxph_treat_4wk <- coxph(Surv(in_, out_, po) ~ treat_4wk, data = port_po_long_1st)
summary(coxph_treat_4wk)

#individuals treated in last 2mo
km_treat8wk <- survfit(Surv(in_wk, out_wk, po) ~ treat_8wk, data = port_po_long_1st)
plot(km_treat8wk)
summary(km_treat8wk)
sum(km_treat8wk$n.event)
treat8wk_surv <- ggsurvplot(km_treat8wk, data = port_po_long_1st, 
                            risk.table = TRUE, 
                            palette = c("#1F78B4","#E31A1C"), 
                            legend.labs = c("Untreated", "Treated in last 2mo"),
                            legend.title = "",
                            conf.int = TRUE,
                            fun = "event",
                            ylim = c(0,1),
                            censor.shape = "",
                            break.x.by = 2,
                            xlab = "Time (w) since first negative test",
                            ylab = "Cumulative incidence")
treat8wk_surv_main <- treat8wk_surv$plot/treat8wk_surv$table + plot_layout(heights = c(8,2))
treat8wk_surv_main
ggsave("figures/surv_long_treat8wk.png", treat8wk_surv_main, height = 6, width = 6, dpi = 600)

#individuals treated in last 10wk
km_treat10wk <- survfit(Surv(in_wk, out_wk, po) ~ treat_10wk, data = port_po_long_1st)
plot(km_treat10wk)
summary(km_treat10wk)
sum(km_treat10wk$n.event)
treat10wk_surv <- ggsurvplot(km_treat10wk, data = port_po_long_1st, 
                            risk.table = TRUE, 
                            palette = c("#1F78B4","#E31A1C"), 
                            legend.labs = c("Untreated", "Treated in last 10w"),
                            legend.title = "",
                            conf.int = TRUE,
                            fun = "event",
                            ylim = c(0,1),
                            censor.shape = "",
                            break.x.by = 2,
                            xlab = "Time (w) since first negative test",
                            ylab = "Cumulative incidence")
treat10wk_surv_main <- treat10wk_surv$plot/treat10wk_surv$table + plot_layout(heights = c(8,2))
treat10wk_surv_main
ggsave("figures/surv_long_treat10wk.png", treat10wk_surv_main, height = 6, width = 6, dpi = 600)


####season as time-varying
km_season <- survfit(Surv(in_wk, out_wk, po) ~ long_wet_6w + short_wet_6w, data = port_po_long_1st)
plot(km_season)
summary(km_season)
sum(km_season$n.event)
season_surv <- ggsurvplot(km_season, data = port_po_long_1st, 
                          risk.table = TRUE, 
                          palette = c("#1F78B4", "#6A3D9A","#FF7F00"), 
                          legend.labs = c("Dry", "Short wet", "Long wet"),
                          legend.title = "",
                          conf.int = TRUE,
                          fun = "event",
                          ylim = c(0,1),
                          censor.shape = "",
                          break.x.by = 2,
                          xlab = "Time (w) since first negative test",
                          ylab = "Cumulative incidence")
season_surv_main <- season_surv$plot/season_surv$table + plot_layout(heights=c(6,2))
season_surv_main
ggsave("figures/season_long.png", season_surv_main, height = 4.8, width = 6, dpi = 600)

#test risk difference at final timepoint
#dry vs short
#subset to final week
season_24w <- summary(km_season, times = 24)
#risk difference
rd_dry_v_short <- season_24w$surv[1] - season_24w$surv[2]
#standard error of RD
diffSE_dry_v_short <- sqrt(season_24w$std.err[1]^2 + season_24w$std.err[2]^2)
#z-test
zStat_dry_v_short <- rd_dry_v_short/diffSE_dry_v_short
2*pnorm(abs(zStat_dry_v_short), lower.tail=FALSE)

#dry vs long
#subset to final week with risk set
season_22w <- summary(km_season, times = 22)
#risk difference
rd_dry_v_long <- season_22w$surv[1] - season_22w$surv[3]
#standard error of RD
diffSE_dry_v_long <- sqrt(season_22w$std.err[1]^2 + season_22w$std.err[3]^2)
#z-test
zStat_dry_v_long <- rd_dry_v_long/diffSE_dry_v_long
2*pnorm(abs(zStat_dry_v_long), lower.tail=FALSE)

#long vs short
#subset to final week
season_22w <- summary(km_season, times = 22)
#risk difference
rd_short_v_long <- season_22w$surv[2] - season_22w$surv[3]
#standard error of RD
diffSE_short_v_long <- sqrt(season_22w$std.err[2]^2 + season_22w$std.err[3]^2)
#z-test
zStat_short_v_long <- rd_short_v_long/diffSE_short_v_long
2*pnorm(abs(zStat_short_v_long), lower.tail=FALSE)

coxph_season <- coxph(Surv(in_wk, out_wk, po) ~ long_wet_6w + short_wet_6w, data = port_po_long_1st)
coxph_season
season_data <- data.frame(long_wet_6w = c(0,1,0), short_wet_6w = c(0,0,1))
season_cox_plot <- ggsurvplot(survfit(coxph_season, newdat = season_data), data = port_po_long_1st,
           palette = c("#1F78B4","#FF7F00", "#6A3D9A"), 
           legend.labs = c("Dry", "Long wet", "Short wet"),
           legend.title = "",
           conf.int = TRUE,
           fun = "event",
           ylim = c(0,1),
           censor.shape = "",
           break.x.by = 2,
           xlab = "Time (w) since first negative test",
           ylab = "Cumulative incidence by Proportional Hazard")
season_cox_plot
ggsave("figures/season_ph.png", season_cox_plot, height = 6, width = 6, dpi = 600)


############################
##### MCC ###########
############################

#make new dataset for mcc
port_mcc_po <- port_po_long

#iterate through dataset to update status variable
#note: counting administrative censoring as censoring, not as a permanent competing event (status = 2)
port_mcc_po <- port_mcc_po %>%
  #group by id
  group_by(screening_id) %>%
  #assign status
  mutate(status = 
           case_when(
             #status = 1 when an interval has po positive following negative test or is at start of follow-up
             po == 1 & (lag(po, n = 1, NA) == 0 | in_ == 0) ~ 1,
             #otherwise, is negative
             TRUE ~ 0)) %>%
  #give all individuals weight of 1
  mutate(wt = 1) %>%
  #create flag for po event status
  mutate(flag = "Po")

#make dataset for seasonal evaluation
port_mcc_po_season <- port_po_long_allvisits

#iterate through dataset to update status variable
#note: counting administrative censoring as censoring, not as a permanent competing event (status = 2)
port_mcc_po_season <- port_mcc_po_season %>%
  #group by id
  group_by(screening_id) %>%
  #assign status
  mutate(status = 
           case_when(
             #status = 1 when an interval has po positive following negative test or is at start of follow-up
             po == 1 & (lag(po, n = 1, NA) == 0 | in_ == 0) ~ 1,
             #otherwise, is negative
             TRUE ~ 0)) %>%
  #give all individuals weight of 1
  mutate(wt = 1) %>%
  #create flag for po event status
  mutate(flag = "Po")

#determine season of incident infections
port_mcc_po %>% tabyl(season_6w, status) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

dry_v_long <- chisq.test(matrix(c(333, 31, 229, 14), nrow = 2))
dry_v_long
dry_v_short <- chisq.test(matrix(c(333, 31, 85, 5), nrow = 2))
dry_v_short

dry_v_wet <- chisq.test(matrix(c(333, 31, 314, 19), nrow = 2))
dry_v_wet

port_mcc_po_season %>% filter(season_6w == "dry") %>%
  tabyl(month_vary, status) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#create additional dataset without early recurrences
port_mcc_po_noearly <- port_mcc_po %>%
  filter(po_wk2_early == 0 | (po_wk2_early == 1 & interval != "wk2")) %>%
  filter(po_wk4_early == 0 | (po_wk4_early == 1 & interval != "wk4")) %>%
  #if early recurrence at wk2 or wk4, then subtract 4 weeks from time
  mutate(in_ = ifelse(po_wk2_early == 1 | po_wk4_early == 1, in_ - 28, in_)) %>%
  mutate(out_ = ifelse(po_wk2_early == 1 | po_wk4_early == 1, out_ - 28, out_)) %>%
  mutate(in_wk = ifelse(po_wk2_early == 1 | po_wk4_early == 1, in_wk - 4, in_wk)) %>%
  mutate(out_wk = ifelse(po_wk2_early == 1 | po_wk4_early == 1, out_wk - 4, out_wk)) %>%
  filter(in_ >= 0) %>%
  #can potentially cut the below
  mutate(idx = paste(screening_id,"_noe", sep = "")) %>%
  mutate(arm = 1)

### create a loop for MCC ###
#no competing risks, so KM survival to any interval is 1 (but may be affected by censoring)
mcc_weighted = function(data_nested) {
  n0 = nrow(data_nested)
  data_nested %>%
    unnest(cols = c(data)) %>%
    group_by(out_wk) %>%
    summarise(e_j = sum(status),
              #r_jplus1 = sum(lead.death*baseline_wt),
              c_jplus1_mcc = sum(lead.censor_mcc)) %>%
              #c_jplus1_km  = sum(lead.censor_km*baseline_wt))
    mutate(
            #r_j = lag(r_jplus1, default = 0),
           c_j_mcc = lag(c_jplus1_mcc, default = 0), 
           #c_j_km  = lag(c_jplus1_km, default = 0),
           n_jminus1_mcc = n0 - cumsum(c_j_mcc),
           #n_jminus1_km  = n0 - cumsum(r_j + c_j_km),
           #km = cumprod(1 - r_j/n_jminus1_km),
           mcc = cumsum(e_j/n_jminus1_mcc)) %>%
    dplyr::select(out_wk, mcc)
}

#created nested data by individual
port_mcc_po_nested = port_mcc_po %>% 
  # group by individual
  group_by(screening_id) %>%
  #determine times
  mutate(lead.censor_mcc = 1*(out_ == max(out_))) %>%
      dplyr::select(screening_id, wt, out_wk, status, lead.censor_mcc, pf_cont) %>%
      group_by(screening_id) %>%
  nest()

mcc_estimates_overall = mcc_weighted(port_mcc_po_nested) %>% rename(estimate = mcc)


### Set number of bootstrap samples for confidence intervals ###
set.seed(59)
n_boots = 1000
n_ids = nrow(port_mcc_po_nested)

mcc_bootstraps = tibble(rep         = 1:n_boots,
                        data_nested = replicate(n        = n_boots,
                                                expr     = port_mcc_po_nested[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(mcc_overall = map(data_nested, mcc_weighted))

### overall data ###
mcc_overall = mcc_bootstraps %>%
  dplyr::select(-data_nested) %>%
  unnest(mcc_overall) %>%
  group_by(out_wk) %>%
  summarise(se = sd(mcc))

mc_ci_overall = mcc_estimates_overall %>%
  left_join(mcc_overall) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate + 1.96*se) %>%
  #create flag for early recurrence exclusions
  mutate(early = "All")


mc_ci_overall

#no early recurrences: created nested data by individual
port_mcc_po_noearly_nested = port_mcc_po_noearly %>% 
  # group by individual
  group_by(screening_id) %>%
  #determine times
  mutate(lead.censor_mcc = 1*(out_ == max(out_))) %>%
  dplyr::select(screening_id, wt, out_wk, status, lead.censor_mcc, pf_cont) %>%
  group_by(screening_id) %>%
  nest()

mcc_noearly_estimates_overall = mcc_weighted(port_mcc_po_noearly_nested) %>% rename(estimate = mcc)

### Set number of bootstrap samples for confidence intervals ###
set.seed(59)
n_boots = 1000
n_ids = nrow(port_mcc_po_noearly_nested)

mcc_bootstraps = tibble(rep         = 1:n_boots,
                        data_nested = replicate(n        = n_boots,
                                                expr     = port_mcc_po_noearly_nested[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(mcc_noearly_overall = map(data_nested, mcc_weighted))

### overall data ###
mcc_noearly_overall = mcc_bootstraps %>%
  dplyr::select(-data_nested) %>%
  unnest(mcc_noearly_overall) %>%
  group_by(out_wk) %>%
  summarise(se = sd(mcc))

mc_noearly_ci_overall = mcc_noearly_estimates_overall %>%
  left_join(mcc_noearly_overall) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate + 1.96*se) %>%
  #add flag for comparison
  mutate(early = "Excluded")


mc_comb_ci_overall <- rbind(mc_ci_overall, mc_noearly_ci_overall) %>%
  #add empty rows for starting point
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, early = "All") %>%
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, early = "Excluded")

#plot MCC for cohort overall
mcc_plot <- ggplot(data = mc_comb_ci_overall) +
  geom_step(aes(x = out_wk, y = estimate, color = as.factor(early))) +
  geom_stepconfint(aes(x = out_wk, ymin = lwr, ymax = upr, fill = as.factor(early)), alpha = 0.3) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,24, by = 2)) +
  scale_y_continuous(breaks = seq(0,1, by = 0.1)) +
  xlab("Time (w) since first negative test") + 
  ylab("Mean cumulative count") + 
  scale_color_manual(values = c("All" = "#1F78B4", "Excluded" = "#B15928")) + 
  scale_fill_manual(values = c("All" = "#1F78B4", "Excluded" = "#B15928")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14),
        legend.position = "inside", legend.position.inside = c(0.5, 0.9),
        panel.grid = element_blank())

mcc_plot

ggsave("figures/mcc_po.png", mcc_plot, height = 4, width = 5, dpi = 600)

###############################
#### Pf cumulative count ######
###############################

#duplicate observations, one per po observation, by po bi-weekly screening variables
port_pf_long_allvisits <- port %>%
  pivot_longer(c(po_enroll, po_wk2, po_wk4, po_wk6, po_wk8, po_wk10, po_wk12, po_wk14, po_wk16, po_wk18, po_wk20, po_wk22, po_wk24,
                 pf_enroll, pf_wk2, pf_wk4, pf_wk6, pf_wk8, pf_wk10, pf_wk12, pf_wk14, pf_wk16, pf_wk18, pf_wk20, pf_wk22, pf_wk24),
               names_to = c(".value", "interval"), names_sep = "_") %>%
  #create time and time2 variables corresponding to each interval
  mutate(in_ = case_when(interval == "enroll" ~ -14, interval == "wk2" ~ 0, interval == "wk4" ~ 14, 
                         interval == "wk6" ~ 28, interval == "wk8" ~ 42, 
                         interval == "wk10" ~ 56, interval == "wk12" ~ 70,
                         interval == "wk14" ~ 84, interval == "wk16" ~ 98,
                         interval == "wk18" ~ 112, interval == "wk20" ~ 126,
                         interval == "wk22" ~ 140, interval == "wk24" ~ 154)) %>%
  mutate(out_ = in_ + 14) %>%
  #add in treatment data, reflecting treatment at start of interval
  mutate(treat_in = case_when(interval == "enroll" ~ NA, interval == "wk2" ~ port_sample_treatyn_do, interval == "wk4" ~ port_sample_treatyn_w2, 
                              interval == "wk6" ~ port_sample_treatyn_w4, interval == "wk8" ~ port_sample_treatyn_w6, 
                              interval == "wk10" ~ port_sample_treatyn_w8, interval == "wk12" ~ port_sample_treatyn_w10,
                              interval == "wk14" ~ port_sample_treatyn_w12, interval == "wk16" ~ port_sample_treatyn_w14,
                              interval == "wk18" ~ port_sample_treatyn_w16, interval == "wk20" ~ port_sample_treatyn_w18,
                              interval == "wk22" ~ port_sample_treatyn_w20, interval == "wk24" ~ port_sample_treatyn_w22)) %>%
  #for missed visits, carry forward missing data
  group_by(screening_id) %>%
  fill(po, .direction = c("down")) %>%
  fill(pf, .direction = c("down")) %>%
  #if enrollment data is still missing, draw from screening
  mutate(po = ifelse(is.na(po) & interval == "enroll", po_screening, po)) %>%
  mutate(pf = ifelse(is.na(pf) & interval == "enroll", pf_screening, pf)) %>%
  #encode lagged po-positivity, reflecting positive at start of interval (includes forward-carrying)
  mutate(po_in = ifelse(interval == "enroll", po_screening, lag(po, n = 1, NA))) %>%
  #individuals with missing treatment data were not treated
  mutate(treat_in = ifelse(is.na(treat_in), 0, treat_in)) %>%
  #encode LTFU as a drop variable representing the final visit in which they were screened
  mutate(drop = ifelse(!is.na(t_visit) & out_ == t_visit, 1, 0)) %>%
  #remove intervals after t_visit (LTFU) indexed from enrollment (not first negative test)
  filter(out_ <= t_visit | is.na(t_visit)) %>%
  #calculate month 6w prior to out_
  mutate(date_6w_lag = as.Date(as.Date(port_sample_date_d0) + out_ - 42)) %>%
  mutate(month_vary = as.numeric(format(date_6w_lag, "%m"))) %>%
  #determine season 6w prior to each test
  mutate(season_6w = case_when(
    month_vary >= 10 & month_vary <= 12 ~ "shortwet",
    month_vary >= 3 & month_vary <= 5 ~ "longwet",
    TRUE ~ "dry")) %>%
  mutate(long_wet_6w = case_when(
    month_vary >= 3 & month_vary <= 5 ~ 1, 
    TRUE ~ 0)) %>%
  mutate(short_wet_6w = case_when(
    month_vary >= 10 & month_vary <= 12 ~ 1, 
    TRUE ~ 0)) %>%
  #calculate the month 6w prior to interval onset
  mutate(month_onset = as.numeric(format(as.Date(as.Date(port_sample_date_d0) + in_ - 42), "%m"))) %>%
  #determine season 6w prior to interval onset
  mutate(season_6w_onset = case_when(
    month_onset >= 10 & month_onset <= 12 ~ "shortwet",
    month_onset >= 3 & month_onset <= 5 ~ "longwet",
    TRUE ~ "dry")) %>%
  mutate(long_wet_6w_onset = case_when(
    month_onset >= 3 & month_onset <= 5 ~ 1, 
    TRUE ~ 0)) %>%
  mutate(short_wet_6w_onset = case_when(
    month_onset >= 10 & month_onset <= 12 ~ 1, 
    TRUE ~ 0)) %>%
  #calculate month no lag prior to out_
  mutate(month_nolag = as.numeric(format(as.Date(as.Date(port_sample_date_d0) + out_), "%m"))) %>%
  #determine season no lag prior to each test
  mutate(season_nolag = case_when(
    month_nolag >= 10 & month_nolag <= 12 ~ "short wet",
    month_nolag >= 3 & month_nolag <= 5 ~ "long wet",
    TRUE ~ "dry")) %>%
  mutate(long_wet_nolag = case_when(
    month_nolag >= 3 & month_nolag <= 5 ~ 1, 
    TRUE ~ 0)) %>%
  mutate(short_wet_nolag = case_when(
    month_nolag >= 10 & month_nolag <= 12 ~ 1, 
    TRUE ~ 0)) %>%
  #subtract 1st poneg from each interval so they reflect time since becoming at risk
  mutate(in_ = in_ - t_1stpfneg) %>%
  mutate(out_ = out_ - t_1stpfneg) %>%
  #recode days as weeks
  mutate(in_wk = in_/7) %>%
  mutate(out_wk = out_/7)


port_pf_long <- port_pf_long_allvisits %>% 
  #remove rows with negative time
  filter(in_ >= 0) %>%
  filter(out_ > 0)

#iterate through dataset to update status variable
#note: counting administrative censoring as censoring, not as a permanent competing event (status = 2)
port_mcc_pf <- port_pf_long %>%
  #group by id
  group_by(screening_id) %>%
  #assign status
  mutate(status = 
           case_when(
             #status = 1 when an interval has po positive following negative test or is at start of follow-up
             pf == 1 & (lag(pf, n = 1, NA) == 0 | in_ == 0) ~ 1,
             #otherwise, is negative
             TRUE ~ 0)) %>%
  #give all individuals weight of 1
  mutate(wt = 1) %>%
  #create flag for po event status
  mutate(flag = "Pf") %>%
  #change screening_id to plot separately from po
  mutate(screening_id = paste(screening_id,"_pf", sep = ""))

#make mcc dataset using all visits
port_mcc_pf_season <- port_pf_long_allvisits %>%
  #group by id
  group_by(screening_id) %>%
  #assign status
  mutate(status = 
           case_when(
             #status = 1 when an interval has po positive following negative test or is at start of follow-up
             pf == 1 & (lag(pf, n = 1, NA) == 0 | in_ == 0) ~ 1,
             #otherwise, is negative
             TRUE ~ 0)) %>%
  #give all individuals weight of 1
  mutate(wt = 1) %>%
  #create flag for po event status
  mutate(flag = "Pf") %>%
  #change screening_id to plot separately from po
  mutate(screening_id = paste(screening_id,"_pf", sep = ""))

#determine season of incident infections
port_mcc_pf %>% tabyl(season_6w, status) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#test season association
dry_v_long <- chisq.test(matrix(c(325, 29, 208, 28), nrow = 2))
dry_v_long
dry_v_short <- chisq.test(matrix(c(325, 29, 68, 16), nrow = 2))
dry_v_short

dry_v_wet <- chisq.test(matrix(c(325, 29, 276, 44), nrow = 2))
dry_v_wet

port_mcc_pf %>% filter(season_6w == "dry") %>%
  tabyl(month_vary, status) %>%
  adorn_totals(c("row")) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(2) %>%
  adorn_ns(position = "front")

#merge po and pf datasets
port_mcc_popf <- rbind(port_mcc_po, port_mcc_pf)

#created nested data by individual
port_mcc_popf_nested = port_mcc_popf %>% 
  # group by individual
  group_by(screening_id) %>%
  #determine times
  mutate(lead.censor_mcc = 1*(out_ == max(out_))) %>%
  dplyr::select(screening_id, wt, out_wk, status, lead.censor_mcc, pf_cont, flag) %>%
  group_by(screening_id, flag) %>%
  nest()

### calculate point estimates ###
mcc_byspecies = function(port_mcc_popf_nested) {
  port_mcc_popf_nested %>% 
    group_by(flag) %>%
    nest() %>%
    mutate(mcc = map(data, mcc_weighted)) %>%
    dplyr::select(-data) %>%
    unnest(mcc) %>%
    pivot_wider(names_from = flag, values_from = mcc) %>%
    mutate(popf_diff = Pf - Po)
}


mcc_estimates_byspecies = mcc_byspecies(port_mcc_popf_nested)


### bootstrap ###
set.seed(59)
n_boots = 1000
n_ids = nrow(port_mcc_popf_nested)

mcc_bootstraps = tibble(rep         = 1:n_boots,
                        data_nested = replicate(n        = n_boots,
                                                expr     = port_mcc_popf_nested[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(mcc_overall = map(data_nested, mcc_weighted),
         mcc_byspecies = map(data_nested, mcc_byspecies))

mcc_errors_byspecies = mcc_bootstraps %>%
  select(-data_nested) %>%
  unnest(mcc_byspecies) %>%
  group_by(out_wk) %>%
  summarise(Po = sd(Po),
            Pf = sd(Pf),
            popf_diff = sd(popf_diff))


mcc_ci_byspecies = mcc_estimates_byspecies %>%
  pivot_longer(-out_wk, names_to = 'parameter', values_to='estimate') %>%
  left_join(mcc_errors_byspecies %>% pivot_longer(-out_wk, names_to='parameter', values_to='se')) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate + 1.96*se)

mcc_ci_byspecies[mcc_ci_byspecies$out_wk == 24 & mcc_ci_byspecies$parameter == "popf_diff",]
zStat <- as.numeric((mcc_ci_byspecies[mcc_ci_byspecies$out_wk == 24 & mcc_ci_byspecies$parameter == "popf_diff", c("estimate")])/mcc_ci_byspecies[mcc_ci_byspecies$out_wk == 24 & mcc_ci_byspecies$parameter == "popf_diff", c("se")])
2*pnorm(abs(zStat), lower.tail=FALSE)

mcc_byspecies_plotting <- mcc_ci_byspecies %>%
  filter(parameter != "popf_diff") %>%
  #add empty rows for starting point
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, parameter = "Po") %>%
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, parameter = "Pf") %>%
  mutate(lwr = case_when(lwr < 0 ~ 0, TRUE ~ lwr))

#plot MCC of Po and Pf
mcc_plot <- ggplot(data = mcc_byspecies_plotting) +
  geom_step(aes(x = out_wk, y = estimate, color = as.factor(parameter))) +
  geom_stepconfint(aes(x = out_wk, ymin = lwr, ymax = upr, fill = as.factor(parameter)), alpha = 0.3) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,24, by = 2)) +
  scale_y_continuous(breaks = seq(0,1.5, by = 0.25)) +
  xlab("Time (w) since first negative test") + 
  ylab("Mean cumulative count") + 
  scale_color_manual(values = c("Po" = "#1F78B4", "Pf" = "#E31A1C"), labels = c("Po" = expression(paste(italic("P. ovale"), " spp.")), "Pf" = expression(paste(italic("P. falciparum"))))) + 
  scale_fill_manual(values = c("Po" = "#1F78B4", "Pf" = "#E31A1C"),  labels = c("Po" = expression(paste(italic("P. ovale"), " spp.")), "Pf" = expression(paste(italic("P. falciparum"))))) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.9),
        panel.grid = element_blank())

mcc_plot

ggsave("figures/mcc_popf.png", mcc_plot, height = 4, width = 5, dpi = 600)

##############################################
###### Incidence Rate #######################
##############################################

#adjust nested dataset to show total events and follow-up time per observation
port_inc_popf_nested <- port_mcc_popf_nested %>%
  unnest(cols = c(data)) %>%
  group_by(screening_id) %>%
  mutate(total_event = sum(status)) %>%
  filter(out_wk == max(out_wk)) %>%
  dplyr::select(screening_id, flag, out_wk, total_event) 

#define function for calculating incidence
incidence_rate = function(data) {
  data %>%
  group_by(flag) %>%
  summarise(rate = sum(total_event)/sum(out_wk)) %>%
  pivot_wider(names_from = flag, values_from = rate) %>%
  mutate(rate_diff = Pf - Po) 
}  


rates <- incidence_rate(port_inc_popf_nested)
rates

### bootstrap ###
set.seed(59)
n_boots = 1000
n_ids = nrow(port_inc_popf_nested)

inc_bootstraps = tibble(rep         = 1:n_boots,
                        data_nested = replicate(n        = n_boots,
                                                expr     = port_inc_popf_nested[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(inc_byspecies = map(data_nested, incidence_rate))

inc_errors_byspecies = inc_bootstraps %>%
  dplyr::select(-data_nested) %>%
  unnest(inc_byspecies) %>%
  reframe(Po = sd(Po),
            Pf = sd(Pf),
            rate_diff = sd(rate_diff)) 

inc_ci_byspecies <- rates %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to='estimate') %>%
  left_join(inc_errors_byspecies %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to='se')) %>%
  mutate(lwr = 52*(estimate - 1.96*se), upr = 52*(estimate + 1.96*se)) %>%
  mutate(rate_year = estimate*52)
inc_ci_byspecies

zStat <- as.numeric((inc_ci_byspecies[inc_ci_byspecies$parameter == "rate_diff", c("estimate")])/inc_ci_byspecies[inc_ci_byspecies$parameter == "rate_diff", c("se")])
2*pnorm(abs(zStat), lower.tail=FALSE)

###############################################
########## Incidence by season ################
###############################################

#define function for calculating incidence by season
incidence_rate_season = function(data) {
  data %>%
    group_by(flag, season_6w) %>%
    summarise(rate = sum(total_event)/sum(out_wk)) %>%
    pivot_wider(names_from = flag, values_from = rate)
}  

#make dataset for incidence by season
popf_season_inc <- port_mcc_popf %>% 
  group_by(season_6w, flag) %>%
  summarise(total_event = sum(status), weeks = n()*2) %>%
  mutate(rate = (total_event/weeks)*52)
popf_season_inc

#calculate ratio of ratio of Pf/Po incident infection identification by season
season_ratio = function(data) {
  data %>%
  mutate(wetdry_6w_name = case_when(season_6w == "longwet" ~ "wet", season_6w == "shortwet" ~ "wet", season_6w == "dry" ~ "dry")) %>%
  group_by(wetdry_6w_name, flag) %>%
  summarise(positivity = sum(status)/n(), .groups = 'keep') %>%
  pivot_wider(names_from = wetdry_6w_name, values_from = positivity) %>%
  mutate(ratio = dry/wet) %>%
  dplyr::select(-dry, -wet) %>%
  pivot_wider(names_from = flag, values_from = ratio) %>%
  mutate(ratioratio = Po/Pf)
}  

ratios <- season_ratio(port_mcc_popf)
ratios

#perform bootstrap to estimate variance
### bootstrap ###
set.seed(59)
n_boots = 1000

n_ids = nrow(port_mcc_popf)

ratio_bootstraps = tibble(rep         = 1:n_boots,
                        data = replicate(n        = n_boots,
                                                expr     = port_mcc_popf[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(season_ratios = map(data, season_ratio)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(season_ratios))

#compile standard deviation of bootstrap simulations
ratio_errors = ratio_bootstraps %>%
  mutate(logPo = log(Po),
         logPf = log(Pf),
         logratioratio = log(ratioratio)) %>%
  reframe(logPo = sd(logPo),
          logPf = sd(logPf),
          logratioratio = sd(logratioratio)) 

ratio_ci <- ratios %>%
  mutate(logPo = log(Po),
         logPf = log(Pf),
         logratioratio = log(ratioratio)) %>%
  dplyr::select(-Po, -Pf, -ratioratio) %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to='logestimate') %>%
  left_join(ratio_errors %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to='logse')) %>%
  mutate(loglwr = (logestimate - 1.96*logse), logupr = (logestimate + 1.96*logse)) %>%
  mutate(estimate = exp(logestimate),
         lwr = exp(loglwr),
         upr = exp(logupr))
ratio_ci

zStat <- as.numeric((ratio_ci[ratio_ci$parameter == "logratioratio", c("logestimate")])/ratio_ci[ratio_ci$parameter == "logratioratio", c("logse")])
2*pnorm(abs(zStat), lower.tail=FALSE)

###############################################
########## MCC by Pf infection ################
###############################################

#created nested data by individual
port_mcc_pfinf_nested = port_mcc_po %>% 
  # group by individual
  group_by(screening_id) %>%
  #determine times
  mutate(lead.censor_mcc = 1*(out_ == max(out_)),
         pf_study_name = case_when(pf_study == 1 ~ "Pfpos", pf_study == 0 ~ "Pfneg")) %>%
  dplyr::select(screening_id, wt, out_wk, status, lead.censor_mcc, pf_study_name) %>%
  group_by(screening_id, pf_study_name) %>%
  nest()

### calculate point estimates ###
mcc_bypfinf = function(port_mcc_pfinf_nested) {
  port_mcc_pfinf_nested %>% 
    group_by(pf_study_name) %>%
    nest() %>%
    mutate(mcc = map(data, mcc_weighted)) %>%
    dplyr::select(-data) %>%
    unnest(mcc) %>%
    pivot_wider(names_from = pf_study_name, values_from = mcc) %>%
    mutate(pf_diff = Pfpos - Pfneg)
}

mcc_estimates_bypfinf = mcc_bypfinf(port_mcc_pfinf_nested)

### bootstrap ###
set.seed(34)
n_boots = 1000

n_ids = nrow(port_mcc_treatpf_nested)

#perform bootstrap
mcc_bootstraps = tibble(rep         = 1:n_boots,
                        data_nested = replicate(n        = n_boots,
                                                expr     = port_mcc_treatpf_nested[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(mcc_bypfinf = map(data_nested, mcc_bypfinf))

#extract errors
mcc_errors_bypfinf = mcc_bootstraps %>%
  select(-data_nested) %>%
  unnest(mcc_bypfinf) %>%
  group_by(out_wk) %>%
  summarise(Pfpos = sd(Pfpos),
            Pfneg = sd(Pfneg),
            pf_diff = sd(pf_diff))

#create confidence intervals
mcc_ci_bypfinf = mcc_estimates_bypfinf %>%
  pivot_longer(-out_wk, names_to = 'parameter', values_to='estimate') %>%
  left_join(mcc_errors_bypfinf %>% pivot_longer(-out_wk, names_to='parameter', values_to='se')) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate + 1.96*se)

#z-test
mcc_ci_bypfinf[mcc_ci_bypfinf$out_wk == 24 & mcc_ci_bypfinf$parameter == "pf_diff",]
zStat <- as.numeric((mcc_ci_bypfinf[mcc_ci_bypfinf$out_wk == 24 & mcc_ci_bypfinf$parameter == "pf_diff", c("estimate")])/mcc_ci_bypfinf[mcc_ci_bypfinf$out_wk == 24 & mcc_ci_bypfinf$parameter == "pf_diff", c("se")])
2*pnorm(abs(zStat), lower.tail=FALSE)

#prepare for plotting
mcc_bypfinf_plotting <- mcc_ci_bypfinf %>%
  filter(parameter != "pf_diff") %>%
  #add empty rows for starting point
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, parameter = "Pfpos") %>%
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, parameter = "Pfneg")

#plot MCC of Pf+ and Pf-
mcc_plot <- ggplot(data = mcc_bypfinf_plotting) +
  geom_step(aes(x = out_wk, y = estimate, color = as.factor(parameter))) +
  geom_stepconfint(aes(x = out_wk, ymin = lwr, ymax = upr, fill = as.factor(parameter)), alpha = 0.3) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,24, by = 2)) +
  scale_y_continuous(breaks = seq(0,1.5, by = 0.25)) +
  xlab("Time (w) since first negative test") + 
  ylab("Mean cumulative count") + 
  scale_color_manual(values = c("Pfneg" = "#1F78B4", "Pfpos" = "#E31A1C"), labels = c("Pfneg" = "Pf-", "Pfpos" = "Pf+")) + 
  scale_fill_manual(values = c("Pfneg" = "#1F78B4", "Pfpos" = "#E31A1C"), labels = c("Pfneg" = "Pf-", "Pfpos" = "Pf+")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.9),
        panel.grid = element_blank())

mcc_plot

ggsave("figures/mcc_po_pfinf.png", mcc_plot, height = 4, width = 5, dpi = 600)

###############################################
########## MCC by Pf treatment ################
###############################################

#created nested data by individual
port_mcc_treatpf_nested = port_mcc_po %>% 
  # group by individual
  group_by(screening_id) %>%
  #determine times
  mutate(lead.censor_mcc = 1*(out_ == max(out_)),
         treat_pf_name = case_when(treat_pf_study == 1 ~ "Pftreat", treat_pf_study == 0 ~ "NoPftreat")) %>%
  dplyr::select(screening_id, wt, out_wk, status, lead.censor_mcc, treat_pf_name) %>%
  group_by(screening_id, treat_pf_name) %>%
  nest()

#investigate incident po infections among individuals treated for pf
km_pftreat_3mo <- survfit(Surv(in_wk, out_wk, status) ~ treat_pf_3mo, data = port_mcc_po)
km_pftreat_3mo
print(port_mcc_po[port_mcc_po$status == 1 & port_mcc_po$treat_pf_3mo == 1, c("screening_id", "interval")])

mcc_bypftreat = function(port_mcc_treatpf_nested) {
  port_mcc_treatpf_nested %>% 
    group_by(treat_pf_name) %>%
    nest() %>%
    mutate(mcc = map(data, mcc_weighted)) %>%
    dplyr::select(-data) %>%
    unnest(mcc) %>%
    pivot_wider(names_from = treat_pf_name, values_from = mcc) %>%
    mutate(treat_diff = Pftreat - NoPftreat)
}

mcc_estimates_bypftreat = mcc_bypftreat(port_mcc_treatpf_nested)

### bootstrap ###
set.seed(34)
n_boots = 1000

n_ids = nrow(port_mcc_treatpf_nested)

#perform bootstrap
mcc_bootstraps = tibble(rep         = 1:n_boots,
                        data_nested = replicate(n        = n_boots,
                                                expr     = port_mcc_treatpf_nested[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(mcc_bypftreat = map(data_nested, mcc_bypftreat))


mcc_errors_bypftreat = mcc_bootstraps %>%
  select(-data_nested) %>%
  unnest(mcc_bypftreat) %>%
  group_by(out_wk) %>%
  summarise(Pftreat = sd(Pftreat),
            NoPftreat = sd(NoPftreat),
            treat_diff = sd(treat_diff))

mcc_ci_bypftreat = mcc_estimates_bypftreat %>%
  pivot_longer(-out_wk, names_to = 'parameter', values_to='estimate') %>%
  left_join(mcc_errors_bypftreat %>% pivot_longer(-out_wk, names_to='parameter', values_to='se')) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate + 1.96*se)

mcc_ci_bypftreat[mcc_ci_bypftreat$out_wk == 24,]
zStat <- as.numeric((mcc_ci_bypftreat[mcc_ci_bypftreat$out_wk == 24 & mcc_ci_bypftreat$parameter == "treat_diff", c("estimate")])/mcc_ci_bypftreat[mcc_ci_bypftreat$out_wk == 24 & mcc_ci_bypftreat$parameter == "treat_diff", c("se")])
2*pnorm(abs(zStat), lower.tail=FALSE)

mcc_bypftreat_plotting <- mcc_ci_bypftreat %>%
  filter(parameter != "treat_diff") %>%
  #add empty rows for starting point
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, parameter = "Pftreat") %>%
  add_row(out_wk = 0, estimate = 0, se = 0, lwr = 0, upr = 0, parameter = "NoPftreat") %>%
  mutate(lwr = case_when(lwr < 0 ~ 0, TRUE ~ lwr))

#plot MCC of Po and Pf
mcc_plot <- ggplot(data = mcc_bypftreat_plotting) +
  geom_step(aes(x = out_wk, y = estimate, color = as.factor(parameter))) +
  geom_stepconfint(aes(x = out_wk, ymin = lwr, ymax = upr, fill = as.factor(parameter)), alpha = 0.3) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,24, by = 2)) +
  scale_y_continuous(breaks = seq(0,1.5, by = 0.25)) +
  xlab("Time (w) since first negative test") + 
  ylab("Mean cumulative count") + 
  scale_color_manual(values = c("NoPftreat" = "#1F78B4", "Pftreat" = "#E31A1C"), labels = c("NoPftreat" = "No Pf Treatment (n = 59)", "Pftreat" = "Received Pf Treatment (n = 15)")) + 
  scale_fill_manual(values = c("NoPftreat" = "#1F78B4", "Pftreat" = "#E31A1C"), labels = c("NoPftreat" = "No Pf Treatment (n = 59)", "Pftreat" = "Received Pf Treatment (n = 15)")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.9),
        panel.grid = element_blank())

mcc_plot

ggsave("figures/mcc_po_treatpf.png", mcc_plot, height = 4, width = 5, dpi = 600)
