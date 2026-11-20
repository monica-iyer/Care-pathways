

#--Monica Iyer
#--August 2026
#--Analysis for care pathways following infections study

#This script contains the code for the multi-state modelling
#Part 1: preparing the data for the functions in the mstate() package
#Part 2: Cox regression (Figure 2)
#Part 3: predicting state occupancy (Figure 3 and Table 2)
#Part 4: Secondary analyses stratified by sex, age, living status (Figure S12, Table S4, Figure S13)
#Part 5: code to create a forest plot of sequentially adjusted models (Figure S3)


#Draws from "Tutorial in biostatistics: competing risks and multi-state models" (Putter, 2007)

library(dplyr)
library(arrow)
library(table1)
library(tidyr)
library(mstate)
library(survival)
library(readr)
library(tidyverse);library(mstate);library(janitor);library(magrittr);
library(gtsummary);library(labelled);library(mstate);
library(purrr);library(glue);library(gt);library(readr);library(readxl);
library(geepack);library(evaluate);library(meta);library(forestplot);
library(cowplot);library(rio);library(mstate);library(Hmisc);
library(rsconnect);library(forestploter);library(ggpubr);library(msm);
library(RColorBrewer);library(viridisLite);library(grid)
library(ggsci)


#----Part 1: prepare data----

##----read in data----

#read in df and outcomes that were cleaned in descriptive_analysis script

#df has patid, indexdate, and all baseline covariates
#outcomes has event-level information
#patients can enter the cohort multiple times; each patid/indexdate combination is unique

df <- read_parquet("")
outcomes <- read_parquet("")

##----format data for msprep()----

###----add transition indicators----
# indicate when people entered each state and at what date
# for simplicity, I am doing this step first among readmission events
#NOTE only including one readmission in the multistate model

readmission_dates <-
  outcomes %>%
  #filter to readmission events
  filter(readmission == 1) %>%
  #add aurum information about observation and prior care home
  left_join(df %>% select(patid, indexdate, observation_end, care_home),
            by=c("patid", "indexdate")) %>%
  group_by(patid, indexdate) %>%
  arrange(spell_start, .by_group = TRUE) %>%   
  # for each patid/index date, keep only the FIRST readmission row
  slice(1) %>%  
  #remove hospitalisations with no spell end
  filter(!is.na(spell_end)) %>%
  mutate(
    #Create indicator for discharged home after readmission
    #definition: discharge destination is home OR it's care home if the admission source or aurum record is care home
    read_discharge_home = if_else(disdest == "Home" | 
                                    (disdest=="Care home" & (care_home==1 | 
                                                               admisorc=="Care home")), 1L, 0L), 
    #add date of discharge home
    date_read_discharge_home = if_else(disdest == "Home" | 
                                         (disdest=="Care home" & (care_home==1 | 
                                                                    admisorc=="Care home")), spell_end, NA),
    #Create indicator for discharged to care home after readmission
    #definition: discharge destination is care home only if admit source and aurum record is no care home
    read_discharge_care_home = if_else(disdest=="Care home" & 
                                         (is.na(care_home) | care_home==0) &
                                         admisorc!="Care home", 1L, 0L), 
    #add date of discharge to care home
    date_read_discharge_care_home = if_else(disdest=="Care home" & 
                                              (is.na(care_home) | care_home==0) &
                                              admisorc!="Care home", spell_end, NA)
  ) %>%
  #indicate which variables to report since after grouping by patid, indexdate
  #each patid/indexdate should have just one of each variable, so reporting the 'first'
  summarise(
    readmission = 1L,
    readmission_date = first(spell_start),
    read_discharge_home = first(read_discharge_home),
    date_read_discharge_home = first(date_read_discharge_home),
    read_discharge_care_home = first(read_discharge_care_home),
    date_read_discharge_care_home = first(date_read_discharge_care_home),
    .groups = "drop"
  )

#create indicators for death (0,1) and set date of death
#using ONS as gold standard, if missing then CPRD deathdate, and if that is missing then HES discharge destination==died
outcome_death <- outcomes %>%
  left_join(df %>% select(patid, indexdate, dementia, cprd_ddate), by=c("patid", "indexdate")) %>% 
  mutate(death=if_else(!is.na(reg_date_of_death) & (reg_date_of_death-indexdate)<=365.25 | 
                         (is.na(reg_date_of_death) & !is.na(cprd_ddate) & (cprd_ddate-indexdate)<=365.25) |
                         is.na(reg_date_of_death) & is.na(cprd_ddate) & disdest=="Died", 1, 0),
         date_death=case_when(death==1 & !is.na(reg_date_of_death) ~ reg_date_of_death,
                              death==1 & is.na(reg_date_of_death) & !is.na(cprd_ddate) ~ cprd_ddate,
                              death==1 & is.na(reg_date_of_death) & is.na(cprd_ddate) ~ spell_end,
                              death==0 ~ NA)) %>%
  group_by(patid, indexdate) %>%
  summarise(death = max(death), date_death = min(date_death, na.rm = TRUE), .groups = "drop") %>%
  mutate(date_death = if_else(is.infinite(date_death), as.Date(NA), date_death))


#create indicators for discharge home/care home
#then join with the df of readmission and death indicators
outcomes_mstate <- 
  outcomes %>% 
  #filter to non readmission events (index hosp)
  filter(readmission==0) %>% 
  left_join(df %>% select(patid, indexdate, observation_end, care_home),
            by=c("patid", "indexdate")) %>%
  mutate(indexhosp=1L,
         discharge_home=if_else(disdest == "Home" | 
                                  (disdest=="Care home" & (care_home==1 | 
                                                             admisorc=="Care home")), 1, 0),
         date_discharge_home=if_else(disdest == "Home" | 
                                       (disdest=="Care home" & (care_home==1 | 
                                                                  admisorc=="Care home")), spell_end, NA),
         discharge_care_home=if_else(disdest=="Care home" & 
                                       (is.na(care_home) | care_home==0) &
                                       admisorc!="Care home", 1, 0),
         date_discharge_care_home=if_else(disdest=="Care home" & 
                                            (is.na(care_home) | care_home==0) &
                                            admisorc!="Care home", spell_end, NA)) %>%
  select(-readmission) %>%
  
  #add the information from the readmission events
  left_join(readmission_dates, by=c("patid", "indexdate")) %>%
  #if they did not have a readmission, replace NA with 0
  mutate(across(c(readmission, read_discharge_home, read_discharge_care_home),
                ~ replace_na(.x, 0L))) %>%
  left_join(outcome_death %>% select(patid, indexdate, death, date_death), by=c("patid", "indexdate")) %>%
  
  #create a variable of the censoring date: either observation end or death date
  mutate(censor_date = pmin(observation_end, date_death, indexdate + 365.25, na.rm = TRUE)) %>%  #if any event happens AFTER censoring date, make the event indicator 0 and the date NA
  mutate(
    discharge_home = if_else(discharge_home == 1 & date_discharge_home > censor_date, 0L, discharge_home),
    date_discharge_home = if_else(discharge_home == 1, date_discharge_home, NA),
    
    discharge_care_home = if_else(discharge_care_home == 1 & date_discharge_care_home > censor_date, 0L, discharge_care_home),
    date_discharge_care_home = if_else(discharge_care_home == 1, date_discharge_care_home, NA),
    
    readmission = if_else(readmission == 1 & readmission_date > censor_date, 0L, readmission),
    readmission_date = if_else(readmission == 1, readmission_date, NA),
    
    read_discharge_home = if_else(read_discharge_home == 1 & date_read_discharge_home > censor_date, 0L, read_discharge_home),
    date_read_discharge_home = if_else(read_discharge_home == 1, date_read_discharge_home, NA),
    
    read_discharge_care_home = if_else(read_discharge_care_home == 1 & date_read_discharge_care_home > censor_date, 0L, read_discharge_care_home),
    date_read_discharge_care_home = if_else(read_discharge_care_home == 1, date_read_discharge_care_home, NA)
  ) %>%
  
  #make death take precedence if it happens on the same day as any other transition
  #if death date is on the same day as a transition, make the transition indicator 0 and the date of transition NA
  mutate(discharge_home=if_else(!is.na(date_death) & !is.na(date_discharge_home) & date_death==date_discharge_home, 0, discharge_home),
         date_discharge_home=if_else(!is.na(date_death) & !is.na(date_discharge_home) & date_death==date_discharge_home, NA, date_discharge_home),
         
         discharge_care_home=if_else(!is.na(date_death) & !is.na(date_discharge_care_home) & date_death==date_discharge_care_home, 0, discharge_care_home),
         date_discharge_care_home=if_else(!is.na(date_death) & !is.na(date_discharge_care_home) & date_death==date_discharge_care_home, NA, date_discharge_care_home),
         
         readmission=if_else(!is.na(date_death) & !is.na(readmission_date) & date_death==readmission_date, 0, readmission),
         readmission_date=if_else(!is.na(date_death) & !is.na(readmission_date) & date_death==readmission_date, NA, readmission_date),
         
         read_discharge_home=if_else(!is.na(date_death) & !is.na(date_read_discharge_home) & date_death==date_read_discharge_home, 0, read_discharge_home),
         date_read_discharge_home=if_else(!is.na(date_death) & !is.na(date_read_discharge_home) & date_death==date_read_discharge_home, NA, date_read_discharge_home),
         
         read_discharge_care_home=if_else(!is.na(date_death) & !is.na(date_read_discharge_care_home) & date_death==date_read_discharge_care_home, 0, read_discharge_care_home),
         date_read_discharge_care_home=if_else(!is.na(date_death) & !is.na(date_read_discharge_care_home) & date_death==date_read_discharge_care_home, NA, date_read_discharge_care_home))


#calculate times between events
outcomes_mstate_times <- outcomes_mstate %>%
  mutate(across(c(discharge_home, discharge_care_home, readmission, read_discharge_home, read_discharge_care_home, death),
                ~ replace_na(.x, 0L)))%>% mutate(time_index=NA,
                                                 #if discharge home==1 and care home=0, time=time between index and discharge home
                                                 #if discharge home==0 and discharge care home==1, time=discharge care home
                                                 #if discharge home==0 and dicharge care home==0, time=censoring
                                                 time_dis_home=if_else(discharge_home==1, 
                                                                       as.numeric(date_discharge_home-indexdate), 
                                                                       if_else(discharge_care_home==0, 
                                                                               as.numeric(censor_date-indexdate), 
                                                                               as.numeric(date_discharge_care_home-indexdate))),
                                                 
                                                 time_dis_carehome=if_else(discharge_care_home==1, 
                                                                           as.numeric(date_discharge_care_home-indexdate),
                                                                           if_else(discharge_home==0, 
                                                                                   as.numeric(censor_date-indexdate), 
                                                                                   as.numeric(date_discharge_home-indexdate))),
                                                 
                                                 #if they were readmitted, time=time from index to readmitted, if not, censoring date
                                                 time_read=if_else(readmission==1, 
                                                                   as.numeric(readmission_date-indexdate), 
                                                                   as.numeric(censor_date-indexdate)),
                                                 
                                                 
                                                 time_read_home=if_else(read_discharge_home==1, 
                                                                        as.numeric(date_read_discharge_home-indexdate), 
                                                                        if_else(read_discharge_care_home==0, 
                                                                                as.numeric(censor_date-indexdate), 
                                                                                as.numeric(date_read_discharge_care_home-indexdate))),
                                                 
                                                 time_read_carehome=if_else(read_discharge_care_home==1, 
                                                                            as.numeric(date_read_discharge_care_home-indexdate),
                                                                            if_else(read_discharge_home==0, 
                                                                                    as.numeric(censor_date-indexdate), 
                                                                                    as.numeric(date_read_discharge_home-indexdate))),
                                                 
                                                 #if they died, time_death== death date, if not: censor date
                                                 time_death=if_else(death==1, 
                                                                    as.numeric(date_death-indexdate), 
                                                                    as.numeric(censor_date-indexdate))) 




#filter to only relevant rows of indicators and time variables
outcomes_mstate_formsprep <- outcomes_mstate_times %>% select(patid, indexdate, time_index, time_dis_home, time_dis_carehome, time_read, time_read_home, time_read_carehome, time_death,
                                                              indexhosp, discharge_home, discharge_care_home, readmission, read_discharge_home, read_discharge_care_home, death)


##----create transition matrix----
#created it manually
tmat <- matrix(NA,7,7)

dimnames(tmat) <- list(from = c("indexhosp", "discharge_home", "discharge_care_home","readmission", "read_discharge_home", "read_discharge_care_home", "death"),
                       to = c("indexhosp", "discharge_home", "discharge_care_home","readmission", "read_discharge_home", "read_discharge_care_home", "death"))

tmat[1,2:3] <- 1:2
tmat[1,7] <- 3
tmat[2,4] <- 4
tmat[2,7] <- 5
tmat[3,4] <- 6
tmat[3,7] <- 7
tmat[4,5:7] <- 8:10
tmat[5,7] <- 11
tmat[6,7] <- 12

##----use msprep()----

#input the time variables and the status (indicators)
msdata <- msprep(time = c(NA,"time_dis_home","time_dis_carehome","time_read", "time_read_home", "time_read_carehome", "time_death"),
                 status = c(NA,"discharge_home","discharge_care_home","readmission", "read_discharge_home", 
                            "read_discharge_care_home", "death"), 
                 data = outcomes_mstate_formsprep %>%
                   left_join(df %>% select(patid, indexdate, dementia, age_c, sex,
                                                      calendar_period, imd, ethnicity_5, efi2,
                                                      diabetes, ihd, depression, ckd, hypertension,
                                                      living_alone, care_home_incl_admisorc,living_status, pracid), by=c("patid","indexdate")),
                 trans = tmat,
                 keep=c("dementia", "age_c", "sex", "calendar_period", "imd", 
                        "ethnicity_5", "efi2",
                        "diabetes", "ihd", "depression", "ckd", "hypertension",
                        "living_alone", "care_home_incl_admisorc","living_status", "pracid")) 

#if there are events occuring on the same day, add 0.5 days between them
#this is for the regression to work later
#if events occur on the same day, make the Tstop 0.5 days later
#for those who transition at time=0, make the Tstop 0.5
msdata <- msdata %>% mutate(Tstop=Tstop + 
                              if_else(Tstart>0 & Tstart==Tstop, 0.5, 0)) %>%
  mutate(Tstop=Tstop + if_else(status==1 & Tstop==0, 0.5, 0))


#mutating after msprep() changes the structure, so need to add it back in

class(msdata) <- c("msdata", "data.frame")
attr(msdata, "trans") <- tmat



#create a version without missing covariate information
msdata_no_missing_cov <- msdata %>% filter(!is.na(ethnicity_5),
                                           !is.na(imd),
                                           sex!="Indeterminate")

##----transition frequencies----

# distribution of all events
events_all <- events(msdata)
events_all_no_missing_cov <- events(msdata_no_missing_cov)


#distribution of events among those with dementia
events_dementia <- events(msdata %>% filter(dementia==1))
events_dementia_no_missing_cov <- events(msdata_no_missing_cov %>% filter(dementia==1))

#distribution of events among those without dementia
events_no_dementia <- events(msdata %>% filter(dementia==0))
events_no_dementia_no_missing_cov <- events(msdata_no_missing_cov %>% filter(dementia==0))


#----Part 2: Estimation----

##----expanding covariates----

#adds the covariates to the different transitions for regression

all_covariates <- c("dementia", "age_c", "sex", "calendar_period", "imd", 
                    "ethnicity_5", "efi2",
                    "diabetes", "ihd", "depression", "ckd", "hypertension",
                    "living_alone", "care_home_incl_admisorc",
                    "living_status")

msdata_expcov <- expand.covs(msdata, all_covariates, append=TRUE, longnames = FALSE)


##----estimation----
#cox regression for each transition, comparing those with and without dementia
#also calculating robust standard errors to account for clustering by practice

###----unadjusted----
c_unadj <- coxph(Surv(Tstart, Tstop, status) ~ dementia.1 + dementia.2+
                   dementia.3 + dementia.4+ dementia.5+
                   dementia.6+ dementia.7 + dementia.8+
                   dementia.9 + dementia.10 + dementia.11 + dementia.12+ strata(as.factor(trans)) + cluster(pracid),
                 method="breslow",
                 data=msdata_expcov)


###----adjustment level 1----
#minimally adjusted for age_c, sex, calendar period
covariates1 <- c("dementia", "age_c", "sex", "calendar_period")

#adds each covariate as "dementia.1+dementia.2, etc
expanded_names_c1 <- grep(
  paste0("^(", paste(covariates1, collapse = "|"), ")[0-9]*\\.[0-9]+$"),
  names(msdata_expcov), value = TRUE
)


c1 <- coxph(as.formula(paste("Surv(Tstart, Tstop, status) ~",
                             paste(expanded_names_c1, collapse = " + "),
                             "+ strata(as.factor(trans)) + cluster(pracid)")), 
            method = "breslow", 
            data = msdata_expcov)


###----adjustment level 2----

#additionally adjusted for deprivation/ethnicity
covariates2 <- c("dementia", "age_c", "sex", "calendar_period", "imd", 
                 "ethnicity_5")

expanded_names_c2 <- grep(
  paste0("^(", paste(covariates2, collapse = "|"), ")[0-9]*\\.[0-9]+$"),
  names(msdata_expcov), value = TRUE)


c2 <- coxph(as.formula(paste("Surv(Tstart, Tstop, status) ~",
                             paste(expanded_names_c2, collapse = " + "),
                             "+ strata(as.factor(trans))+cluster(pracid)")), 
            method = "breslow", 
            data = msdata_expcov)


###----adjustment level 3----

#additionally adjusted for efi2 score

covariates3 <- c("dementia", "age_c", "sex", "calendar_period", "imd", 
                 "ethnicity_5", "efi2")

expanded_names_c3 <- grep(
  paste0("^(", paste(covariates3, collapse = "|"), ")[A-Za-z0-9]*\\.[0-9]+$"),
  names(msdata_expcov), value = TRUE)
keep_cols <- c("Tstart", "Tstop", "status", "trans", "pracid", expanded_names_c3)
msdata_c3 <- msdata_expcov %>% select(all_of(keep_cols))
rm(msdata_expcov)

c3 <- coxph(as.formula(paste("Surv(Tstart, Tstop, status) ~",
                             paste(expanded_names_c3, collapse = " + "),
                             "+ strata(as.factor(trans)) +cluster(pracid)")), 
            method = "breslow", 
            data = msdata_c3,
            model  = FALSE,   # don't store the model frame
            x      = FALSE,
            y      = FALSE,
            robust = TRUE     # still gives cluster-robust var via cluster(pracid)
)



tidy(c3, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high)

###----fully adjusted----

#fully adjusted model needs additional adjustment for hypertension, ckd, ihd, diabetes, depression

#to save computing power and not run out of memory, need to calculate the hazard ratios for each transition separately
#this gives the same transition-specific HRs as running one large model, as instead of stratifying the model by transitions, I am just running each transition separately
#don't need expanded covariates for this version


covariates4 <- c("dementia", "age_c", "sex", "calendar_period", "imd",
                 "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                 "ckd", "hypertension")


#create a function to fit a cox model for each transition
fit_one_trans <- function(k, data, covs) {
  print(k)
  #filter the data by transition
  sub <- data[data$trans == k, ]
  f <- as.formula(paste("Surv(Tstart, Tstop, status) ~",
                        paste(covs, collapse = " + "),
                        "+ cluster(pracid)"))
  fit <- coxph(f, data = sub, method = "breslow",
               x = FALSE, y = FALSE, model = FALSE)
  rm(sub, fit); gc()
}

#run the function for each transition
for (k in seq_len(12)) fit_one_trans(k, msdata, covariates4)



##----Figure 2 - forest plot----

#takes a list of cox models and pulls out the HR (95% CI)
#uses the event counts from the events() function

create_fp_c4 <- function(models_list, trans_order, term_name) {
  
  # Pull the single row for `term_name` out of each transition's model,
  # in the order given by trans_order, and stack them into one data frame
  plot_df <- map_dfr(trans_order, function(k) {
    tidy(models_list[[as.character(k)]], exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == term_name)
  }) |>
    select(term, estimate, conf.low, conf.high) |>
    # Add n. events/at risk number of events and N
    mutate(events_exp = c(events_dementia_no_missing_cov$Frequencies[1,2],events_dementia_no_missing_cov$Frequencies[1,3],
                          events_dementia_no_missing_cov$Frequencies[1,7],events_dementia_no_missing_cov$Frequencies[2,4],
                          events_dementia_no_missing_cov$Frequencies[2,7],events_dementia_no_missing_cov$Frequencies[3,4],
                          events_dementia_no_missing_cov$Frequencies[3,7],events_dementia_no_missing_cov$Frequencies[4,5],
                          events_dementia_no_missing_cov$Frequencies[4,6],events_dementia_no_missing_cov$Frequencies[4,7],
                          events_dementia_no_missing_cov$Frequencies[5,7],events_dementia_no_missing_cov$Frequencies[6,7]),
           sample_size_exp = c(rep(events_dementia_no_missing_cov$Frequencies[1,9],3),
                               rep(events_dementia_no_missing_cov$Frequencies[2,9],2),
                               rep(events_dementia_no_missing_cov$Frequencies[3,9],2),
                               rep(events_dementia_no_missing_cov$Frequencies[4,9],3),
                               events_dementia_no_missing_cov$Frequencies[5,9],
                               events_dementia_no_missing_cov$Frequencies[6,9])) |>
    mutate(events_un_exp = c(events_no_dementia_no_missing_cov$Frequencies[1,2],events_no_dementia_no_missing_cov$Frequencies[1,3],
                             events_no_dementia_no_missing_cov$Frequencies[1,7],events_no_dementia_no_missing_cov$Frequencies[2,4],
                             events_no_dementia_no_missing_cov$Frequencies[2,7],events_no_dementia_no_missing_cov$Frequencies[3,4],
                             events_no_dementia_no_missing_cov$Frequencies[3,7],events_no_dementia_no_missing_cov$Frequencies[4,5],
                             events_no_dementia_no_missing_cov$Frequencies[4,6],events_no_dementia_no_missing_cov$Frequencies[4,7],
                             events_no_dementia_no_missing_cov$Frequencies[5,7],events_no_dementia_no_missing_cov$Frequencies[6,7]),
           sample_size_un_exp = c(rep(events_no_dementia_no_missing_cov$Frequencies[1,9],3),
                                  rep(events_no_dementia_no_missing_cov$Frequencies[2,9],2),
                                  rep(events_no_dementia_no_missing_cov$Frequencies[3,9],2),
                                  rep(events_no_dementia_no_missing_cov$Frequencies[4,9],3),
                                  events_no_dementia_no_missing_cov$Frequencies[5,9],
                                  events_no_dementia_no_missing_cov$Frequencies[6,9])) |>
    mutate(n_N_exp = paste0(events_exp,"/",sample_size_exp),
           n_N_un_exp = paste0(events_un_exp,"/",sample_size_un_exp)) |>
    mutate(n_N_exp = paste0(if_else(events_exp    %in% 1:50, "<50", format(events_exp, big.mark = ",")),
                            "/", format(sample_size_exp, big.mark = ",")),
           n_N_un_exp = paste0(if_else(events_un_exp %in% 1:50, "<50", format(events_un_exp, big.mark = ",")),
                               "/", format(sample_size_un_exp, big.mark = ","))) |>
    mutate(term_description = c("Usual residence","New care home",
                                "Death","Usual residence → Readmitted",
                                "Usual residence → Death","New care home → Readmitted",
                                "New care home → Death","Readmitted → Usual residence",
                                "Readmitted → New care home","Readmitted → Death",
                                "Readmitted → Usual residence → Death", "Readmitted → New care home → Death"),.after = 1) |>
    mutate(term_description = factor(term_description,
                                     levels = c("Usual residence","New care home",
                                                "Death","Usual residence → Readmitted",
                                                "Usual residence → Death","New care home → Readmitted",
                                                "New care home → Death","Readmitted → Usual residence",
                                                "Readmitted → New care home","Readmitted → Death",
                                                "Readmitted → Usual residence → Death", "Readmitted → New care home → Death"))) |>
    mutate(across(c(conf.low,conf.high),~if_else(.x < 0.0001 | .x > 100,0,.x))) |>
    mutate(est_label = if_else(is.na(estimate),NA_character_,
                               paste0(formatC(estimate,format = "f",digits = 2),
                                      " (",formatC(conf.low,format = "f",digits = 2),", ",
                                      formatC(conf.high,format = "f",digits = 2),")")),
           across(c(estimate,conf.low,conf.high),~round(.x,2)) ) |>
    mutate(across(c(estimate,conf.low,conf.high),~if_else(.x < 0.0001 | .x > 100,NA,.x))) |>
    mutate(est_label_middot=paste0(str_replace_all(est_label,"[.]","\U00B7")),
           estimate_middot = str_replace(estimate,"[.]","\U00B7")) |>
    mutate(across(c(conf.low,conf.high,est_label,estimate),~ if_else(events_exp==0|events_un_exp==0,NA,.x)))
  
  fp <- forestplot(plot_df,
                   labeltext=c(term_description,n_N_exp,n_N_un_exp,est_label),
                   mean = estimate,
                   lower = conf.low,
                   upper = conf.high,
                   zero=1,
                   xlab      = "Hazard Ratio (95% CI)", 
                   xticks= c(0.7, 1, 1.4, 2),
                   txt_gp   = fpTxtGp(
                     ticks = gpar(cex = 1),                   
                     xlab  = gpar(cex = 1),
                     legend=gpar(cex=0.85)
                   ),
                   boxsize=0.15,
                   xlog=TRUE,
                   colgap=unit(8, "mm")) %>% #added title below
    fp_add_header(n_N_exp = "Dementia\n(n. events/at risk)",
                  n_N_un_exp = "No dementia\n(n. events/at risk)",
                  est_label= "Adjusted HR (95% CI)",
                  term_description = "Transition from index hospitalisation:") |>
    fp_set_style(line = "black",
                 box = "black",
                 zero =gpar(col = "black", lty = "solid")) |>
    fp_decorate_graph(graph.pos = 4)
  print(fp)
  
  
}

#list of models for each transition, fully adjusted

#models_list <- list of the transition-specific HRs

create_fp_c4(models_list, trans_order = 1:12, term_name = "dementia")



#----Part 3 - Prediction----

#first need to use estimated parameters and baseline transition hazards and covariate values to obtain patient-specific transition hazards
#create a new dataset containing expanded covariates, one line per transition

#using the model adjusted for age, sex, calendar period, ethnicity, IMD, frailty
#this is the last model I was able to run as one cox model, rather than transition-specific models
#fully-adjusted model failed when running as one model, as it had too many parameters
#msfit() and downstream functions require one cox model, which is why I'm using this model
#this model vs the fully-adjusted model gave nearly identical HRs, so I am comfortable using this model for prediction

covariates_for_pred <- c("dementia", "age_c", "sex", "calendar_period",
                         "ethnicity_5", "imd", "efi2")

# un_exp (no dementia)
#chose covariates for referent person
newdata_un_exp <- data.frame(
  trans           = 1:12,
  strata          = 1:12,
  dementia        = 0,
  age_c           = 0,
  sex             = factor("Female", levels = levels(df$sex)),
  calendar_period = factor("2018-2019", levels = levels(df$calendar_period)),
  ethnicity_5 = factor("White", levels=levels(as.factor(df$ethnicity_5))),
  imd = factor("5", levels=levels(df$imd)),
  efi2=median(df$efi2))


attr(newdata_un_exp,"trans") <- tmat
class(newdata_un_exp) <- c("msdata","data.frame")
newdata_un_exp <- expand.covs(newdata_un_exp,covariates_for_pred,longnames=F)


# exposed (dementia)
newdata_exp <- data.frame(
  trans           = 1:12,
  strata          = 1:12,
  dementia        = 1, #dementia=1 for exposed group
  age_c= 0,
  sex             = factor("Female", levels = levels(df$sex)),
  calendar_period = factor("2018-2019", levels = levels(df$calendar_period)),
  ethnicity_5 = factor("White", levels=levels(as.factor(df$ethnicity_5))),
  imd = factor("5", levels=levels(df$imd)),
  efi2=median(df$efi2))

attr(newdata_exp,"trans") <- tmat
class(newdata_exp) <- c("msdata","data.frame")
newdata_exp <- expand.covs(newdata_exp,covariates_for_pred,longnames=F)



#use msfit function to obtain patient-specific transition hazards
#the patient of interest is someone with the referent covariate values, and either dementia=0 or dementia=1
#haz=estimated baseline cumulative hazard for each of the transitions at each time


#c3 is the cox model adjusted for age, sex, calendar period, ethnicity, imd, frailty
msf_un_exp <- msfit(c3, newdata=newdata_un_exp, trans=tmat)
msf_exp <- msfit(c3, newdata=newdata_exp, trans=tmat)


#now calculate the probability of each transition
#predt= time to start predictions (time 0)
#probtrans() uses the baseline hazards (from msfit) to predict the patient-specific transition probabilities

probtrans_un_exp <- probtrans(msf_un_exp,predt = 0)
probtrans_exp <- probtrans(msf_exp,predt = 0)



##----probability of occupying each state at 6 months----
#everyone starts at initial state [[1]]
#calculating probabilities of being in each state at 183 days

state_occu_prob_6_mo <- rbind(
  slice(filter(probtrans_exp[[1]],time>=183),1),
  slice(filter(probtrans_un_exp[[1]],time>=183),1))|>
  rename(Index_hosp = pstate1,
         Usual_res = pstate2,
         Care_home = pstate3,
         Readmission = pstate4,
         Read_usual_res = pstate5,
         Read_care_home = pstate6,
         Death = pstate7,
         Index_hosp_se = se1,
         Usual_res_se = se2,
         Care_home_se = se3,
         Readmission_se = se4,
         Read_usual_res_se = se5,
         Read_care_home_se = se6,
         Death_se = se7,) |>
  mutate(across(Index_hosp:Death_se,~round(.x,3)*100)) |>
  mutate(dementia=c(1,0)) %>%

  # calculate 95% CI
  mutate(across(matches("(Index_hosp|Usual_res|Care_home|Readmission|Read_usual_res|Read_care_home|Death)$"),
                ~paste0(.x," (",.x-round(1.96*get(glue("{cur_column()}_se")),1),", ",
                        .x+round(1.96*get(glue("{cur_column()}_se")),1),")"),
                .names = "{.col}_ci")) |>
  select(time,dementia,matches("ci")) |>
  mutate(across(matches("ci"),~if_else(.x=="0 (0-0)",NA_character_,.x)))

print(state_occu_prob_6_mo)

state_occu_prob_6_mo |>
  select(-c(time)) |>
  gt(caption = "Probability of transition at 5 years")


rbind(
  slice(filter(probtrans_exp[[1]],time>=183),1),
  slice(filter(probtrans_un_exp[[1]],time>=183),1))|>
  mutate(Hospital = pstate1+pstate4,
         Usual_res = pstate2+pstate5,
         Care_home = pstate3+pstate6,
         Death = pstate7) |>
  mutate(dementia=c(1,0)) %>%
  select(time, Hospital, Usual_res, Care_home, Death, dementia)

##----stacked occupancy plots----


#add the states of hospital, discharge to home, discharge to care home
modified_un_exp <- probtrans_un_exp[[1]] %>%
  mutate(hospital=pstate1+pstate4,
         home=pstate2+pstate5,
         care_home=pstate3+pstate6,
         death=pstate7) %>%
  select(time,hospital,home,care_home,death) %>%
  pivot_longer(cols = c(hospital, home, care_home, death),
               names_to = "state",
               values_to = "prob") %>%
  mutate(state = factor(state,
                        levels = c("death", "care_home", "home", "hospital"),
                        labels = c("Death", "New care home", "Usual residence", "Hospital")))

modified_exp <- probtrans_exp[[1]] %>%
  mutate(hospital=pstate1+pstate4,
         home=pstate2+pstate5,
         care_home=pstate3+pstate6,
         death=pstate7) %>%
  select(time,hospital,home,care_home,death) %>%
  pivot_longer(cols = c(hospital, home, care_home, death),
               names_to = "state",
               values_to = "prob") %>%
  mutate(state = factor(state,
                        levels = c("death", "care_home", "home", "hospital"),
                        labels = c("Death", "New care home", "Usual residence", "Hospital")))

#create a function to plot the stacked occupancy probability plots
plot_collapsed_categories <- function(data, title) {
  
  ggplot(data, aes(x = time, y = prob, fill = state)) +
    geom_area(linewidth = 0.75, colour = "white") +
    scale_fill_viridis_d(direction = -1,
                         breaks = c("Hospital", "Usual residence", "New care home", "Death")) +
    labs(title = title,
         x = "Time since index date (days)",
         y = "Probability",
         fill = NULL) +
    coord_cartesian(xlim = c(0, 365)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0, size = 20),
          plot.title.position = "plot",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.position = "bottom",
          panel.grid = element_blank()) +
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))
}

collapsed_stacked_plot_dementia<- plot_collapsed_categories(modified_exp, "A. Dementia")
collapsed_stacked_plot_no_dementia <- plot_collapsed_categories(modified_un_exp, "B. No dementia")



##----expected length of stay in each state----

#ELOS() function uses the probtrans() object to estimate number of days in each state
#using the diagonal values gives the estimate CONDITIONAL on having reached that state

elos_exp <- round(ELOS(probtrans_exp),1)
elos_un_exp <- round(ELOS(probtrans_un_exp),1)


#----Part 4. Stratified analyses----
#assessing effect modification by stratifying the analysis by sex, age group, living status


#function to fit adjusted cox models for each transition AND each subset
#data= msdata subsetted by sex, age, etc
#k=transitions (12)
#covs=list of covariates to adjust for

fit_one_trans_subset <- function(k, data, covs, subset) {
  print(k)
  sub <- data[data$trans == k, ]
  f <- as.formula(paste("Surv(Tstart, Tstop, status) ~",
                        paste(covs, collapse = " + "),
                        "+ cluster(pracid)"))
  fit <- coxph(f, data = sub, method = "breslow",
               x = FALSE, y = FALSE, model = FALSE)
  rm(sub, fit); gc()
}

##----sex----

#list of covariates excluding sex
covariates_exc_sex <- c("dementia", "age_c", "calendar_period", "imd",
                        "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                        "ckd", "hypertension")

#msdata subsetted by sex

msdata_male<- msdata %>%
  filter(sex=="Male")

msdata_female<- msdata %>%
  filter(sex=="Female")


#run the cox models on each transition per subset
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_male, covariates_exc_sex, subset="male")
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_female, covariates_exc_sex, subset="female")

#create a list of models by sex which will go into the forest plot function

#models_list_male <- list of models for males
#models_list_female <- list of models for females

###----create forest plots for EM----

#function to pull out the HR and CI for each model list
#models_list = list of models
#trans_order = 1:12
#term_name=dementia
#trans_labels = labels for each transition in the forest plot

build_subset_df <- function(models_list, trans_order, term_name, trans_labels) {
  
  map_dfr(trans_order, function(k) {
    mod <- models_list[[as.character(k)]]
    
    if (is.null(mod) || !inherits(mod, "coxph")) {
      return(tibble(term = term_name,
                    estimate = NA_real_,
                    conf.low = NA_real_,
                    conf.high = NA_real_))
    }
    
    result <- tidy(mod, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == term_name)
    
    if (nrow(result) == 0) {
      return(tibble(term = term_name,
                    estimate = NA_real_,
                    conf.low = NA_real_,
                    conf.high = NA_real_))
    }
    
    result
  }) |>
    select(term, estimate, conf.low, conf.high) |>
    mutate(term_description = factor(trans_labels[trans_order],
                                     levels = trans_labels[trans_order])) |>
    mutate(across(c(conf.low, conf.high), ~ if_else(.x < 0.0001 | .x > 100, 0, .x))) |>
    mutate(
      est_label = if_else(is.na(estimate), NA_character_,
                          paste0(formatC(estimate, format = "f", digits = 2),
                                 " (", formatC(conf.low, format = "f", digits = 2), ", ",
                                 formatC(conf.high, format = "f", digits = 2), ")")),
      across(c(estimate, conf.low, conf.high), ~ round(.x, 2))
    ) |>
    mutate(across(c(estimate, conf.low, conf.high), ~ if_else(.x < 0.0001 | .x > 100, NA, .x)))
}

#create forest plots for EM, by subset

create_fp_c4_subset <- function(models_lists,
                                trans_order,
                                term_name,
                                title,
                                group_names,
                                group_colors = NULL,
                                graph.pos    = 2,
                                strata,
                                trans_labels) {
  
  n_groups <- length(models_lists)
  if (is.null(group_names)) group_names <- paste0("Group", seq_len(n_groups))
  
  if (is.null(group_colors)) {
    jama_pal  <- ggsci::pal_jama("default")(7)
    group_colors <- jama_pal[seq_len(n_groups)]
  }
  
  subset_dfs <- map(models_lists, function(ml) build_subset_df(ml, trans_order, term_name, trans_labels))
  
  mean_mat  <- do.call(cbind, map(subset_dfs, "estimate"))
  lower_mat <- do.call(cbind, map(subset_dfs, "conf.low"))
  upper_mat <- do.call(cbind, map(subset_dfs, "conf.high"))
  colnames(mean_mat) <- colnames(lower_mat) <- colnames(upper_mat) <- group_names
  
  labeltext <- list(term_description = trans_labels[trans_order])
  for (i in seq_len(n_groups)) {
    labeltext[[group_names[i]]] <- subset_dfs[[i]]$est_label
  }
  
  hr_headers <- setNames(
    lapply(seq_len(n_groups), function(i) group_names[i]),
    group_names
  )
  
  total_cols <- n_groups + 2
  hrzl_lines <- setNames(
    lapply(2:13, function(i) gpar(lwd = 4, columns = 1:total_cols, col = "#99999922")),
    as.character(2:13)
  )
  
  fp <- forestplot(
    labeltext = labeltext,
    mean      = mean_mat,
    lower     = lower_mat,
    upper     = upper_mat,
    zero      = 1,
    legend    = group_names,
    title     = title,
    boxsize   = 0.15,
    xlog      = TRUE,
    xticks    = c(0.7, 1, 1.4, 2),
    xlab      = "Hazard Ratio (95% CI)",   
    align     = c("l", rep("c", n_groups)), 
    txt_gp    = fpTxtGp(
      ticks  = gpar(cex = 1),                   
      xlab   = gpar(cex = 1),
      legend = gpar(cex = 1)
    ),
    graphwidth = unit(200, "mm"),
    hrzl_lines = hrzl_lines
  ) |>
    fp_set_style(
      box  = group_colors,
      line = group_colors,
      zero = gpar(col = "black", lty = "solid")
    ) |>
    fp_decorate_graph(graph.pos = graph.pos)   # <- removed the trailing fp_add_header here
  
  fp <- do.call(
    fp_add_header,
    c(list(fp, term_description = "Transition from index hospitalisation:"), hr_headers)
  )
  
  fp
}

transition_labels <- c(
  "Usual residence","New care home",
  "Death","Usual residence → Readmitted",
  "Usual residence → Death","New care home → Readmitted",
  "New care home → Death","Readmitted → Usual residence",
  "Readmitted → New care home","Readmitted → Death",
  "Readmitted → Usual residence → Death", "Readmitted → New care home → Death"
)

#run the function to create the forest plot
create_fp_c4_subset(
  models_list=list(Men=models_list_male, Women=models_list_female),
  trans_order=1:12,
  term_name="dementia",
  title="Figure 2. Adjusted hazard ratios for each transition by sex",
  group_names  = c("Men", "Women"),
  strata="sex",
  trans_labels=transition_labels)


##----age group----

#create age groups
msdata <- msdata %>% mutate(age=age_c+mean(df$age),
                            age_group=case_when(age <70 ~ "65-69",
                                                age >=70 & age <75 ~ "70-74",
                                                age >=75 & age <80 ~ "75-79",
                                                age >=80 & age <85 ~ "80-84",
                                                age >=85 & age <90 ~ "85-89",
                                                age >=90 ~ "90+"))

covariates_exc_age <- c("dementia", "sex", "calendar_period", "imd",
                        "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                        "ckd", "hypertension")

msdata_age1<- msdata %>% filter(age_group=="65-69")
msdata_age2<- msdata %>% filter(age_group=="70-74")
msdata_age3<- msdata %>% filter(age_group=="75-79")
msdata_age4<- msdata %>% filter(age_group=="80-84")
msdata_age5<- msdata %>% filter(age_group=="85-89")
msdata_age6<- msdata %>% filter(age_group=="90+")

#run the cox models on each transition per subset
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_age1, covariates_exc_age, subset="age1")
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_age2, covariates_exc_age, subset="age2")
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_age3, covariates_exc_age, subset="age3")
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_age4, covariates_exc_age, subset="age4")
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_age5, covariates_exc_age, subset="age5")
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_age6, covariates_exc_age, subset="age6")

#models_list_age1 <- list
#models_list_age2 <- list
#models_list_age3 <- list
#models_list_age4 <- list
#models_list_age5 <- list
#models_list_age6 <- list

#create a table of HR by age group for each transition


# extract HR 95% CI for one model
extract_hr_ci <- function(model, term) {
  if (is.null(model)) return(NA_character_)
  tt <- tryCatch(tidy(model, exponentiate = TRUE, conf.int = TRUE),
                 error = function(e) NULL)
  if (is.null(tt)) return(NA_character_)
  row <- tt %>% filter(term == !!term)
  if (nrow(row) == 0) return(NA_character_)   
  sprintf("%.2f (%.2f-%.2f)", row$estimate, row$conf.low, row$conf.high)
}

# find all models_list_ageN objects in the environment, sorted numerically
age_list_names <- ls(envir = .GlobalEnv, pattern = "^models_list_age[0-9]+$")
age_list_names <- age_list_names[order(as.numeric(str_extract(age_list_names, "[0-9]+")))]


# build one column per age group
hr_table <- map_dfc(age_list_names, function(nm) {
  model_list <- get(nm, envir = .GlobalEnv)
  col <- map_chr(model_list, extract_hr_ci, term = "dementia")
  tibble(!!nm := col)
})

hr_table <- bind_cols(Transition = transition_labels, hr_table) 
hr_table_age <- hr_table %>% 
  rename("Age 65-69"=models_list_age1,
         "Age 70-74"=models_list_age2,
         "Age 75-79"=models_list_age3,
         "Age 80-84"=models_list_age4,
         "Age 85-89"=models_list_age5,
         "Age 90+"=models_list_age6)


##----living status----

covariates <- c("dementia","sex", "age_c", "calendar_period", "imd",
                "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                "ckd", "hypertension")

###----living alone----
msdata_livingalone<- msdata %>%
  filter(living_alone==1)

#run the cox models on each transition per subset
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_livingalone, covariates, subset="livingalone")

#create a list of models by sex which will go into the forest plot function
#models_list_livingalone <- list


###----care home----

msdata_carehome<- msdata %>%
  filter(care_home_incl_admisorc==1)

#only showing possible transitions for care home
k_carehome <- c(1,2,3, 4,5, 8, 10, 11)

#run the cox models on each transition per subset
for (k in k_carehome) fit_one_trans_subset(k, msdata_carehome, covariates, subset="carehome")

#models_list_carehome <- list


###----community----
#combining living alone and living in a care home into one plot using living_status_variable

#creating a subset of those labelled as 'community dwelling'
msdata_community<- msdata %>%
  filter(living_status=="Community-dwelling with others")

#run the cox models on each transition per subset
for (k in seq_len(12)) fit_one_trans_subset(k, msdata_community, covariates, subset="community")

#models_list_community <- list

create_fp_c4_living_status <- function(models_lists,
                                trans_order,
                                term_name,
                                title,
                                group_names,
                                group_colors = NULL,
                                graph.pos    = 2,
                                strata,
                                trans_labels) {
  
  n_groups <- length(models_lists)
  if (is.null(group_names)) group_names <- paste0("Group", seq_len(n_groups))
  
  if (is.null(group_colors)) {
    jama_pal  <- ggsci::pal_jama("default")(7)
    group_colors <- jama_pal[seq_len(n_groups)]
  }
  
  subset_dfs <- map(models_lists, function(ml) build_subset_df(ml, trans_order, term_name, trans_labels))
  
  mean_mat  <- do.call(cbind, map(subset_dfs, "estimate"))
  lower_mat <- do.call(cbind, map(subset_dfs, "conf.low"))
  upper_mat <- do.call(cbind, map(subset_dfs, "conf.high"))
  colnames(mean_mat) <- colnames(lower_mat) <- colnames(upper_mat) <- group_names
  
  labeltext <- list(term_description = trans_labels[trans_order])

  hrzl_lines <- setNames(
    lapply(2:13, function(i) gpar(lwd = 4, columns = 1:2, col = "#99999922")),
    as.character(2:13)
  )
  
  fp <- forestplot(
    labeltext = labeltext,
    mean      = mean_mat,
    lower     = lower_mat,
    upper     = upper_mat,
    zero      = 1,
    legend    = group_names,
    title     = title,
    boxsize   = 0.15,
    xlog      = TRUE,
    xticks    = c(0.7, 1, 1.4, 2, 2.3),
    xlab      = "Hazard Ratio (95% CI)",   
    align     = c("l", rep("c", n_groups)), 
    txt_gp    = fpTxtGp(
      ticks  = gpar(cex = 1),                   
      xlab   = gpar(cex = 1),
      legend = gpar(cex = 1)
    ),
    graphwidth = unit(200, "mm"),
    hrzl_lines = hrzl_lines
  ) |>
    fp_set_style(
      box  = group_colors,
      line = group_colors,
      zero = gpar(col = "black", lty = "solid")
    ) |>
    fp_decorate_graph(graph.pos = graph.pos) 
  
  fp <- do.call(
    fp_add_header,
    c(list(fp, term_description = "Transition from index hospitalisation:"))
  )
  
  fp
}


create_fp_c4_living_status(
  models_list=list(livingalone=models_list_livingalone, carehome=models_list_carehome, community=models_list_community),
  trans_order=1:12,
  term_name="dementia",
  title="Adjusted hazard ratios for each transition by living status",
  group_names=c("Living alone", "Living in a care home", "Community-dwelling with others"),
  strata="living_status",
  trans_labels=transition_labels)


#----Part 5 - sequential adjustment plot----

#creating a plot to show the the sequentially adjusted models
#c_unadj, c1, c2, c3 are all one big model with transitions incorporated
#c4 (fully adjusted model) is a list of transition-specific models


#function to extract HR and 95% CI from the models where all transitions are included in a single odel

extract_single_model <- function(mod, trans_order, term_prefix) {
  map_dfr(trans_order, function(k) {
    term_name <- paste0(term_prefix, ".", k)
    
    if (is.null(mod) || !inherits(mod, "coxph")) {
      return(tibble(term = term_name, estimate = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_))
    }
    
    result <- tidy(mod, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == term_name)
    
    if (nrow(result) == 0) {
      return(tibble(term = term_name, estimate = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_))
    }
    
    result |> select(term, estimate, conf.low, conf.high)
  })
}

#function to extract HR and CI from the models list
extract_model_list <- function(models_list, trans_order, term_name) {
  map_dfr(trans_order, function(k) {
    mod <- models_list[[as.character(k)]]
    
    if (is.null(mod) || !inherits(mod, "coxph")) {
      return(tibble(term = term_name, estimate = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_))
    }
    
    result <- tidy(mod, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == term_name)
    
    if (nrow(result) == 0) {
      return(tibble(term = term_name, estimate = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_))
    }
    
    result |> select(term, estimate, conf.low, conf.high)
  })
}

#function to combine the above model information
build_subset_df <- function(source, trans_order, trans_labels) {
  
  df <- switch(source$type,
               single = extract_single_model(source$model, trans_order, source$term_prefix),
               list   = extract_model_list(source$models, trans_order, source$term_name),
               stop("Unknown source type: ", source$type)
  )
  
  df |>
    mutate(term_description = factor(trans_labels[trans_order],
                                     levels = trans_labels[trans_order])) |>
    mutate(across(c(conf.low, conf.high), ~ if_else(.x < 0.0001 | .x > 100, 0, .x))) |>
    mutate(
      est_label = if_else(is.na(estimate), NA_character_,
                          paste0(formatC(estimate, format = "f", digits = 2),
                                 " (", formatC(conf.low, format = "f", digits = 2), ", ",
                                 formatC(conf.high, format = "f", digits = 2), ")")),
      across(c(estimate, conf.low, conf.high), ~ round(.x, 2))
    ) |>
    mutate(across(c(estimate, conf.low, conf.high), ~ if_else(.x < 0.0001 | .x > 100, NA, .x)))
}

create_fp_adjustment_subset <- function(sources,          # named list of source specs
                                        trans_order,
                                        title,
                                        group_colors = NULL,
                                        graph.pos = 2,  
                                        strata,
                                        trans_labels) {
  
  group_names <- names(sources)
  n_groups <- length(sources)
  
  if (is.null(group_colors)) {
    jama_pal  <- ggsci::pal_jama("default")(7)
    group_colors <- jama_pal[seq_len(n_groups)]
  }
  
  subset_dfs <- map(sources, ~ build_subset_df(.x, trans_order, trans_labels))
  
  mean_mat  <- do.call(cbind, map(subset_dfs, "estimate"))
  lower_mat <- do.call(cbind, map(subset_dfs, "conf.low"))
  upper_mat <- do.call(cbind, map(subset_dfs, "conf.high"))
  colnames(mean_mat) <- colnames(lower_mat) <- colnames(upper_mat) <- group_names
  
  labeltext <- list(term_description = trans_labels[trans_order])
  
  hrzl <- setNames(
    map(2:13, ~ gpar(lwd = 4, columns = 1:2, col = "#99999922")),
    as.character(2:13)
  )
  
  fp <- forestplot(
    labeltext = labeltext,
    mean      = mean_mat,
    lower     = lower_mat,
    upper     = upper_mat,
    zero      = 1,
    legend    = group_names,
    title     = title,
    boxsize   = 0.15,
    xlog      = TRUE,
    xticks    = c(0.67, 1, 2, 3),
    xlab      = "Hazard Ratio (95% CI)",        
    txt_gp   = fpTxtGp(
      ticks = gpar(cex = 1),                   
      xlab  = gpar(cex = 1),
      legend=gpar(cex=0.85)
    ),
    graphwidth = unit(200, "mm"),
    hrzl_lines = hrzl
  ) |>
    fp_set_style(
      box  = group_colors,
      line = group_colors,
      zero = gpar(col = "black", lty = "solid")
    ) |>
    fp_decorate_graph(graph.pos = graph.pos) |>
    fp_add_header(term_description = "Transition from index hospitalisation:")
  

  fp
}


sources <- list(
  "Unadjusted"      = list(type = "single", model = c_unadj, term_prefix = "dementia"),
  "Age, sex, calendar period"         = list(type = "single", model = c1,      term_prefix = "dementia"),
  "+ Ethnicity, deprivation"         = list(type = "single", model = c2,      term_prefix = "dementia"),
  "+ eFI2 frailty score"         = list(type = "single", model = c3,      term_prefix = "dementia"),
  "+ Additional conditions* (Fully adjusted)"  = list(type = "list",   models = models_list, term_name = "dementia")
)

transition_labels <- c(
  "Usual residence","New care home",
  "Death","Usual residence → Readmitted",
  "Usual residence → Death","New care home → Readmitted",
  "New care home → Death","Readmitted → Usual residence",
  "Readmitted → New care home","Readmitted → Death",
  "Readmitted → Usual residence → Death", "Readmitted → New care home → Death"
)

seq_adjustment_plot <- create_fp_adjustment_subset(
  sources      = sources,
  trans_order  = 1:12,
  title        = "Sequential adjustment",
  trans_labels = transition_labels
)

