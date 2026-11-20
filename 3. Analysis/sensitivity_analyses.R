
#--Monica Iyer
#--August 2026
#--Analysis for care pathways following infections study

#This script contains code for the sensitivity analyses@

#Part 1: Set up data for sensitivity analyses:
#--End follow-up on March 1, 2020 (Figure S4)
#--Restrict cohort to individuals' first infection-related hospitalisation (Figure S5)
#--Only include hospital admissions where infection is in first episode (Figure S9)
#--Restrict to those who were not previously in a care home (Figure S10)
#--Additional covariate adjustment (Figure S11)
#Part 2: Run mstate analysis pipeline on the data for sensitivity analyses and plot results
#Part 3: Run sensitivity analysis with clock-reset timescale (Figure S8)
#Part 4: Run sensitivity analysis where we allow up to 3 readmissions in the model (Figures S6-S7)




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

#----read in parquet files for the baseline variables and outcomes----

##----read in data----

#read in df and outcomes that were cleaned in descriptive_analysis script

#df has patid, indexdate, and all baseline covariates
#outcomes has event-level information
#patients can enter the cohort multiple times; each patid/indexdate combination is unique

df <- read_parquet("")
outcomes <- read_parquet("")

#----Part 1: create dfs for sensitivity analyses----


##----end followup in march 1 2020----

#restrict follow-up to March 1, 2020
#only include indexdates on or before that date;
#if follow-up ends after that date, say it ends on that date

sens_fu_df <- df %>% filter(indexdate<= as.Date("2020-03-01")) %>%
  mutate(observation_end=if_else(observation_end>as.Date("2020-03-01"), as.Date("2020-03-01"), observation_end))

#create outcomes data for sensitivity analysis
sens_fu_outcomes <- outcomes %>% 
  inner_join(sens_fu_df %>% select(patid, indexdate), by=c("patid", "indexdate")) %>%
  filter(spell_start<=as.Date("2020-03-01"))

##----first infection----

#limit to cohort to first admission for infection
sens_firstinfec_df <- df %>% arrange(patid, indexdate) %>% distinct(patid, .keep_all = TRUE)
sens_firstinfec_outcomes <- outcomes %>% 
  inner_join(sens_firstinfec_df %>% select(patid, indexdate), by=c("patid", "indexdate"))

##----only include infections within first episode----

#needs a separate pipeline: study_eligibility_criteria_sens --> _targets to create variables, then read in
sens_infecfirstep_df <- read_parquet("")
sens_infecfirstep_df <- sens_infecfirstep_df %>% distinct() %>%
  mutate(ethnicity_5=case_when(ethnicity_5=="0. White" ~ "White",
                             ethnicity_5=="1. South Asian" ~ "South Asian",
                             ethnicity_5=="2. Black" ~ "Black",
                             ethnicity_5=="3. Other" ~ "Other",
                             ethnicity_5=="4. Mixed" ~ "Mixed",
                             ethnicity_5=="5. Not Stated" ~ NA)) %>%
  rename(efi2=efi2_score_excl_dementia) %>%
  mutate(imd=as.factor(imd)) %>%
  subset(sex!="Indeterminate") %>%
  mutate(sex=as.factor(sex)) %>%
  mutate(calendar_period=as.factor(calendar_period)) %>%
  mutate(age_c = age - mean(age, na.rm=TRUE))

outcomes_sens <- read_parquet("")
sens_infecfirstep_outcomes <- outcomes_sens %>% distinct()


##----restrict to those without a prior code for care home----

sens_nocarehome_df <- df %>% filter(care_home==0)
sens_nocarehome_outcomes <- outcomes %>% 
  inner_join(sens_nocarehome_df %>% select(patid, indexdate), by=c("patid", "indexdate"))

##----additional adjustment----
#recategorize smoking: current/ex=1, none/missing=0
sens_adj_df <- df %>% mutate(alcohol_cat=case_when(alcohol=="harmful" | alcohol=="higher" ~ 1,
                                                   TRUE ~ 0),
                             smoking_cat=case_when(smoking=="current" | smoking=="ex" ~ 1,
                                                   TRUE ~ 0))


#----Part 2: run analyses----

#create massive functions to run the whole fully adjusted analyses on different datasets for each sensitivity analysis

function_sens_analysis <- function(df, outcomes, sens, covs) {
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
  
  #----create transition matris----
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
  
  
  #----use msprep()----
  
  
  #filter to only relevant rows of indicators and time variables
  outcomes_mstate_formsprep <- outcomes_mstate_times %>% select(patid, indexdate, time_index, time_dis_home, time_dis_carehome, time_read, time_read_home, time_read_carehome, time_death,
                                                                indexhosp, discharge_home, discharge_care_home, readmission, read_discharge_home, read_discharge_care_home, death)
  
  #----use msprep()----
  
  

  #input the time variables and the status (indicators)
  msdata <- msprep(time = c(NA,"time_dis_home","time_dis_carehome","time_read", "time_read_home", "time_read_carehome", "time_death"),
                   status = c(NA,"discharge_home","discharge_care_home","readmission", "read_discharge_home", 
                              "read_discharge_care_home", "death"), 
                   data = outcomes_mstate_formsprep %>%
                     left_join(df %>% select(patid, indexdate, covs, pracid), by=c("patid","indexdate")),
                   trans = tmat,
                   keep=c(covs, "pracid")) 
  
  
  #for those with events occuring on the same day, make the Tstop 0.5 days later
  #for those who transition at time=0, make the Tstop 0.5
  msdata <- msdata %>% mutate(Tstop=Tstop + 
                                if_else(Tstart>0 & Tstart==Tstop, 0.5, 0)) %>%
    mutate(Tstop=Tstop + if_else(status==1 & Tstop==0, 0.5, 0))
  
  class(msdata) <- c("msdata", "data.frame")
  
  #mutating after msprep() loses the trans argument, so need to add it back in
  attr(msdata, "trans") <- tmat
  
  msdata_no_missing_cov <- msdata %>% filter(!is.na(ethnicity_5),
                                             !is.na(imd))
  
  #----transition frequencies----
  
  #save distribution of all events
  events_all_no_missing_cov <- events(msdata_no_missing_cov)
  
  #distribution of events among those with dementia
  events_dementia_no_missing_cov <- events(msdata_no_missing_cov %>% filter(dementia==1))
  

  #distribution of events among those without dementia
  events_no_dementia_no_missing_cov <- events(msdata_no_missing_cov %>% filter(dementia==0))

  n_trans <- length(unique(msdata$trans))
  
  #if there aren't enough events, it saves that transition as NA
  fit_one_trans <- function(k, data, covs) {
    print(k)
    
    fit <- tryCatch({
      sub <- data[data$trans == k, ]
      f <- as.formula(paste("Surv(Tstart, Tstop, status) ~",
                            paste(covs, collapse = " + "),
                            "+ cluster(pracid)"))
      coxph(f, data = sub, method = "breslow",
            x = FALSE, y = FALSE, model = FALSE)
    }, error = function(e) {
      message("Transition ", k, " did not fit: ", conditionMessage(e))
      NULL
    })
    
    rm(fit); gc()
  }
  
  for (k in seq_len(n_trans)) fit_one_trans(k, msdata, covs)

    
  
}

##----running main function----
covs <- c("dementia", "age_c", "sex", "calendar_period", "imd",
          "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
          "ckd", "hypertension")

function_sens_analysis(df=sens_fu_df,
                       outcomes=sens_fu_outcomes,
                       sens="sens_fu",
                       covs=covs)

function_sens_analysis(df=sens_firstinfec_df,
                       outcomes=sens_firstinfec_outcomes,
                       sens="sens_firstinfec",
                       covs=covs)

function_sens_analysis(df=sens_infecfirstep_df,
                       outcomes=sens_infecfirstep_outcomes,
                       sens="sens_infecfirstep",
                       covs=covs)

function_sens_analysis(df=sens_nocarehome_df,
                       outcomes=sens_nocarehome_outcomes,
                       sens="sens_nocarehome",
                       covs=covs)

#Additionally adjusting for bmi, alcohol, smoking, social vulnerability
covs_sens_adj <- c("dementia","age_c","sex","calendar_period","imd","ethnicity_5",
               "efi2","diabetes", "ihd","depression","ckd","hypertension",
               "bmi_classification", "alcohol_cat", "smoking", "social_vulnerability")


function_sens_analysis(df=sens_adj_df,
                       outcomes=outcomes,
                       sens="sens_adj",
                       covs=covs_sens_adj)


##----create forest plot----
#use this function to just change the plot
#reads in the previously-saved numbers of events and cox regression outputs for each sensitivity analysis

create_fp_c4_sens <- function(sens) {
  
  #models_list <- list
  #events_dementia_no_missing_cov <- 
  #events_no_dementia_no_missing_cov <-
  
  # Pull the single row for `term_name` out of each transition's model,
  # in the order given by trans_order, and stack them into one data frame
  plot_df <- map_dfr(1:12, function(k) {
    tidy(models_list[[as.character(k)]], exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == "dementia")
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
                   xlog=TRUE) %>% #added title below
    fp_add_header(n_N_exp = "Dementia\n(n. events/at risk)",
                  n_N_un_exp = "No dementia\n(n. events/at risk)",
                  est_label= "Adjusted HR (95% CI)",
                  term_description = "Transition from index hospitalisation:") |>
    fp_set_style(line = "black",
                 box = "black",
                 zero =gpar(col = "black", lty = "solid")) |>
    fp_decorate_graph(graph.pos = 4)
  fp
  
}
 

create_fp_c4_sens("sens_fu")
create_fp_c4_sens("sens_firstinfec")
create_fp_c4_sens("sens_infecfirstep")
create_fp_c4_sens("sens_nocarehome")
create_fp_c4_sens("sens_adj")



#----Part 3: Clock reset----

#add time in current state as a confounder for to test Markov assumption

#Tstart=time since entering current state
msdata$Tstart_cov <- msdata$Tstart/365.25  

covariates <- c("dementia", "age_c", "sex", "calendar_period", "imd",
                "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                "ckd", "hypertension")

#assessing Markov assumptions:
#fitting model with time since entering current state as a covariate
#running it on the usual residence -> readmission transition (trans==4)
Tstart_cov_model <- coxph(Surv(Tstart, Tstop, status) ~ dementia + Tstart_cov +
                            age_c+ sex+ calendar_period+ imd+
                            ethnicity_5+ efi2+ diabetes+ ihd+ depression+
                            ckd+ hypertension,
                          method = "breslow", 
                          data = msdata %>% filter(trans==4))

#interaction between time since entering state and dementia
Tstart_cov_int_model <- coxph(Surv(Tstart, Tstop, status) ~ dementia*Tstart_cov +
                            age_c+ sex+ calendar_period+ imd+
                          ethnicity_5+ efi2+ diabetes+ ihd+ depression+
                          ckd+ hypertension,
                          method = "breslow", 
                          data = msdata %>% filter(trans==4))

tidy(Tstart_cov_model, exp=T, conf.int = T) 
tidy(Tstart_cov_int_model, exp=T, conf.int = T) 


covariates <- c("dementia", "age_c", "sex", "calendar_period", "imd",
                 "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                 "ckd", "hypertension")

#create a function to fit a cox model for each transition using clock reset

fit_one_trans_reset <- function(k, data, covs) {
  print(k)
  #filter the data by transition
  sub <- data[data$trans == k, ]
  f <- as.formula(paste("Surv(Tstop-Tstart, status) ~",
                        paste(covs, collapse = " + "),
                        "+ cluster(pracid)"))
  fit <- coxph(f, data = sub, method = "breslow",
               x = FALSE, y = FALSE, model = FALSE)
  rm(sub, fit); gc()
}

#run the function for each transition
for (k in seq_len(12)) fit_one_trans_reset(k, msdata, covariates)


#plot the clock reset forest plot
create_fp <- function(models_list) {
  
  
  #events_dementia_no_missing_cov <- 
  #events_no_dementia_no_missing_cov <-
  
  # Pull the single row for `term_name` out of each transition's model,
  # in the order given by trans_order, and stack them into one data frame
  plot_df <- map_dfr(1:12, function(k) {
    tidy(models_list[[as.character(k)]], exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == "dementia")
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
                   xlog=TRUE) %>% #added title below
    fp_add_header(n_N_exp = "Dementia\n(n. events/at risk)",
                  n_N_un_exp = "No dementia\n(n. events/at risk)",
                  est_label= "Adjusted HR (95% CI)",
                  term_description = "Transition from index hospitalisation:") |>
    fp_set_style(line = "black",
                 box = "black",
                 zero =gpar(col = "black", lty = "solid")) |>
    fp_decorate_graph(graph.pos = 4)
  fp
  
}

#models_list <- list
plot_clockreset_sens <- create_fp(models_list)


#----Part 4: Up to 3 readmissions in the model----

#essentially re-doing the pipeline in mstate_analysis script, but adding additional readmissions into the model

##----prepare data----
# indicate when people entered each state and at what date

#among those with no readmission
#create indicators and dates for the index hospitalisation
# Death is derived from ALL spells (index + readmissions) per patid/indexdate,
# using registry death date, then CPRD death date, then spell_end if discharged "Died"
outcome_death <- outcomes %>%
  left_join(df %>% select(patid, indexdate, cprd_ddate), by = c("patid", "indexdate")) %>%
  mutate(death = if_else(!is.na(reg_date_of_death) & (reg_date_of_death - indexdate) <= 365.25 |
                           (is.na(reg_date_of_death) & !is.na(cprd_ddate) & (cprd_ddate - indexdate) <= 365.25) |
                           is.na(reg_date_of_death) & is.na(cprd_ddate) & disdest == "Died", 1, 0),
         date_death = case_when(death == 1 & !is.na(reg_date_of_death) ~ reg_date_of_death,
                                death == 1 & is.na(reg_date_of_death) & !is.na(cprd_ddate) ~ cprd_ddate,
                                death == 1 & is.na(reg_date_of_death) & is.na(cprd_ddate) ~ spell_end,
                                death == 0 ~ NA)) %>%
  group_by(patid, indexdate) %>%
  summarise(death = max(death), date_death = min(date_death, na.rm = TRUE), .groups = "drop") %>%
  mutate(date_death = if_else(is.infinite(date_death), as.Date(NA), date_death))


#index admission
# death is no longer defined here; it is joined from outcome_death below
index_admission <- 
  outcomes %>% 
  #filter to non readmission events (index hosp)
  filter(readmission == 0) %>%
  left_join(df %>% select(patid, indexdate, observation_end, care_home),
            by = c("patid", "indexdate")) %>%
  mutate(indexhosp = 1L,
         discharge_home = if_else(disdest == "Home" | 
                                    (disdest == "Care home" & (care_home == 1 | 
                                                                 admisorc == "Care home")), 1L, 0L),
         date_discharge_home = if_else(disdest == "Home" | 
                                         (disdest == "Care home" & (care_home == 1 | 
                                                                      admisorc == "Care home")), spell_end, NA),
         discharge_care_home = if_else(disdest == "Care home" & 
                                         (is.na(care_home) | care_home == 0) &
                                         admisorc != "Care home", 1L, 0L),
         date_discharge_care_home = if_else(disdest == "Care home" & 
                                              (is.na(care_home) | care_home == 0) &
                                              admisorc != "Care home", spell_end, NA)) %>%
  select(patid, indexdate, indexhosp, observation_end,
         discharge_home, date_discharge_home, discharge_care_home, date_discharge_care_home)


#readmissions 1 to 3
readmission_1_dates <-
  outcomes %>%
  filter(spell_n == 2) %>%
  filter(!is.na(spell_end)) %>% #remove hospitalisations with no spell end
  left_join(df %>% select(patid, indexdate, observation_end, care_home),
            by = c("patid", "indexdate")) %>%
  group_by(patid, indexdate) %>%
  mutate(
    read_1_home = if_else(disdest == "Home" | 
                            (disdest == "Care home" & (care_home == 1 | 
                                                         admisorc == "Care home")), 1L, 0L), 
    date_read_1_home = if_else(disdest == "Home" | 
                                 (disdest == "Care home" & (care_home == 1 | 
                                                              admisorc == "Care home")), spell_end, NA),
    read_1_care_home = if_else(disdest == "Care home" & 
                                 (is.na(care_home) | care_home == 0) &
                                 admisorc != "Care home", 1L, 0L), 
    date_read_1_care_home = if_else(disdest == "Care home" & 
                                      (is.na(care_home) | care_home == 0) &
                                      admisorc != "Care home", spell_end, NA)
  ) %>%
  summarise(
    read_1 = 1L,
    read_1_date = first(spell_start),
    read_1_home = first(read_1_home),
    date_read_1_home = first(date_read_1_home),
    read_1_care_home = first(read_1_care_home),
    date_read_1_care_home = first(date_read_1_care_home),
    .groups = "drop"
  )

readmission_2_dates <-
  outcomes %>%
  filter(spell_n == 3) %>%
  filter(!is.na(spell_end)) %>%
  left_join(df %>% select(patid, indexdate, observation_end, care_home),
            by = c("patid", "indexdate")) %>%
  group_by(patid, indexdate) %>%
  mutate(
    read_2_home = if_else(disdest == "Home" | 
                            (disdest == "Care home" & (care_home == 1 | 
                                                         admisorc == "Care home")), 1L, 0L), 
    date_read_2_home = if_else(disdest == "Home" | 
                                 (disdest == "Care home" & (care_home == 1 | 
                                                              admisorc == "Care home")), spell_end, NA),
    read_2_care_home = if_else(disdest == "Care home" & 
                                 (is.na(care_home) | care_home == 0) &
                                 admisorc != "Care home", 1L, 0L), 
    date_read_2_care_home = if_else(disdest == "Care home" & 
                                      (is.na(care_home) | care_home == 0) &
                                      admisorc != "Care home", spell_end, NA)
  ) %>%
  summarise(
    read_2 = 1L,
    read_2_date = first(spell_start),
    read_2_home = first(read_2_home),
    date_read_2_home = first(date_read_2_home),
    read_2_care_home = first(read_2_care_home),
    date_read_2_care_home = first(date_read_2_care_home),
    .groups = "drop"
  )

readmission_3_dates <-
  outcomes %>%
  filter(spell_n == 4) %>%
  filter(!is.na(spell_end)) %>%
  left_join(df %>% select(patid, indexdate, observation_end, care_home),
            by = c("patid", "indexdate")) %>%
  group_by(patid, indexdate) %>%
  mutate(
    read_3_home = if_else(disdest == "Home" | 
                            (disdest == "Care home" & (care_home == 1 | 
                                                         admisorc == "Care home")), 1L, 0L), 
    date_read_3_home = if_else(disdest == "Home" | 
                                 (disdest == "Care home" & (care_home == 1 | 
                                                              admisorc == "Care home")), spell_end, NA),
    read_3_care_home = if_else(disdest == "Care home" & 
                                 (is.na(care_home) | care_home == 0) &
                                 admisorc != "Care home", 1L, 0L), 
    date_read_3_care_home = if_else(disdest == "Care home" & 
                                      (is.na(care_home) | care_home == 0) &
                                      admisorc != "Care home", spell_end, NA)
  ) %>%
  summarise(
    read_3 = 1L,
    read_3_date = first(spell_start),
    read_3_home = first(read_3_home),
    date_read_3_home = first(date_read_3_home),
    read_3_care_home = first(read_3_care_home),
    date_read_3_care_home = first(date_read_3_care_home),
    .groups = "drop"
  )


#combine
outcomes_mstate <- index_admission %>%
  #add the information from the readmission events
  left_join(readmission_1_dates, by = c("patid", "indexdate")) %>%
  left_join(readmission_2_dates, by = c("patid", "indexdate")) %>%
  left_join(readmission_3_dates, by = c("patid", "indexdate")) %>% 
  #if they did not have a given readmission, replace NA indicators with 0
  mutate(across(c(read_1, read_1_home, read_1_care_home,
                  read_2, read_2_home, read_2_care_home,
                  read_3, read_3_home, read_3_care_home),
                ~ replace_na(.x, 0L))) %>%
  #add death 
  left_join(outcome_death %>% select(patid, indexdate, death, date_death),
            by = c("patid", "indexdate")) %>%
  
  #censoring date: earliest of observation end, death date, or 1 year after index
  mutate(censor_date = pmin(observation_end, date_death, indexdate + 365.25, na.rm = TRUE)) %>%
  
  #if any event happens AFTER censoring date, make the event indicator 0 and the date NA
  mutate(
    discharge_home = if_else(discharge_home == 1 & date_discharge_home > censor_date, 0L, discharge_home),
    date_discharge_home = if_else(discharge_home == 1, date_discharge_home, NA),
    
    discharge_care_home = if_else(discharge_care_home == 1 & date_discharge_care_home > censor_date, 0L, discharge_care_home),
    date_discharge_care_home = if_else(discharge_care_home == 1, date_discharge_care_home, NA),
    
    read_1 = if_else(read_1 == 1 & read_1_date > censor_date, 0L, read_1),
    read_1_date = if_else(read_1 == 1, read_1_date, NA),
    
    read_1_home = if_else(read_1_home == 1 & date_read_1_home > censor_date, 0L, read_1_home),
    date_read_1_home = if_else(read_1_home == 1, date_read_1_home, NA),
    
    read_1_care_home = if_else(read_1_care_home == 1 & date_read_1_care_home > censor_date, 0L, read_1_care_home),
    date_read_1_care_home = if_else(read_1_care_home == 1, date_read_1_care_home, NA),
    
    read_2 = if_else(read_2 == 1 & read_2_date > censor_date, 0L, read_2),
    read_2_date = if_else(read_2 == 1, read_2_date, NA),
    
    read_2_home = if_else(read_2_home == 1 & date_read_2_home > censor_date, 0L, read_2_home),
    date_read_2_home = if_else(read_2_home == 1, date_read_2_home, NA),
    
    read_2_care_home = if_else(read_2_care_home == 1 & date_read_2_care_home > censor_date, 0L, read_2_care_home),
    date_read_2_care_home = if_else(read_2_care_home == 1, date_read_2_care_home, NA),
    
    read_3 = if_else(read_3 == 1 & read_3_date > censor_date, 0L, read_3),
    read_3_date = if_else(read_3 == 1, read_3_date, NA),
    
    read_3_home = if_else(read_3_home == 1 & date_read_3_home > censor_date, 0L, read_3_home),
    date_read_3_home = if_else(read_3_home == 1, date_read_3_home, NA),
    
    read_3_care_home = if_else(read_3_care_home == 1 & date_read_3_care_home > censor_date, 0L, read_3_care_home),
    date_read_3_care_home = if_else(read_3_care_home == 1, date_read_3_care_home, NA)
  ) %>%
  
  #make death take precedence if it happens on the same day as any other transition
  mutate(discharge_home = if_else(!is.na(date_death) & !is.na(date_discharge_home) & date_death == date_discharge_home, 0L, discharge_home),
         date_discharge_home = if_else(!is.na(date_death) & !is.na(date_discharge_home) & date_death == date_discharge_home, NA, date_discharge_home),
         
         discharge_care_home = if_else(!is.na(date_death) & !is.na(date_discharge_care_home) & date_death == date_discharge_care_home, 0L, discharge_care_home),
         date_discharge_care_home = if_else(!is.na(date_death) & !is.na(date_discharge_care_home) & date_death == date_discharge_care_home, NA, date_discharge_care_home),
         
         read_1 = if_else(!is.na(date_death) & !is.na(read_1_date) & date_death == read_1_date, 0L, read_1),
         read_1_date = if_else(!is.na(date_death) & !is.na(read_1_date) & date_death == read_1_date, NA, read_1_date),
         
         read_1_home = if_else(!is.na(date_death) & !is.na(date_read_1_home) & date_death == date_read_1_home, 0L, read_1_home),
         date_read_1_home = if_else(!is.na(date_death) & !is.na(date_read_1_home) & date_death == date_read_1_home, NA, date_read_1_home),
         
         read_1_care_home = if_else(!is.na(date_death) & !is.na(date_read_1_care_home) & date_death == date_read_1_care_home, 0L, read_1_care_home),
         date_read_1_care_home = if_else(!is.na(date_death) & !is.na(date_read_1_care_home) & date_death == date_read_1_care_home, NA, date_read_1_care_home),
         
         read_2 = if_else(!is.na(date_death) & !is.na(read_2_date) & date_death == read_2_date, 0L, read_2),
         read_2_date = if_else(!is.na(date_death) & !is.na(read_2_date) & date_death == read_2_date, NA, read_2_date),
         
         read_2_home = if_else(!is.na(date_death) & !is.na(date_read_2_home) & date_death == date_read_2_home, 0L, read_2_home),
         date_read_2_home = if_else(!is.na(date_death) & !is.na(date_read_2_home) & date_death == date_read_2_home, NA, date_read_2_home),
         
         read_2_care_home = if_else(!is.na(date_death) & !is.na(date_read_2_care_home) & date_death == date_read_2_care_home, 0L, read_2_care_home),
         date_read_2_care_home = if_else(!is.na(date_death) & !is.na(date_read_2_care_home) & date_death == date_read_2_care_home, NA, date_read_2_care_home),
         
         read_3 = if_else(!is.na(date_death) & !is.na(read_3_date) & date_death == read_3_date, 0L, read_3),
         read_3_date = if_else(!is.na(date_death) & !is.na(read_3_date) & date_death == read_3_date, NA, read_3_date),
         
         read_3_home = if_else(!is.na(date_death) & !is.na(date_read_3_home) & date_death == date_read_3_home, 0L, read_3_home),
         date_read_3_home = if_else(!is.na(date_death) & !is.na(date_read_3_home) & date_death == date_read_3_home, NA, date_read_3_home),
         
         read_3_care_home = if_else(!is.na(date_death) & !is.na(date_read_3_care_home) & date_death == date_read_3_care_home, 0L, read_3_care_home),
         date_read_3_care_home = if_else(!is.na(date_death) & !is.na(date_read_3_care_home) & date_death == date_read_3_care_home, NA, date_read_3_care_home))


#add times between states

#use the indicators and dates to create numeric time variables
outcomes_mstate_times <- outcomes_mstate %>%
  #replace any remaining NA indicators with 0 (NEW)
  mutate(across(c(discharge_home, discharge_care_home, 
                  read_1, read_1_home, read_1_care_home,
                  read_2, read_2_home, read_2_care_home,
                  read_3, read_3_home, read_3_care_home, 
                  death),
                ~ replace_na(.x, 0L))) %>%
  mutate(time_index = NA,
         time_dis_home = if_else(discharge_home == 1, 
                                 as.numeric(date_discharge_home - indexdate), 
                                 if_else(discharge_care_home == 0, 
                                         as.numeric(censor_date - indexdate), 
                                         as.numeric(date_discharge_care_home - indexdate))),
         
         time_dis_carehome = if_else(discharge_care_home == 1, 
                                     as.numeric(date_discharge_care_home - indexdate),
                                     if_else(discharge_home == 0, 
                                             as.numeric(censor_date - indexdate), 
                                             as.numeric(date_discharge_home - indexdate))),
         
         time_read_1 = if_else(read_1 == 1, 
                               as.numeric(read_1_date - indexdate), 
                               as.numeric(censor_date - indexdate)),
         
         time_read_1_home = if_else(read_1_home == 1, 
                                    as.numeric(date_read_1_home - indexdate), 
                                    if_else(read_1_care_home == 0, 
                                            as.numeric(censor_date - indexdate), 
                                            as.numeric(date_read_1_care_home - indexdate))),
         
         time_read_1_carehome = if_else(read_1_care_home == 1, 
                                        as.numeric(date_read_1_care_home - indexdate),
                                        if_else(read_1_home == 0, 
                                                as.numeric(censor_date - indexdate), 
                                                as.numeric(date_read_1_home - indexdate))),
         
         time_read_2 = if_else(read_2 == 1, 
                               as.numeric(read_2_date - indexdate), 
                               as.numeric(censor_date - indexdate)),
         
         time_read_2_home = if_else(read_2_home == 1, 
                                    as.numeric(date_read_2_home - indexdate), 
                                    if_else(read_2_care_home == 0, 
                                            as.numeric(censor_date - indexdate), 
                                            as.numeric(date_read_2_care_home - indexdate))),
         
         time_read_2_carehome = if_else(read_2_care_home == 1, 
                                        as.numeric(date_read_2_care_home - indexdate),
                                        if_else(read_2_home == 0, 
                                                as.numeric(censor_date - indexdate), 
                                                as.numeric(date_read_2_home - indexdate))),
         
         time_read_3 = if_else(read_3 == 1, 
                               as.numeric(read_3_date - indexdate), 
                               as.numeric(censor_date - indexdate)),
         
         time_read_3_home = if_else(read_3_home == 1, 
                                    as.numeric(date_read_3_home - indexdate), 
                                    if_else(read_3_care_home == 0, 
                                            as.numeric(censor_date - indexdate), 
                                            as.numeric(date_read_3_care_home - indexdate))),
         
         time_read_3_carehome = if_else(read_3_care_home == 1, 
                                        as.numeric(date_read_3_care_home - indexdate),
                                        if_else(read_3_home == 0, 
                                                as.numeric(censor_date - indexdate), 
                                                as.numeric(date_read_3_home - indexdate))),
         
         #if they died, time_death = death date, if not: censor date
         time_death = if_else(death == 1, 
                              as.numeric(date_death - indexdate), 
                              as.numeric(censor_date - indexdate)))


##----create transition matrix----
#created it manually
tmat <- matrix(NA, 13, 13)

dimnames(tmat) <- list(from = c("indexhosp", 
                                "discharge_home", 
                                "discharge_care_home",
                                "read_1", 
                                "read_1_home", 
                                "read_1_care_home", 
                                "read_2", 
                                "read_2_home", 
                                "read_2_care_home", 
                                "read_3", 
                                "read_3_home", 
                                "read_3_care_home", 
                                "death"),
                       to = c("indexhosp", 
                              "discharge_home", 
                              "discharge_care_home",
                              "read_1", 
                              "read_1_home", 
                              "read_1_care_home", 
                              "read_2", 
                              "read_2_home", 
                              "read_2_care_home", 
                              "read_3", 
                              "read_3_home", 
                              "read_3_care_home", 
                              "death"))

# indexhosp
tmat["indexhosp", "discharge_home"]        <- 1
tmat["indexhosp", "discharge_care_home"]   <- 2
tmat["indexhosp", "death"]                 <- 3

# discharge_home
tmat["discharge_home", "read_1"]           <- 6
tmat["discharge_home", "death"]            <- 4

# discharge_care_home
tmat["discharge_care_home", "read_1"]      <- 7
tmat["discharge_care_home", "death"]       <- 5

# read_1
tmat["read_1", "read_1_home"]              <- 9
tmat["read_1", "read_1_care_home"]         <- 10
tmat["read_1", "death"]                    <- 8

# read_1_home
tmat["read_1_home", "read_2"]              <- 13
tmat["read_1_home", "death"]               <- 11

# read_1_care_home
tmat["read_1_care_home", "read_2"]         <- 14
tmat["read_1_care_home", "death"]          <- 12

# read_2
tmat["read_2", "read_2_home"]              <- 16
tmat["read_2", "read_2_care_home"]         <- 17
tmat["read_2", "death"]                    <- 15

# read_2_home
tmat["read_2_home", "read_3"]              <- 20
tmat["read_2_home", "death"]               <- 18

# read_2_care_home
tmat["read_2_care_home", "read_3"]         <- 21
tmat["read_2_care_home", "death"]          <- 19

# read_3
tmat["read_3", "read_3_home"]              <- 23
tmat["read_3", "read_3_care_home"]         <- 24
tmat["read_3", "death"]                    <- 22

# read_3_home
tmat["read_3_home", "death"]               <- 25

# read_3_care_home
tmat["read_3_care_home", "death"]          <- 26

tmat


###----use msprep()----

#input the time variables and the status (indicators)
msdata <- msprep(time = c(NA,
                          "time_dis_home",
                          "time_dis_carehome",
                          "time_read_1", 
                          "time_read_1_home", 
                          "time_read_1_carehome", 
                          "time_read_2", 
                          "time_read_2_home", 
                          "time_read_2_carehome",
                          "time_read_3", 
                          "time_read_3_home", 
                          "time_read_3_carehome",
                          "time_death"),
                 status = c(NA,
                            "discharge_home",
                            "discharge_care_home",
                            "read_1", 
                            "read_1_home", 
                            "read_1_care_home", 
                            "read_2", 
                            "read_2_home", 
                            "read_2_care_home",
                            "read_3", 
                            "read_3_home", 
                            "read_3_care_home",
                            "death"), 
                 data = outcomes_mstate_times %>%
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






##----HR estimation----


#function fit a cox model for each transition
fit_one_transition <- function(i, msdata, covariates) {
  print(i)
  dat_i <- msdata[msdata$trans == i, ]
  
  if (sum(dat_i$status == 1, na.rm = TRUE) == 0) {
    return(NULL)   # no events -> no model to fit
  }
  
  keep <- covariates[sapply(covariates, function(v) {
    x <- dat_i[[v]]
    if (is.numeric(x)) sd(x, na.rm = TRUE) > 0 else length(unique(x)) > 1
  })]
  
  form_i <- as.formula(paste(
    "Surv(Tstart, Tstop, status) ~",
    paste(keep, collapse = " + ")
  ))
  
  coxph(form_i, data = dat_i, method = "breslow",
        cluster = pracid, x = FALSE, y = FALSE, model = FALSE)
}




# Define all 26 transitions as (from, to) state-name pairs
trans_states <- tibble::tribble(
  ~k,  ~from,                  ~to,                    ~term_description,
  1,   "indexhosp",            "discharge_home",       "Usual residence",
  2,   "indexhosp",            "discharge_care_home",  "New care home",
  3,   "indexhosp",            "death",                "Died in hospital",
  4,   "discharge_home",       "death",                "Usual residence → Died",
  5,   "discharge_care_home",  "death",                "New care home → Died",
  6,   "discharge_home",       "read_1",               "Usual residence → Readmitted (1)",
  7,   "discharge_care_home",  "read_1",               "New care home → Readmitted (1)",
  8,   "read_1",               "death",                "Readmitted (1) → Died in hospital",
  9,   "read_1",               "read_1_home",          "Readmitted (1) → Usual residence",
  10,  "read_1",               "read_1_care_home",     "Readmitted (1) → New care home",
  11,  "read_1_home",          "death",                "Readmitted (1) → Usual residence → Died",
  12,  "read_1_care_home",     "death",                "Readmitted (1) → New care home → Died",
  13,  "read_1_home",          "read_2",               "Readmitted (1) → Usual residence → Readmitted (2)",
  14,  "read_1_care_home",     "read_2",               "Readmitted (1) → New care home → Readmitted (2)",
  15,  "read_2",               "death",                "Readmitted (2) → Died in hospital",
  16,  "read_2",               "read_2_home",          "Readmitted (2) → Usual residence",
  17,  "read_2",               "read_2_care_home",     "Readmitted (2) → New care home",
  18,  "read_2_home",          "death",                "Readmitted (2) → Usual residence → Died",
  19,  "read_2_care_home",     "death",                "Readmitted (2) → New care home → Died",
  20,  "read_2_home",          "read_3",               "Readmitted (2) → Usual residence → Readmitted (3)",
  21,  "read_2_care_home",     "read_3",               "Readmitted (2) → New care home → Readmitted (3)",
  22,  "read_3",               "death",                "Readmitted (3) → Died in hospital",
  23,  "read_3",               "read_3_home",          "Readmitted (3) → Usual residence",
  24,  "read_3",               "read_3_care_home",     "Readmitted (3) → New care home",
  25,  "read_3_home",          "death",                "Readmitted (3) → Usual residence → Died",
  26,  "read_3_care_home",     "death",                "Readmitted (3) → New care home → Died"
)

#function to create forest plot of HRs
create_fp_extra_states <- function(models_list, trans_order, term_name) {
  
  # --- 1. Pull HR/CI per transition, defensively (handles missing models
  #        and zero-event transitions where `term_name` was dropped) ---
  plot_df <- purrr::map_dfr(1:26, function(k) {
    mod <- models_list[[(k)]]
    
    if (is.null(mod) || inherits(mod, "try-error")) {
      return(tibble::tibble(term = "dementia", estimate = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_))
    }
    
    tt <- tryCatch(
      broom::tidy(mod, exponentiate = TRUE, conf.int = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(tt) || !"term" %in% names(tt) || !"dementia" %in% tt$term) {
      return(tibble::tibble(term = "dementia", estimate = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_))
    }
    
    tt |>
      dplyr::filter(term == "dementia") |>
      dplyr::select(term, estimate, conf.low, conf.high)
  })
  
  # --- 2. Pull events/N by state name, not position ---
  get_events <- function(freq_table) {
    trans_states |>
      dplyr::rowwise() |>
      dplyr::mutate(
        events = as.numeric(freq_table[from, to]),
        N      = as.numeric(freq_table[from, "total entering"])
      ) |>
      dplyr::ungroup() |>
      dplyr::select(k, events, N)
  }
  
  events_exp_df   <- get_events(events_dementia_no_missing_cov$Frequencies)
  events_unexp_df <- get_events(events_no_dementia_no_missing_cov$Frequencies)
  
  # --- 3. Assemble plot_df: bind term_description + both event sets ---
  plot_df <- plot_df |>
    dplyr::mutate(k = 1:26) |>
    dplyr::left_join(trans_states |> dplyr::select(k, term_description), by = "k") |>
    dplyr::left_join(events_exp_df   |> dplyr::rename(events_exp = events, sample_size_exp = N), by = "k") |>
    dplyr::left_join(events_unexp_df |> dplyr::rename(events_un_exp = events, sample_size_un_exp = N), by = "k") |>
    dplyr::mutate(
      term_description = factor(term_description, levels = trans_states$term_description)
    ) |>
    dplyr::mutate(n_N_exp = paste0(dplyr::if_else(events_exp %in% 1:50, "<50", format(events_exp, big.mark = ",")),
                                   "/", format(sample_size_exp, big.mark = ",")),
                  n_N_un_exp = paste0(dplyr::if_else(events_un_exp %in% 1:50, "<50", format(events_un_exp, big.mark = ",")),
                                      "/", format(sample_size_un_exp, big.mark = ","))) |>
    dplyr::mutate(dplyr::across(c(conf.low, conf.high), ~ dplyr::if_else(.x < 0.0001 | .x > 100, 0, .x))) |>
    dplyr::mutate(est_label = dplyr::if_else(is.na(estimate), NA_character_,
                                             paste0(formatC(estimate, format = "f", digits = 2),
                                                    " (", formatC(conf.low, format = "f", digits = 2), ", ",
                                                    formatC(conf.high, format = "f", digits = 2), ")")),
                  dplyr::across(c(estimate, conf.low, conf.high), ~ round(.x, 2))) |>
    dplyr::mutate(dplyr::across(c(estimate, conf.low, conf.high), ~ dplyr::if_else(.x < 0.0001 | .x > 100, NA, .x))) |>
    dplyr::mutate(est_label_middot = stringr::str_replace_all(est_label, "[.]", "\U00B7"),
                  estimate_middot  = stringr::str_replace(estimate, "[.]", "\U00B7")) |>
    # Zero events in either group -> show N but blank out the HR
    dplyr::mutate(dplyr::across(c(conf.low, conf.high, est_label, estimate),
                                ~ dplyr::if_else(events_exp == 0 | events_un_exp == 0, NA, .x)))
  
  fp <- forestplot::forestplot(plot_df,
                               labeltext = c(term_description, n_N_exp, n_N_un_exp, est_label),
                               mean = estimate,
                               lower = conf.low,
                               upper = conf.high,
                               zero = 1,
                               xlab = "Hazard Ratio (95% CI)",
                               xticks = c(0.7, 1, 1.4, 2),
                               txt_gp = forestplot::fpTxtGp(
                                 ticks = grid::gpar(cex = 1),
                                 xlab  = grid::gpar(cex = 1),
                                 legend = grid::gpar(cex = 0.85)
                               ),
                               boxsize = 0.15,
                               xlog = TRUE) |>
    forestplot::fp_add_header(n_N_exp = "Dementia\n(n. events/at risk)",
                              n_N_un_exp = "No dementia\n(n. events/at risk)",
                              est_label = "Adjusted HR (95% CI)",
                              term_description = "Transition from index hospitalisation:") |>
    forestplot::fp_set_style(line = "black",
                             box = "black",
                             zero = grid::gpar(col = "black", lty = "solid")) |>
    forestplot::fp_decorate_graph(graph.pos = 4)
  
  print(fp)
}


#fully adjusted model covariates
covariates <- c("dementia", "age_c", "sex", "calendar_period", "imd",
                "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                "ckd", "hypertension")
#run each transition-specific HR and save as a list
models_list <- lapply(seq_len(n_trans), fit_one_transition,
                      msdata = msdata, covariates = covariates)

#run the function to plot the HR including extra readmissions
extra_states_hr_plot <- create_fp_extra_states(models_list, trans_order = 1:26, term_name = "dementia")


##----prediction----


#model is too big to run in one nodel, so fitting separate models for each transition then manually combining them for msprep()


covariates_for_pred <- c("dementia", "age_c", "sex", "calendar_period",
                         "ethnicity_5", "imd", "efi2")

n_trans <- max(msdata$trans)

fits      <- vector("list", n_trans)
basehaz_l <- vector("list", n_trans)

fits <- lapply(seq_len(n_trans), fit_one_transition,
               msdata = msdata, covariates = covariates_for_pred)


basehaz_l <- lapply(seq_along(fits), function(i) {
  if (is.null(fits[[i]])) {
    # zero-event transition: hazard is 0 at a single nominal time point
    return(data.frame(time = 0, Haz = 0, trans = i))
  }
  bh <- basehaz(fits[[i]], centered = FALSE)
  data.frame(time = bh$time, Haz = bh$hazard, trans = i)
})


# 1. Define the two covariate profiles — no trans/strata columns, no expand.covs
newdata_un_exp <- data.frame(
  dementia        = 0,
  age_c           = 0,
  sex             = factor("Female", levels = levels(df$sex)),
  calendar_period = factor("2018-2019", levels = levels(df$calendar_period)),
  ethnicity_5     = factor("White", levels = levels(as.factor(df$ethnicity_5))),
  imd             = factor("5", levels = levels(df$imd)),
  efi2            = median(df$efi2)
)

newdata_exp <- newdata_un_exp
newdata_exp$dementia <- 1   # exposed group

newdata_un_exp$pracid <- msdata$pracid[1]   # any valid existing pracid value
newdata_exp$pracid    <- msdata$pracid[1]

# 2. Baseline Haz data frame (built once, reused for both profiles)
Haz_base <- do.call(rbind, basehaz_l)

get_lp <- function(fit, newdata) {
  tt <- delete.response(terms(fit))            # drop the Surv() response part
  mf <- model.frame(tt, data = newdata, xlev = fit$xlevels)
  mm <- model.matrix(tt, mf, contrasts.arg = fit$contrasts)
  mm <- mm[, names(coef(fit)), drop = FALSE]    # align to exactly the fitted coefficients
  as.numeric(mm %*% coef(fit))
}

build_msfit <- function(newdata_profile, fits, Haz_base, tmat) {
  Haz <- Haz_base
  for (i in seq_along(fits)) {
    if (is.null(fits[[i]])) next   # zero-hazard transition, leave as 0
    lp_i <- get_lp(fits[[i]], newdata_profile)
    Haz$Haz[Haz$trans == i] <- Haz$Haz[Haz$trans == i] * exp(lp_i)
  }
  msf <- list(Haz = Haz, trans = tmat)
  class(msf) <- "msfit"
  msf
}

msf_un_exp <- build_msfit(newdata_un_exp, fits, Haz_base, tmat)
msf_exp    <- build_msfit(newdata_exp,    fits, Haz_base, tmat)



#now calculate the probability of each transition
#predt= time to start predictions (time 0)
#probtrans() uses the baseline hazards (from msfit) to predict the patient-specific transition probabilities

probtrans_un_exp <- probtrans(msf_un_exp,predt = 0, variance=FALSE)
probtrans_exp <- probtrans(msf_exp,predt = 0, variance=FALSE)


###----stacked occupancy plots----

modified_un_exp <- probtrans_un_exp[[1]] %>%
  mutate(hospital=pstate1+pstate4+pstate7+pstate10,
         home=pstate2+pstate5+pstate8+pstate11,
         care_home=pstate3+pstate6+pstate9+pstate12,
         death=pstate13) %>%
  select(time,hospital,home,care_home,death) %>%
  pivot_longer(cols = c(hospital, home, care_home, death),
               names_to = "state",
               values_to = "prob") %>%
  mutate(state = factor(state,
                        levels = c("death", "care_home", "home", "hospital"),
                        labels = c("Death", "Care home", "Usual residence", "Hospital")))

modified_exp <- probtrans_exp[[1]] %>%
  mutate(hospital=pstate1+pstate4+pstate7+pstate10,
         home=pstate2+pstate5+pstate8+pstate11,
         care_home=pstate3+pstate6+pstate9+pstate12,
         death=pstate13) %>%
  select(time,hospital,home,care_home,death) %>%
  pivot_longer(cols = c(hospital, home, care_home, death),
               names_to = "state",
               values_to = "prob") %>%
  mutate(state = factor(state,
                        levels = c("death", "care_home", "home", "hospital"),
                        labels = c("Death", "Care home", "Usual residence", "Hospital")))

plot_collapsed_categories <- function(data, title) {
  
  ggplot(data, aes(x = time, y = prob, fill = state)) +
    geom_area(linewidth = 0.75, colour = "white") +
    scale_fill_viridis_d(direction = -1,
                         breaks = c("Hospital", "Usual residence", "Care home", "Death")) +
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


dementia_sensplot <- plot_collapsed_categories(modified_exp, "A. Dementia")
nodementia_sensplot <- plot_collapsed_categories(modified_un_exp, "B. No dementia")


####----6 month occupancy----

state_occu_prob_6_mo <- rbind(
  slice(filter(probtrans_exp[[1]],time>=183),1),
  slice(filter(probtrans_un_exp[[1]],time>=183),1))|>
  mutate(Hospital = pstate1 + pstate4 + pstate7+pstate10,
         Usual_res = pstate2+ pstate5 + pstate8 + pstate11,
         Care_home = pstate3 + pstate6 + pstate9+pstate12,
         Death = pstate13) |>
  select(time, Hospital, Usual_res, Care_home, Death) %>%
  mutate(dementia=c(1,0))




