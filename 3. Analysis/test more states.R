


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


#setwd("J:/EHR-Working/Monica I/Care pathways/Care pathways code/R/Analysis/")

#----read in data----

#df has patid, indexdate, and all baseline covariates
#patients can enter the cohort multiple times

df <- read_parquet("J:/EHR-Working/Monica I/Care pathways/Parquet files/study_pop_with_variables_6.parquet")

#additional data cleaning, ensuring there is no duplicate rows
df <- df %>% distinct() %>%
  #create follow-up variable
  mutate(follow_up=ifelse(observation_end-indexdate<365, observation_end-indexdate, 365)) %>% 
  #rename ethnicity for ease later one
  mutate(ethnicity_5=case_when(ethnicity_5=="0. White" ~ "White",
                               ethnicity_5=="1. South Asian" ~ "South Asian",
                               ethnicity_5=="2. Black" ~ "Black",
                               ethnicity_5=="3. Other" ~ "Other",
                               ethnicity_5=="4. Mixed" ~ "Mixed",
                               ethnicity_5=="5. Not Stated" ~ NA)) %>%
  rename(efi2=efi2_score_excl_dementia) %>%
  #create quintiles of efi2
  mutate(efi2_quintile=ntile(efi2, 5)) %>%
  mutate(imd=as.factor(imd)) %>%
  #remove people who have missing sex
  subset(sex!="Indeterminate") %>%
  mutate(sex=as.factor(sex)) %>%
  mutate(calendar_period=as.factor(calendar_period)) %>%
  #creating a centered age variable
  mutate(age_c = age - median(age, na.rm=TRUE)) %>%
  #create age group variable
  mutate(age_group=case_when(age <70 ~ "65-69",
                             age >=70 & age <75 ~ "70-74",
                             age >=75 & age <80 ~ "75-79",
                             age >=80 & age <85 ~ "80-84",
                             age >=85 & age <90 ~ "85-89",
                             age >=90 ~ "90+"))

#outcomes has rows for each hospital spell, and indicators for the outcomes (readmission, failed discharge)
#multiple rows per patid/indexdate 
outcomes <- read_parquet("J:/EHR-Working/Monica I/Care pathways/Parquet files/outcomes_6.parquet")
outcomes <- outcomes %>% distinct()

#create a variable to indicate that people lived in a care home if they had GP code OR admisorc was care home
df <- df %>% left_join(outcomes %>% filter(spell_n==1), by=c("patid", "indexdate")) %>%
  mutate(care_home_incl_admisorc=if_else(care_home==1 | admisorc=="Care home", 1, 0))

#create a variable to lump in living alone and care home as a 3 tiered variable
df <- df %>%
  mutate(living_status=case_when(living_alone==1 ~ "Living alone",
                                 care_home_incl_admisorc==1 ~ "Care home",
                                 TRUE ~ "Community-dwelling with others"))



#----format data for msprep----

##----add transition indicators----
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


##---- index admission ----
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


##---- readmissions 1 to 3 ----
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


##---- combine ----
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


##----add times between states----

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


#----create transition matris----
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


#----use msprep()----

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


covariates3 <- c("dementia", "age_c", "sex", "calendar_period", "imd", 
                 "ethnicity_5", "efi2")

#the below code tries to fit the whole model in one but it is too big and doesn't have enough memory
#msdata_expcov <- expand.covs(msdata, covariates3, append=TRUE, longnames = FALSE)

#expanded_names_c3 <- grep(
#  paste0("^(", paste(covariates3, collapse = "|"), ")[A-Za-z0-9]*\\.[0-9]+$"),
#  names(msdata_expcov), value = TRUE)

#c3 <- coxph(as.formula(paste("Surv(Tstart, Tstop, status) ~",
#                             paste(expanded_names_c3, collapse = " + "),
#                             "+ strata(as.factor(trans))+cluster(pracid)")), 
#            method = "breslow", 
#            data = msdata_expcov)

#saveRDS(msdata_expcov, "./test_more_states/msdata_expcov.rds")
#msdata_expcov <- readRDS("./test_more_states/msdata_expcov.rds")

#---trying to manually get the baseline hazard for each tranision separately
#with help of AI

n_trans <- max(msdata$trans)

fits      <- vector("list", n_trans)
basehaz_l <- vector("list", n_trans)

# Modify the fitting loop/lapply to skip zero-event transitions
fit_one_transition <- function(i, msdata, covariates3) {
  print(i)
  dat_i <- msdata[msdata$trans == i, ]
  
  if (sum(dat_i$status == 1, na.rm = TRUE) == 0) {
    return(NULL)   # no events -> no model to fit
  }
  
  keep <- covariates3[sapply(covariates3, function(v) {
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

fits <- lapply(seq_len(n_trans), fit_one_transition,
               msdata = msdata, covariates3 = covariates3)


basehaz_l <- lapply(seq_along(fits), function(i) {
  if (is.null(fits[[i]])) {
    # zero-event transition: hazard is 0 at a single nominal time point
    return(data.frame(time = 0, Haz = 0, trans = i))
  }
  bh <- basehaz(fits[[i]], centered = FALSE)
  data.frame(time = bh$time, Haz = bh$hazard, trans = i)
})

#----prediction----
#first need to use estimated parameters and baseline transition hazards and covaraite values to obtain patient-specific transition hazards
#create a new dataset contaninig expanded covariates, one line per transition


covariates_for_pred <- c("dementia", "age_c", "sex", "calendar_period",
                         "ethnicity_5", "imd", "efi2")
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
head(probtrans_exp[[1]])
head(probtrans_un_exp[[1]])


#----stacked occupancy plots----


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

saveRDS(dementia_sensplot, "./test_more_states/dementia_sensplot.rds")
saveRDS(nodementia_sensplot, "./test_more_states/nodementia_sensplot.rds")

dementia_sensplot <- readRDS("./test_more_states/dementia_sensplot.rds")
nodementia_sensplot <- readRDS("./test_more_states/nodementia_sensplot.rds")

#----expected length of stay in each state----

#ELOS() function uses the probtrans() object to estimate number of days in each state
#using the diagonal values gives the estimate CONDITIONAL on having reached that state

elos_exp <- round(ELOS(probtrans_exp),1)
elos_un_exp <- round(ELOS(probtrans_un_exp),1)

#----6 month occupancy----

state_occu_prob_6_mo <- rbind(
  slice(filter(probtrans_exp[[1]],time>=183),1),
  slice(filter(probtrans_un_exp[[1]],time>=183),1))|>
  mutate(Hospital = pstate1 + pstate4 + pstate7+pstate10,
         Usual_res = pstate2+ pstate5 + pstate8 + pstate11,
         Care_home = pstate3 + pstate6 + pstate9+pstate12,
         Death = pstate13) |>
  select(time, Hospital, Usual_res, Care_home, Death) %>%
  mutate(dementia=c(1,0))


#----looking at HRs----

msdata <- readRDS("./test_more_states/msdata.rds")
msdata_no_missing_cov <- msdata %>% filter(!is.na(ethnicity_5),
                                           !is.na(imd))

events_all <- events(msdata)
events_all_no_missing_cov <- events(msdata_no_missing_cov)

events_dementia <- events(msdata %>% filter(dementia==1))
events_dementia_no_missing_cov <- events(msdata_no_missing_cov %>% filter(dementia==1))

events_no_dementia <- events(msdata %>% filter(dementia==0))
events_no_dementia_no_missing_cov <- events(msdata_no_missing_cov %>% filter(dementia==0))

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

covariates <- c("dementia", "age_c", "sex", "calendar_period", "imd",
                "ethnicity_5", "efi2", "diabetes", "ihd", "depression",
                "ckd", "hypertension")
models_list <- lapply(seq_len(n_trans), fit_one_transition,
               msdata = msdata, covariates = covariates)



extra_states_hr_plot <- create_fp_extra_states(models_list, trans_order = 1:26, term_name = "dementia")

saveRDS(extra_states_hr_plot, "./test_more_states/extra_states_hr_plot.rds")
