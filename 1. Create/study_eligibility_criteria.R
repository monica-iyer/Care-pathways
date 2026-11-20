

#Monica Iyer
#Care pathways CPRD study
#April 9, 2026
#Study eligibility criteria



library(tidyverse)
library(stringr)
library(haven)
library(data.table)
library(arrow)
library(dplyr)
library(duckdb)

###read in the data######


#reading in hes hospital data
#this data is from a Type 2 linkage request with CPRD and contains HES data for the study population
#placeholders where the data should be

#hospital-level data
raw_hes_hosp_parquet <- open_dataset("")

#episode level data with primary diagnosis only

raw_hes_primary_diagnosis_epi_parquet <- open_dataset("")

#infection code list
infection_codelist <- readLines("")

#read in denominator file
#September 2025 data cut
denom_0925_parquet <- open_dataset("")

#read in eligibility linkage file
elig_parquet <- open_dataset("")

#read in linked ONS data
ons_parquet <- open_dataset("")

###applying eligibility criteria####

#step 1: apply first phase of eligibility criteria


elig_65 <- denom_0925_parquet %>%
  to_duckdb() %>%
  inner_join(to_duckdb(elig_parquet), by="patid") %>%
  #require eligibility to be linked to ONS and HES APC
  filter(hes_apc_e==1 & ons_death_e==1) %>%
  inner_join(to_duckdb(raw_hes_primary_diagnosis_epi_parquet), by="patid") %>%
  #only include events with infections in the codelist
  mutate(infection_dx=as.integer(ICD %in% infection_codelist)) %>%
  filter(infection_dx==1) %>%
  inner_join(to_duckdb(raw_hes_hosp_parquet) %>% select(-patid), by="spno") %>%
  mutate(admidate   = as.Date(strptime(admidate, "%Y/%m/%d")),
         discharged = as.Date(strptime(discharged, "%Y/%m/%d")),
         epistart   = as.Date(strptime(epistart, "%Y/%m/%d")),
         epiend     = as.Date(strptime(epiend, "%Y/%m/%d"))) %>%
  mutate(time_between_admit_epi= epistart-admidate) %>%
  #infection dx within a week of hospital admission
  filter(time_between_admit_epi <8 & time_between_admit_epi >=0) %>% 
  #infection dx on or after Jan 1, 2010 (study start)
  filter(admidate >= as.Date("2010-01-01")) %>%
  #clean the data - no missing discharge date; discharge date has to be on or after admitdate
  filter(!is.na(discharged) & 
           discharged>=admidate) %>%
  collect() %>%
  #estimate date of birth as June 1, year of birth
  mutate(yob=as.Date(paste0(yob, "-06-01"))) %>%
  #at least 65 years old
  mutate(age_at_event = as.numeric(difftime(admidate, yob, units = "days")/365.25)) %>%
  filter(age_at_event >= 65) %>%
  as_arrow_table()

#Step 2: add ONS data

elig_65_dates <- elig_65 %>%
  to_duckdb() %>%                                              
  mutate(regstartdate = as.Date(strptime(regstartdate, "%d/%m/%Y")),
         regenddate   = as.Date(strptime(regenddate,   "%d/%m/%Y")),
         cprd_ddate   = as.Date(strptime(cprd_ddate,   "%d/%m/%Y")),
         lcd          = as.Date(strptime(lcd,          "%d/%m/%Y"))) %>%
  select(-any_of("reg_date_of_death")) %>%                     
  left_join(to_duckdb(ons_parquet), by = "patid") %>%
  mutate(reg_date_of_death = as.Date(strptime(reg_date_of_death, "%Y-%m-%d"))) %>%
  collect() %>%                                                
  as_arrow_table()                                             

#Step 3: require 6 months minimum registration; set observation_end

study_pop_obs <- elig_65_dates %>%
  collect() %>%                                  
  # use ONS as death date; if missing, use CPRD death date
  mutate(death_date = coalesce(reg_date_of_death, cprd_ddate)) %>%
  # study end = March 31, 2025 - end of HES linkage 
  mutate(observation_end = pmin(regenddate, death_date, lcd,
                                as.Date("2025-03-31"), na.rm = TRUE)) %>%
  filter(as.numeric(admidate - regstartdate) >= 182.5,
         observation_end >= admidate)

#Step 4: use the prior dataset, which has every eligible observation, to keep only observations a year apart
#participants can re-enter the cohort if they have another eligible event more than a year later


study_pop <- study_pop_obs %>%
  arrange(patid, admidate) %>% #arrange by admitdate
  group_by(patid) %>%
  filter({
    keep <- logical(n()) #creates an empty true/false vector for every event per patid
    last_kept_date <- as.Date(NA)
    
    for (i in seq_len(n())) { #for every event per patient, keep the first event
      if (i==1) {
        keep[i] <- TRUE
        last_kept_date <- admidate[i]
      } else { #if it's not their first event but the event is at least a year after the last event then keep it
        if(as.numeric(difftime(admidate[i], last_kept_date, units="days")) >= 365.25) {
          keep[i] <- TRUE
          last_kept_date <- admidate[i]
        }
      }
    }
    keep
  }) %>%
  ungroup() %>%
  
  select(-pracid.y) %>%
  rename(pracid=pracid.x) %>%
  rename(indexdate=admidate)



study_pop <- study_pop %>% select(patid,
                       indexdate,
                       pracid,
                       gender,
                       yob,
                       regstartdate,
                       regenddate,
                       lcd,
                       region,
                       spno,
                       epikey,
                       epistart,
                       epiend,
                       ICD,
                       d_order,
                       infection_dx,
                       discharged,
                       admisorc,
                       disdest,
                       age_at_event,
                       cprd_ddate,
                       reg_date_of_death,
                       death_date,
                       observation_end)

#study_pop has the eligible events for the study population
#this will now be used in the targets script to create the study variables

###sensitivity analysis####
#one of the sensitivity analyses requires including only individuals who have infection as the primary diagnosis on the first episode of admission
#repeating the above steps, with one change



#step 1: apply first phase of eligibility criteria


elig_65 <- denom_0925_parquet %>%
  to_duckdb() %>%
  inner_join(to_duckdb(elig_parquet), by="patid") %>%
  filter(hes_apc_e==1 & ons_death_e==1) %>%
  inner_join(to_duckdb(raw_hes_primary_diagnosis_epi_parquet), by="patid") %>%
  mutate(infection_dx=as.integer(ICD %in% infection_codelist)) %>%
  filter(infection_dx==1) %>%
  inner_join(to_duckdb(raw_hes_hosp_parquet) %>% select(-patid), by="spno") %>%
  mutate(admidate   = as.Date(strptime(admidate, "%Y/%m/%d")),
         discharged = as.Date(strptime(discharged, "%Y/%m/%d")),
         epistart   = as.Date(strptime(epistart, "%Y/%m/%d")),
         epiend     = as.Date(strptime(epiend, "%Y/%m/%d"))) %>%
  filter(admidate==epistart) %>% #this is where the change is, require episode start=admidate
  filter(admidate >= as.Date("2010-01-01")) %>%
  filter(!is.na(discharged) & 
           discharged>=admidate) %>%
  collect() %>%
  mutate(yob=as.Date(paste0(yob, "-06-01"))) %>%
  mutate(age_at_event = as.numeric(difftime(admidate, yob, units = "days")/365.25)) %>%
  filter(age_at_event >= 65) %>%
  as_arrow_table()

elig_65_dates <- elig_65 %>%
  to_duckdb() %>%                                              
  mutate(regstartdate = as.Date(strptime(regstartdate, "%d/%m/%Y")),
         regenddate   = as.Date(strptime(regenddate,   "%d/%m/%Y")),
         cprd_ddate   = as.Date(strptime(cprd_ddate,   "%d/%m/%Y")),
         lcd          = as.Date(strptime(lcd,          "%d/%m/%Y"))) %>%
  select(-any_of("reg_date_of_death")) %>%                     
  left_join(to_duckdb(ons_parquet), by = "patid") %>%
  mutate(reg_date_of_death = as.Date(strptime(reg_date_of_death, "%Y-%m-%d"))) %>%
  collect() %>%                                                
  as_arrow_table()                                             

study_pop_obs <- elig_65_dates %>%
  collect() %>%                                  
  mutate(death_date = coalesce(reg_date_of_death, cprd_ddate)) %>%
  mutate(observation_end = pmin(regenddate, death_date, lcd,
                                as.Date("2025-03-31"), na.rm = TRUE)) %>%
  filter(as.numeric(admidate - regstartdate) >= 182.5,
         observation_end >= admidate)

sens_study_pop <- study_pop_obs %>%
  arrange(patid, admidate) %>% 
  group_by(patid) %>%
  filter({
    keep <- logical(n()) 
    last_kept_date <- as.Date(NA)
    
    for (i in seq_len(n())) { 
      if (i==1) {
        keep[i] <- TRUE
        last_kept_date <- admidate[i]
      } else { 
        if(as.numeric(difftime(admidate[i], last_kept_date, units="days")) >= 365.25) {
          keep[i] <- TRUE
          last_kept_date <- admidate[i]
        }
      }
    }
    keep
  }) %>%
  ungroup() %>%
  
  select(-pracid.y) %>%
  rename(pracid=pracid.x) %>%
  rename(indexdate=admidate)

sens_study_pop <- sens_study_pop %>% select(patid,
                                  indexdate,
                                  pracid,
                                  gender,
                                  yob,
                                  regstartdate,
                                  regenddate,
                                  lcd,
                                  region,
                                  spno,
                                  epikey,
                                  epistart,
                                  epiend,
                                  ICD,
                                  d_order,
                                  infection_dx,
                                  discharged,
                                  admisorc,
                                  disdest,
                                  age_at_event,
                                  cprd_ddate,
                                  reg_date_of_death,
                                  death_date,
                                  observation_end)





