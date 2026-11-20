
#--Monica Iyer
#--June 2026
#--Analysis for care pathways following infections study

#This script contains code to obtain descriptive statistics, both of baseline information and outcomes
#Script corresponds to Table 1 (baseline characteristics), Table S2 (baseline characteristics by missing ethnicity), Table S3 (outcomes)
#note that the Table S3 outcomes of discharge destination after readmission are taken from events() in the mstate_analysis script


library(dplyr)
library(arrow)
library(table1)
library(tidyr)
library(mstate)
library(survival)
library(readr)
library(tidyverse);library(janitor);library(magrittr);
library(gtsummary);library(labelled);
library(purrr);library(glue);library(gt);library(readxl);
library(geepack);library(evaluate);library(meta);library(forestplot);
library(cowplot);library(rio);library(Hmisc);
library(rsconnect);library(forestploter);library(ggpubr);library(msm);
library(RColorBrewer);library(grid)



#----read in parquet files for the baseline variables and outcomes----

#read in study_pop_with_variables
df <- read_parquet("./study_pop_with_variables.parquet")

#clean dataset
df <- df %>% distinct() %>%
  mutate(follow_up=ifelse(observation_end-indexdate<365, observation_end-indexdate, 365)) %>% 
  mutate(ethnicity_5=case_when(ethnicity_5=="0. White" ~ "White",
                               ethnicity_5=="1. South Asian" ~ "South Asian",
                               ethnicity_5=="2. Black" ~ "Black",
                               ethnicity_5=="3. Other" ~ "Other",
                               ethnicity_5=="4. Mixed" ~ "Mixed",
                               ethnicity_5=="5. Not Stated" ~ NA)) %>%
  rename(efi2=efi2_score_excl_dementia) %>%
  mutate(efi2_quintile=ntile(efi2, 5)) %>%
  mutate(imd=as.factor(imd)) %>%
  mutate(sex=as.factor(sex)) %>%
  mutate(calendar_period=as.factor(calendar_period)) %>%
  #creating a centered age variable
  mutate(age_c = age - mean(age, na.rm=TRUE)) %>%
  mutate(age_group=case_when(age <70 ~ "65-69",
                      age >=70 & age <75 ~ "70-74",
                      age >=75 & age <80 ~ "75-79",
                      age >=80 & age <85 ~ "80-84",
                      age >=85 & age <90 ~ "85-89",
                      age >=90 ~ "90+")) %>%
  mutate(follow_up=ifelse(observation_end-indexdate<365, observation_end-indexdate, 365))




#read in the outcomes

outcomes <- read_parquet("./outcomes.parquet")
outcomes <- outcomes %>% distinct()


#create variable (new_care_home) to indicate that people lived in a care home if they had GP code OR admisorc was care home
new_care_home <- outcomes %>% 
  filter(spell_n==1) %>%
  inner_join(df %>% select(patid, indexdate, care_home, living_alone), by=c("patid", "indexdate")) %>%
  mutate(new_care_home=if_else(disdest=="Care home" & 
                                 (is.na(care_home) | care_home==0) &
                                 admisorc!="Care home", 1, 0)) %>%
  mutate(admisorc_carehome=if_else(admisorc=="Care home", 1, 0))


#create a variable to indicate that people lived in a prior care home if they had GP code OR admisorc was care home
df <- df %>% left_join(outcomes %>% filter(spell_n==1), by=c("patid", "indexdate")) %>%
  mutate(care_home_incl_admisorc=if_else(care_home==1 | admisorc=="Care home", 1, 0))

#edit the living alone/care home to be one variable
df <- df %>%
  mutate(living_status=case_when(living_alone==1 ~ "Living alone",
                                 care_home_incl_admisorc==1 ~ "Care home",
                                 TRUE ~ "Community-dwelling with others"))


#----descriptives----

##----table 1----

#create table 1, characteristics by dementia status
table1(~age + sex + calendar_period + ethnicity_5+ as.factor(imd) + as.factor(ihd)+ as.factor(diabetes)+
         as.factor(ckd)+as.factor(hypertension)+as.factor(depression) + 
         efi2  +as.factor(efi2_quintile)+as.factor(gp_consultation_prior_year) +
         (delirium_index) + as.factor(living_alone) +
         as.factor(care_home)+as.factor(care_home_incl_admisorc) +
         as.factor(living_status)| dementia, df)


#getting IQR of age
summary(df$age)
summary((df%>% filter(dementia==1))$age)
summary((df%>% filter(dementia==0))$age)


#number of unique participants in dementia/no dementia groups
df %>% filter(dementia==1) %>% distinct(patid) %>% nrow()
df %>% filter(dementia==0) %>% distinct(patid) %>% nrow()



##----exploring missing ethnicity; Table S2----#

missing_eth_test <- df %>% mutate(eth_missing=case_when(ethnicity_5=="5. Not Stated" | 
                                                          is.na(ethnicity_5) ~ "Ethnicity missing",
                                                        TRUE ~ "Ethnicity recorded"))

#checking the distribution of missing ethnicity by all varibles
table1(~age + sex + calendar_period + as.factor(dementia)+ as.factor(imd) + as.factor(ihd)+ as.factor(diabetes)+
         as.factor(ckd)+as.factor(hypertension)+as.factor(depression) + efi2+
         as.factor(efi2_quintile) +as.factor(living_status)| eth_missing, missing_eth_test)

#IQR of age
missing_eth_test %>% filter(eth_missing=="Ethnicity missing") %>% select(age) %>% summary()
missing_eth_test %>% filter(eth_missing=="Ethnicity recorded") %>% select(age) %>% summary()

#number of unique participants in missing ethnicity/no missing groups
df %>% filter(is.na(ethnicity_5)) %>% distinct(patid) %>% nrow()
df %>% filter(!is.na(ethnicity_5)) %>% distinct(patid) %>% nrow()

##----outcomes - table s3----

###----follow-up----

#getting IQR of follow-up
summary(df$follow_up)
summary((df%>% filter(dementia==1))$follow_up)
summary((df%>% filter(dementia==0))$follow_up)

###----readmissions----

#readmission, overall
outcomes %>% filter(readmission==1) %>% 
  left_join(df %>% select(patid, indexdate, dementia), by=c("patid", "indexdate")) %>%
  distinct(patid, indexdate) %>%
  nrow()

#readmission among those without dementia
outcomes %>% filter(readmission==1) %>% 
  left_join(df %>% select(patid, indexdate, dementia), by=c("patid", "indexdate")) %>%
  filter(dementia==0) %>%
  distinct(patid, indexdate) %>%
  nrow()

#readmission among those with dementia
outcomes %>% filter(readmission==1) %>% 
  left_join(df %>% select(patid, indexdate, dementia), by=c("patid", "indexdate")) %>%
  filter(dementia==1) %>%
  distinct(patid, indexdate) %>%
  nrow()

###----discharge destinations----

#discharge destination
outcomes %>%left_join(df %>% select(patid, indexdate, dementia), by=c("patid", "indexdate")) %>%
  filter(spell_n==1) %>% 
  group_by(disdest) %>%
  summarise(n=n(), dementia=sum(dementia))

#discharged to NEW care home
new_care_home %>%left_join(df %>% select(patid, indexdate, dementia), by=c("patid", "indexdate")) %>%
  filter(spell_n==1) %>% 
  group_by(new_care_home) %>%
  summarise(n=n(), dementia=sum(dementia))



###----length of index admission----

#overall
outcomes %>% filter(spell_n==1) %>%
  left_join(df %>% select(patid, indexdate, dementia,), by=c("patid", "indexdate")) %>%
  select(los) %>%
  summary()

#among those with dementia
outcomes %>% filter(spell_n==1) %>%
  left_join(df %>% select(patid, indexdate, dementia,), by=c("patid", "indexdate")) %>%
  filter(dementia==1) %>%
  select(los) %>%
  summary()

#among those without dementia
outcomes %>% filter(spell_n==1) %>%
  left_join(df %>% select(patid, indexdate, dementia,), by=c("patid", "indexdate")) %>%
  filter(dementia==0) %>%
  select(los) %>%
  summary()



###----length of readmission----

#length of readmission; spell_n==2
#overall
outcomes %>% filter(spell_n==2) %>%
  left_join(df %>% select(patid, indexdate, dementia,), by=c("patid", "indexdate")) %>%
  select(los) %>%
  summary()

#among those with dementia
outcomes %>% filter(spell_n==2) %>%
  left_join(df %>% select(patid, indexdate, dementia,), by=c("patid", "indexdate")) %>%
  filter(dementia==1) %>%
  select(los) %>%
  summary()

#among those without dementia
outcomes %>% filter(spell_n==2) %>%
  left_join(df %>% select(patid, indexdate, dementia,), by=c("patid", "indexdate")) %>%
  filter(dementia==0) %>%
  select(los) %>%
  summary()


###----number of readmissions----

#when reporting number of readmissions, subtract 1 becuase the output of the following code includes index admission
n_spells_df <- outcomes %>%
  left_join(df %>% select(patid, indexdate, dementia), by=c("patid", "indexdate")) %>%
  group_by(patid, indexdate) %>% summarise(n_spells=max(spell_n),
                                           dementia=max(dementia))

#overall
summary(n_spells_df$n_spells)
sd(n_spells_df$n_spells)

#among those without dementia
summary(n_spells_df %>% filter(dementia == 0) %>% pull(n_spells))
sd(n_spells_df %>% filter(dementia == 0) %>% pull(n_spells))

#among those with dementia
summary(n_spells_df %>% filter(dementia == 1) %>% pull(n_spells))
sd(n_spells_df %>% filter(dementia == 1) %>% pull(n_spells))

###---- 2+ readmissions----

#overall
outcomes %>% filter(spell_n>2) %>% 
  distinct(patid, indexdate) %>%  nrow()

#with dementia
outcomes %>% filter(spell_n>2) %>% 
  distinct(patid, indexdate) %>% 
  left_join(df %>% select(patid, indexdate, dementia), by=c("patid","indexdate")) %>%
  filter(dementia==1) %>% nrow()

#without dementia
outcomes %>% filter(spell_n>2) %>% 
  distinct(patid, indexdate) %>% 
  left_join(df %>% select(patid, indexdate, dementia), by=c("patid","indexdate")) %>%
  filter(dementia==0) %>% nrow()



###----death----

#overall
outcomes %>%
  left_join(df %>% select(patid, indexdate, dementia, cprd_ddate), by=c("patid", "indexdate")) %>% 
  mutate(death_in_fu=if_else((reg_date_of_death-indexdate)<=365.25 | 
                         (is.na(reg_date_of_death) & (cprd_ddate-indexdate)<=365.25) |
                         is.na(reg_date_of_death) & is.na(cprd_ddate) & disdest=="Died", 1, 0)) %>%
  mutate(death_in_fu=if_else(is.na(death_in_fu), 0, death_in_fu)) %>% 
  group_by(patid, indexdate) %>%
  summarise(dementia=max(dementia),
            death_in_fu=max(death_in_fu)) %>%
  filter(death_in_fu==1) %>%nrow()

#among those without dementia
outcomes %>%
  left_join(df %>% select(patid, indexdate, dementia, cprd_ddate), by=c("patid", "indexdate")) %>% 
  mutate(death_in_fu=if_else((reg_date_of_death-indexdate)<=365.25 | 
                               (is.na(reg_date_of_death) & (cprd_ddate-indexdate)<=365.25) |
                               is.na(reg_date_of_death) & is.na(cprd_ddate) & disdest=="Died", 1, 0)) %>%
  mutate(death_in_fu=if_else(is.na(death_in_fu), 0, death_in_fu)) %>% 
  group_by(patid, indexdate) %>%
  summarise(dementia=max(dementia),
            death_in_fu=max(death_in_fu)) %>%
  filter(death_in_fu==1) %>%
  filter(dementia==0) %>%nrow()

#among those with dementia
outcomes %>%
  left_join(df %>% select(patid, indexdate, dementia, cprd_ddate), by=c("patid", "indexdate")) %>% 
  mutate(death_in_fu=if_else((reg_date_of_death-indexdate)<=365.25 | 
                               (is.na(reg_date_of_death) & (cprd_ddate-indexdate)<=365.25) |
                               is.na(reg_date_of_death) & is.na(cprd_ddate) & disdest=="Died", 1, 0)) %>%
  mutate(death_in_fu=if_else(is.na(death_in_fu), 0, death_in_fu)) %>% 
  group_by(patid, indexdate) %>%
  summarise(dementia=max(dementia),
            death_in_fu=max(death_in_fu)) %>%
  filter(death_in_fu==1) %>%
  filter(dementia==1) %>%nrow()




  