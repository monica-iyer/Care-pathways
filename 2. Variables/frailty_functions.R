
## Care pathway following hospitalisation study
## CPRD
##
## Author: Monica Iyer
## Date:   May 2026
##
## This script contains functions to create variables used in the electronic frailty index 2 (efi2) algorithm,
## The functions here are for the frailty components with specific instructions
## all other components use either create_var_pre_index_medcode or create_var_pre_index_medcode_5 functions in functions.R

## All functions return a dataframe with (patid, indexdate, [variable])
##
## Component scores and weights are taken directly from the eFI2 algorithm instructions.
## The final score function (alg_efi2_excl_dementia) excludes dementia as it is the primary exposure in this study.
##


##------- ACTIVITY LIMITATION----

# Flags patients with an activity limitation code on or before index date.
# Two code types:
#   Simple codes : no other instructions, any time on or before index
#   Value codes  : have numeric value requirements (Barthel index ranges,
#                  per eFI2 instructions: 0-18 or 21-90)

alg_activity_limitation <- function(study_pop,
                                    obs_files,
                                    codelists_aurum,
                                    study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "activity_limitation") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Value codes: Barthel index ranges per eFI2 instructions (0-18 or 21-90)
  codes_with_instruction <- codelist %>%
    filter(!is.na(Otherinstructions) & Otherinstructions != "")
  
  flag_with_instruction <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_instruction$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    filter(value >= 0 & value <= 18 | value >= 21 & value <= 90) %>%  # Barthel index ranges per eFI2
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_instruction = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,           by = c("patid", "indexdate")) %>%
    left_join(flag_with_instruction, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple           = replace_na(flag_simple,           0L),
      flag_with_instruction = replace_na(flag_with_instruction, 0L),
      activity_limitation   = if_else(
        flag_simple == 1L | flag_with_instruction == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, activity_limitation)
  
}


##---- ALCOHOL----

# Assigns alcohol category based on the highest-risk category recorded
# within 5 years on or before index date.
# Categories (ordered lowest to highest risk):
#   zero, previous_harmful, lower, higher, harmful
# Patients with no qualifying codes are assigned "missing".

alg_alcohol <- function(study_pop,
                        obs_files,
                        codelists_aurum,
                        study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "alcohol") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  category_order <- c("zero", "previous_harmful", "lower", "higher", "harmful")
  
  # Split codelist into numeric codes and category codes
  numeric_codes <- codelist %>%
    filter(grepl("numeric", Otherinstructions, ignore.case = TRUE)) %>%
    pull(medcodeid) %>%
    as.character()
  
  category_codes <- codelist %>%
    filter(!grepl("numeric", Otherinstructions, ignore.case = TRUE)) %>%
    mutate(
      medcodeid = as.character(medcodeid),
      category  = case_when(
        grepl("harmful drinking", Otherinstructions, ignore.case = TRUE) ~ "harmful",
        grepl("higher risk",      Otherinstructions, ignore.case = TRUE) ~ "higher",
        grepl("lower risk",       Otherinstructions, ignore.case = TRUE) ~ "lower",
        grepl("previous",         Otherinstructions, ignore.case = TRUE) ~ "previous_harmful",
        grepl("zero",             Otherinstructions, ignore.case = TRUE) ~ "zero",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(category)) %>%
    select(medcodeid, category, term)
  
  all_relevant_codes <- c(numeric_codes, category_codes$medcodeid)
  
  # loading the entire observation file into memory
  events_raw <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    
    filter(medcodeid %in% all_relevant_codes) %>%
    mutate(obsdate   = as.Date(obsdate, "%d/%m/%Y"),
           medcodeid = as.character(medcodeid)) %>%
    collect()
  
  # Apply 5-year window on or before index date
  events_windowed <- events_raw %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(obsdate <= indexdate,
           indexdate - obsdate <= 1826.25,
           obsdate <= study_end)
  
  # Category codes: assign alcohol category from codelist
  cat_events <- events_windowed %>%
    filter(medcodeid %in% category_codes$medcodeid) %>%
    left_join(category_codes, by = "medcodeid") %>%
    select(patid, indexdate, obsdate, category)
  
  # Numeric codes: derive category from weekly units
  numeric_events <- events_windowed %>%
    filter(medcodeid %in% numeric_codes) %>%
    left_join(
      codelist %>%
        mutate(medcodeid = as.character(medcodeid)) %>%
        select(medcodeid, Otherinstructions),
      by = "medcodeid"
    ) %>%
    mutate(
      weekly_units = if_else(
        grepl("multiply by 7", Otherinstructions, ignore.case = TRUE),
        as.numeric(value) * 7,
        as.numeric(value)
      ),
      category = case_when(
        weekly_units == 0  ~ "zero",
        weekly_units <= 20 ~ "lower",
        weekly_units <= 48 ~ "higher",
        weekly_units >= 49 ~ "harmful",
        TRUE               ~ NA_character_
      )
    ) %>%
    filter(!is.na(category)) %>%
    select(patid, indexdate, obsdate, category)
  
  # Combine all categorised events and take highest-risk category per patient-episode
  all_events <- bind_rows(cat_events, numeric_events) %>%
    mutate(category = factor(category, levels = category_order, ordered = TRUE))
  
  highest_category <- all_events %>%
    group_by(patid, indexdate) %>%
    summarise(alcohol = as.character(max(category)), .groups = "drop")
  
  # Join back to study population; assign "missing" where no codes found
  study_pop %>%
    select(patid, indexdate) %>%
    left_join(highest_category, by = c("patid", "indexdate")) %>%
    mutate(alcohol = if_else(is.na(alcohol), "missing", alcohol))
  
}



#---- ATRIAL FIBRILLATION----

# Flags patients with an atrial fibrillation code on or before index date.
# Three code types:
#   Simple codes    : no other instructions, any time on or before index
#   Value codes     : have numeric value requirements (value > 0, per eFI2 instructions)
#   Age-limit codes : only flagged if patient was at or above the specified age at the event

alg_atrial_fibrillation <- function(study_pop,
                                    obs_files,
                                    codelists_aurum,
                                    study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "atrial_fibrillation") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) | Otherinstructions == "") %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Value codes: flag if value > 0 on or before index (per eFI2 instructions)
  codes_with_instruction <- codelist %>%
    filter(!is.na(Otherinstructions) & Otherinstructions != "")
  
  flag_with_instruction <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_instruction$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    filter(value > 0) %>%   # per eFI2 instructions
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_instruction = 1L)
  
  # Age-limit codes: flag only if patient was at or above the specified age at the event
  codes_with_age_limit <- codelist %>%
    filter(!is.na(AgeLimit) & AgeLimit != "")
  
  flag_with_age_limit <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_age_limit$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    select(patid, event_date = obsdate, medcodeid) %>%
    left_join(study_pop %>% select(patid, indexdate, yob), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    mutate(age_at_deficit_event = as.numeric(event_date - yob) / 365.25) %>%
    left_join(codes_with_age_limit %>% select(medcodeid, AgeLimit), by = "medcodeid") %>%
    filter(age_at_deficit_event >= as.numeric(AgeLimit)) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_age_limit = 1L)
  
  # Combine: flag if any criterion met
  study_pop %>%
    left_join(flag_simple,           by = c("patid", "indexdate")) %>%
    left_join(flag_with_instruction, by = c("patid", "indexdate")) %>%
    left_join(flag_with_age_limit,   by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple           = replace_na(flag_simple,           0L),
      flag_with_instruction = replace_na(flag_with_instruction, 0L),
      flag_age_limit        = replace_na(flag_age_limit,        0L),
      atrial_fibrillation   = if_else(
        flag_simple == 1L | flag_with_instruction == 1L | flag_age_limit == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, atrial_fibrillation)
  
}


#==== BMI====

# BMI algorithm from JM code with updates from GGL to add prioritisation:
# https://github.com/julianmatthewman/psoriasis_dementia_public/blob/main/R/variables.R
#
# BMI defined pragmatically using the status recorded closest to cohort entry date.
# Prioritisation (Krishnan B group):
#   1st best : nearest status from -1y to +1 month from index
#   2nd best : nearest up to 1y after index
#   3rd best : any before index
#   4th best : any after index
#
# Add medcodes with a BMI category in Otherinstructions (in the codelist), used only where a numeric BMI cannot be calculated.

alg_bmi <- function(study_pop,
                    obs_files,
                    codelists_aurum,
                    study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "bmi") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Pull all BMI-related observations.
  bmi_eventdata <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, numunitid, value) %>%
    filter(medcodeid %in% codelist$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    left_join(codelist %>% select(medcodeid, Otherinstructions), by = "medcodeid") %>%
    collect()
  
  # Separate into height, weight, and recorded BMI subsets by unit
  height_in_cm <- bmi_eventdata %>% fsubset(numunitid == 122)               %>% fsubset(value > 100 & value < 250)
  height_in_m  <- bmi_eventdata %>% fsubset(numunitid == 173)               %>% fsubset(value > 1   & value < 2.5)
  height_in_m  <- bind_rows(height_in_m, mutate(height_in_cm, value = value / 100))
  
  weight_in_kg <- bmi_eventdata %>% fsubset(numunitid %in% c(156, 827, 10185))          %>% fsubset(value > 30  & value < 300)
  bmi_recorded <- bmi_eventdata %>% fsubset(numunitid %in% c(157, 657, 907, 108, 568, 1309)) %>% fsubset(value > 10  & value < 60)
  
  # Pick height: prioritise by time since index date, take closest
  # Note: height assumed constant so no grouping by episode (fgroup_by on patid, indexdate
  # but first() takes the closest record regardless of episode)
  height <- study_pop[c("patid", "indexdate")] %>%
    join(height_in_m[c("patid", "value", "obsdate")], on = "patid", how = "left", multiple = TRUE) %>%
    ftransform(distance = as.numeric(obsdate - indexdate)) %>%
    mutate(
      priority    = case_when(
        distance >= -365 & distance <= 30  ~ 1,
        distance < -365                    ~ 2,
        distance > 30   & distance <= 365  ~ 3,
        distance > 365                     ~ 4,
        TRUE                               ~ NA_integer_
      ),
      absdistance = abs(distance)
    ) %>%
    roworder(patid, priority, absdistance) %>%
    fgroup_by(patid, indexdate) %>%
    ffirst() %>%
    fselect(patid, indexdate, height_in_m = value)
  
  # Pick weight: same prioritisation as height
  weight <- study_pop[c("patid", "indexdate")] %>%
    join(weight_in_kg[c("patid", "value", "obsdate")], on = "patid", how = "left", multiple = TRUE) %>%
    ftransform(distance = as.numeric(obsdate - indexdate)) %>%
    mutate(
      priority    = case_when(
        distance >= -365 & distance <= 30  ~ 1,
        distance < -365                    ~ 2,
        distance > 30   & distance <= 365  ~ 3,
        distance > 365                     ~ 4,
        TRUE                               ~ NA_integer_
      ),
      absdistance = abs(distance)
    ) %>%
    roworder(patid, priority, absdistance) %>%
    fgroup_by(patid, indexdate) %>%
    fsummarise(value = first(value)) %>%
    fungroup() %>%
    fselect(patid, indexdate, weight_in_kg = value)
  
  # Pick recorded BMI: same prioritisation
  bmi_recorded <- study_pop[c("patid", "indexdate")] %>%
    join(bmi_recorded[c("patid", "value", "obsdate")], on = "patid", how = "left", multiple = TRUE) %>%
    ftransform(distance = as.numeric(obsdate - indexdate)) %>%
    mutate(
      priority    = case_when(
        distance >= -365 & distance <= 30  ~ 1,
        distance < -365                    ~ 2,
        distance > 30   & distance <= 365  ~ 3,
        distance > 365                     ~ 4,
        TRUE                               ~ NA_integer_
      ),
      absdistance = abs(distance)
    ) %>%
    roworder(patid, priority, absdistance) %>%
    fgroup_by(patid, indexdate) %>%
    fsummarise(value = first(value)) %>%
    fungroup() %>%
    select(patid, indexdate, bmi_recorded = value)
  
  # Calculate BMI from height and weight
  bmi_calc <- height %>%
    join(weight, on = c("patid", "indexdate"), how = "left", multiple = TRUE) %>%
    ftransform(bmi_calc = weight_in_kg / (height_in_m^2)) %>%
    fsubset(!is.na(bmi_calc)) %>%
    fselect(patid, indexdate, bmi_calc)
  
  # Take calculated BMI first; use recorded BMI as second best
  bmi_final <- full_join(bmi_recorded, bmi_calc, by = c("patid", "indexdate")) %>%
    mutate(bmi = coalesce(bmi_calc, bmi_recorded)) %>%
    select(patid, indexdate, bmi)
  
  # Derive BMI category from medcodes with a category in Otherinstructions
  # (used only where numeric BMI cannot be calculated — same prioritisation)
  bmi_from_medcode <- study_pop[c("patid", "indexdate")] %>%
    join(
      bmi_eventdata %>%
        filter(!is.na(Otherinstructions) & Otherinstructions != "") %>%
        select(patid, obsdate, Otherinstructions),
      on = "patid", how = "left", multiple = TRUE
    ) %>%
    ftransform(distance = as.numeric(obsdate - indexdate)) %>%
    mutate(
      priority    = case_when(
        distance >= -365 & distance <= 30  ~ 1,
        distance < -365                    ~ 2,
        distance > 30   & distance <= 365  ~ 3,
        distance > 365                     ~ 4,
        TRUE                               ~ NA_integer_
      ),
      absdistance = abs(distance)
    ) %>%
    roworder(patid, priority, absdistance) %>%
    fgroup_by(patid, indexdate) %>%
    fsummarise(bmi_classification_medcode = first(Otherinstructions)) %>%
    fungroup()
  
  # Merge medcode-derived category; only used where numeric BMI is missing
  bmi_final <- bmi_final %>%
    left_join(bmi_from_medcode, by = c("patid", "indexdate")) %>%
    mutate(bmi_classification_source = case_when(
      !is.na(bmi)                        ~ "calculated",
      !is.na(bmi_classification_medcode) ~ "medcode",
      TRUE                               ~ NA_character_
    ))
  
  # Derive final BMI classification
  study_pop %>%
    left_join(bmi_final, by = c("patid", "indexdate")) %>%
    mutate(bmi_classification = case_when(
      bmi < 18.5              ~ "Underweight",
      bmi >= 18.5 & bmi < 25 ~ "Recommended weight",
      bmi >= 25   & bmi < 30 ~ "Overweight",
      bmi >= 30               ~ "Obese",
      !is.na(bmi) ~ bmi_classification_medcode,
      TRUE                               ~ NA_character_
    )) %>%
    select(patid, indexdate, bmi, bmi_classification)
  
}



#---- COGNITIVE IMPAIRMENT----

# Flags patients with a cognitive impairment code on or before index date.
# Two code types:
#   Simple codes : no other instructions, any time on or before index
#   Value codes  : cognitive test score >= 8, per eFI2 instructions

alg_cognitive_impairment <- function(study_pop,
                                     obs_files,
                                     codelists_aurum,
                                     study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "cognitive_impairment") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Value codes: cognitive test score >= 8, per eFI2 instructions
  codes_with_instruction <- codelist %>%
    filter(!is.na(Otherinstructions) & Otherinstructions != "")
  
  flag_with_instruction <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_instruction$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    filter(value >= 8) %>%   # cognitive test score cutoff per eFI2 instructions
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_instruction = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,           by = c("patid", "indexdate")) %>%
    left_join(flag_with_instruction, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple           = replace_na(flag_simple,           0L),
      flag_with_instruction = replace_na(flag_with_instruction, 0L),
      cognitive_impairment  = if_else(
        flag_simple == 1L | flag_with_instruction == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, cognitive_impairment)
  
}



#---- COPD----

# Flags patients with a COPD code on or before index date.
# Two code types:
#   Simple codes : no other instructions, any time on or before index
#   Value codes  : flag if value > 0, per eFI2 instructions

alg_copd <- function(study_pop,
                     obs_files,
                     codelists_aurum,
                     study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "copd") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Value codes: flag if value > 0, per eFI2 instructions
  codes_with_instruction <- codelist %>%
    filter(!is.na(Otherinstructions) & Otherinstructions != "")
  
  flag_with_instruction <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_instruction$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    filter(value > 0) %>%   # per eFI2 instructions
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_instruction = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,           by = c("patid", "indexdate")) %>%
    left_join(flag_with_instruction, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple           = replace_na(flag_simple,           0L),
      flag_with_instruction = replace_na(flag_with_instruction, 0L),
      copd                  = if_else(
        flag_simple == 1L | flag_with_instruction == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, copd)
  
}


#---- FALLS ----


# Flags patients with a falls code on or before index date.
#all codes must be on or after age 55
# Three code types:
#   Simple codes         : no other instructions (other than at least 55 years old), any time on or before index
#   Time-constrained     : within 5 years on or before index (per eFI2 instructions)
#   Value codes (">0")   : flag if recorded value > 0, any time on or before index
#                          (per eFI2 instructions)
#
#


alg_falls <- function(study_pop,
                           obs_files,
                           codelists_aurum,
                           study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "falls_for_efi2") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    filter(medcodeid %in% codes_simple) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate, yob), by = "patid") %>%
    collect() %>%
    mutate(age_at_deficit_event = as.numeric(event_date - yob) / 365.25) %>%
    filter(event_date <= indexdate) %>%
    filter(age_at_deficit_event>=55) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Time-constrained codes: within 5 years on or before index date
  codes_with_time_constraints <- codelist %>%
    filter(!is.na(TimeConstraintyears) & TimeConstraintyears != "")
  
  flag_with_time_constraints <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_time_constraints$medcodeid) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    left_join(codes_with_time_constraints %>% select(medcodeid, TimeConstraintyears), by = "medcodeid") %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate, yob), by = "patid") %>%
    collect() %>%
    filter(event_date <= indexdate & as.numeric(indexdate - event_date) <= 1826.25) %>%
    mutate(age_at_deficit_event = as.numeric(event_date - yob) / 365.25) %>%
    filter(age_at_deficit_event>=55) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_time_constraints = 1L)
  
  # Value codes (Otherinstructions == ">0"): flag if value > 0, on or before index
  codes_with_value <- codelist %>%
    filter(Otherinstructions == ">0") %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_with_value <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    filter(medcodeid %in% codes_with_value) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(value > 0) %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate, yob), by = "patid") %>%
    collect() %>%
  filter(event_date <= indexdate) %>%
    mutate(age_at_deficit_event = as.numeric(event_date - yob) / 365.25) %>%
    filter(age_at_deficit_event>=55) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_value = 1L)
  
  # Combine: flag if any criterion met
  study_pop %>%
    left_join(flag_simple,                by = c("patid", "indexdate")) %>%
    left_join(flag_with_time_constraints, by = c("patid", "indexdate")) %>%
    left_join(flag_with_value,            by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple                = replace_na(flag_simple,                0L),
      flag_with_time_constraints = replace_na(flag_with_time_constraints, 0L),
      flag_with_value            = replace_na(flag_with_value,            0L),
      falls                 = if_else(
        flag_simple == 1L | flag_with_time_constraints == 1L | flag_with_value == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, falls)
  
}




#---- HOUSEBOUND----

# Flags patients with a housebound code on or before index date.
# Two code types:
#   Simple codes     : no other instructions, any time on or before index
#   Time-constrained : within 5 years on or before index (per eFI2 instructions)

alg_housebound <- function(study_pop,
                           obs_files,
                           codelists_aurum,
                           study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "housebound") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Time-constrained codes: within 5 years on or before index date
  codes_with_time_constraints <- codelist %>%
    filter(!is.na(TimeConstraintyears) & TimeConstraintyears != "")
  
  flag_with_time_constraints <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_time_constraints$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    left_join(codes_with_time_constraints %>% select(medcodeid, TimeConstraintyears), by = "medcodeid") %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate & as.numeric(indexdate - event_date) <= 1826.25) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_time_constraints = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,                by = c("patid", "indexdate")) %>%
    left_join(flag_with_time_constraints, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple                = replace_na(flag_simple,                0L),
      flag_with_time_constraints = replace_na(flag_with_time_constraints, 0L),
      housebound                 = if_else(
        flag_simple == 1L | flag_with_time_constraints == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, housebound)
  
}


#---- HYPOTENSION / SYNCOPE----

# Flags patients meeting any of four criteria on or before index date:
#   1. Direct binary medcode match (no time constraint): any time on or before index
#   2. Direct binary medcode match (time constrained):   within 5 years on or before index
#   3. Systolic BP  : 3+ readings < 90 mmHg on separate dates, within 5 years
#   4. Diastolic BP : 3+ readings < 60 mmHg on separate dates, within 5 years

alg_hypotension_syncope <- function(study_pop,
                                    obs_files,
                                    codelists_aurum,
                                    study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "hypotension_syncope") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Split codelist into groups
  direct_codes_no_time_const <- codelist %>%
    filter(is.na(Otherinstructions) | Otherinstructions == "") %>%
    filter(is.na(TimeConstraintyears)) %>%
    pull(medcodeid)
  
  direct_codes_time_const <- codelist %>%
    filter(is.na(Otherinstructions) | Otherinstructions == "") %>%
    filter(!is.na(TimeConstraintyears)) %>%
    pull(medcodeid)
  
  systolic_codes  <- codelist %>% filter(grepl("systolic",  Otherinstructions, ignore.case = TRUE)) %>% pull(medcodeid)
  diastolic_codes <- codelist %>% filter(grepl("diastolic", Otherinstructions, ignore.case = TRUE)) %>% pull(medcodeid)
  
  all_relevant_codes <- c(direct_codes_no_time_const, direct_codes_time_const,
                          systolic_codes, diastolic_codes)
  
  events_raw <- open_dataset(obs_files) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    
    filter(medcodeid %in% all_relevant_codes) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect()
  
  # Direct medcode match — no time constraint: any time on or before index
  direct_flag_no_time_const <- events_raw %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(medcodeid %in% direct_codes_no_time_const) %>%
    filter(obsdate <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(direct_hypotension = TRUE)
  
  # Direct medcode match — time constrained: within 5 years on or before index
  direct_flag_time_const <- events_raw %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(medcodeid %in% direct_codes_time_const) %>%
    filter(obsdate <= indexdate & as.numeric(indexdate - obsdate) <= 1826.25) %>%
    distinct(patid, indexdate) %>%
    mutate(direct_hypotension = TRUE)
  
  direct_flag <- bind_rows(direct_flag_no_time_const, direct_flag_time_const) %>%
    distinct(patid, indexdate) %>%
    mutate(direct_hypotension = TRUE)
  
  # Systolic: 3+ readings < 90 mmHg on separate dates, within 5 years on or before index
  systolic_flag <- events_raw %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(medcodeid %in% systolic_codes) %>%
    mutate(value = as.numeric(value)) %>%
    filter(!is.na(value), value < 90) %>%
    filter(obsdate <= indexdate & as.numeric(indexdate - obsdate) <= 1826.25) %>%
    distinct(patid, indexdate, obsdate) %>%
    count(patid, indexdate, name = "n_systolic") %>%
    filter(n_systolic >= 3) %>%
    mutate(systolic_hypotension = TRUE)
  
  # Diastolic: 3+ readings < 60 mmHg on separate dates, within 5 years on or before index
  
  diastolic_flag <- events_raw %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(medcodeid %in% diastolic_codes) %>%
    mutate(value = as.numeric(value)) %>%
    filter(!is.na(value), value < 60) %>%
    filter(obsdate <= indexdate & as.numeric(indexdate - obsdate) <= 1826.25) %>%
    distinct(patid, indexdate, obsdate) %>%
    count(patid, indexdate, name = "n_diastolic") %>%
    filter(n_diastolic >= 3) %>%
    mutate(diastolic_hypotension = TRUE)
  
  study_pop %>%
    select(patid, indexdate) %>%
    left_join(direct_flag,    by = c("patid", "indexdate")) %>%
    left_join(systolic_flag,  by = c("patid", "indexdate")) %>%
    left_join(diastolic_flag, by = c("patid", "indexdate")) %>%
    mutate(
      hypotension_syncope = if_else(
        coalesce(direct_hypotension,   FALSE) |
          coalesce(systolic_hypotension,  FALSE) |
          coalesce(diastolic_hypotension, FALSE),
        1L, 0L
      )
    ) %>%
    select(patid, indexdate, hypotension_syncope)
  
}


#----- LIVER PROBLEMS----

# Flags patients with a liver problems code on or before index date.
# Two code types:
#   Simple codes     : no other instructions, any time on or before index
#   Time-constrained : within 5 years on or before index (per eFI2 instructions)
#

alg_liver_problems <- function(study_pop,
                               obs_files,
                               codelists_aurum,
                               study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "liver_problems") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Time-constrained codes: within 5 years on or before index date
  codes_with_time_constraints <- codelist %>%
    filter(!is.na(TimeConstraintyears) & TimeConstraintyears != "")
  
  flag_with_time_constraints <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_time_constraints$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    left_join(codes_with_time_constraints %>% select(medcodeid, TimeConstraintyears), by = "medcodeid") %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate & as.numeric(indexdate - event_date) <= 1826.25) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_time_constraints = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,                by = c("patid", "indexdate")) %>%
    left_join(flag_with_time_constraints, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple                = replace_na(flag_simple,                0L),
      flag_with_time_constraints = replace_na(flag_with_time_constraints, 0L),
      liver_problems             = if_else(
        flag_simple == 1L | flag_with_time_constraints == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, liver_problems)
  
}



#---- PERIPHERAL VASCULAR DISEASE----

# Flags patients with a peripheral vascular disease code on or before index date.
# Two code types:
#   Simple codes : no other instructions, any time on or before index
#   Value codes  : ankle-brachial index < 0.95, per eFI2 instructions

alg_peripheral_vascular_disease <- function(study_pop,
                                            obs_files,
                                            codelists_aurum,
                                            study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "peripheral_vascular_disease") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Value codes: ankle-brachial index < 0.95, per eFI2 instructions
  codes_with_instruction <- codelist %>%
    filter(!is.na(Otherinstructions) & Otherinstructions != "")
  
  flag_with_instruction <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_instruction$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    filter(value < 0.95) %>%   # ankle-brachial index cutoff per eFI2 instructions
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_instruction = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,           by = c("patid", "indexdate")) %>%
    left_join(flag_with_instruction, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple                 = replace_na(flag_simple,           0L),
      flag_with_instruction       = replace_na(flag_with_instruction, 0L),
      peripheral_vascular_disease = if_else(
        flag_simple == 1L | flag_with_instruction == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, peripheral_vascular_disease)
  
}


#---- POLYPHARMACY----

# Counts distinct BNF chapters of medications prescribed in the 90 days
# on or before index date, and assigns a polypharmacy category:
#   "0-4"  : fewer than 5 distinct chapters (no polypharmacy)
#   "5-9"  : 5-9 distinct chapters
#   "10+"  : 10 or more distinct chapters
#
# Uses the Aurum product codebrowser to map prodcodeids to BNF chapters.
# Patients with no matching prescriptions are recoded to "0-4".

alg_polypharmacy <- function(study_pop,
                             drug_files,
                             polypharmacy_codelist,
                             study_end) {
  
  polypharmacy_codelist <- read_tsv(polypharmacy_codelist,
                                    col_types = cols(prodcodeid = col_character(),
                                                     bnfchapter=col_character()))
  polypharmacy_codelist_arrow <- arrow_table(polypharmacy_codelist)
  
  #break study_pop into chunks of 500,000 and loop through each chunk
  chunks <- split(study_pop, ceiling(seq_len(nrow(study_pop)) / 500000))
  
  process_chunk <- function(chunk) {
    open_dataset(drug_files) %>%
      select(patid, issuedate, prodcodeid) %>%
      filter(prodcodeid %in% polypharmacy_codelist$prodcodeid) %>%
      left_join(polypharmacy_codelist_arrow, by = "prodcodeid") %>%
      mutate(issuedate = as.Date(issuedate, "%d/%m/%Y")) %>%  
      filter(issuedate <= study_end, issuedate >= as.Date("2009-09-30")) %>% #filtering to 3 months before first possible index
      inner_join(arrow_table(chunk %>% select(patid, indexdate)), by = "patid") %>%
      filter(issuedate <= indexdate) %>%
      mutate(broad_chapter = str_sub(bnfchapter, 1, 4)) %>%
      collect() %>%
      # Filter to records within 90 days on or before index date
      filter(as.numeric(indexdate - issuedate) < 90) %>%
      # Create broad chapter variable matching the codelist level (first 4 characters)
      group_by(patid, indexdate) %>%
      summarise(n_chapters = n_distinct(broad_chapter), .groups = "drop") %>%
      mutate(polypharmacy = case_when(
        n_chapters < 5                    ~ "0-4",
        n_chapters >= 5 & n_chapters < 10 ~ "5-9",
        n_chapters >= 10                  ~ "10+"))
  }
  temp <- purrr::map(chunks, process_chunk) %>%
    bind_rows()
  
  study_pop %>%
    left_join(temp, by = c("patid", "indexdate")) %>%
    # Recode missing to "0-4" (no prescriptions found = no polypharmacy)
    mutate(polypharmacy = ifelse(is.na(polypharmacy), "0-4", polypharmacy)) %>%
    select(patid, indexdate, polypharmacy)
  
}



#---- RESPIRATORY DISEASE----

# Flags patients with a respiratory disease code on or before index date.
# Two code types:
#   Simple codes     : no other instructions, any time on or before index
#   Time-constrained : within 5 years on or before index (per eFI2 instructions)

alg_respiratory_disease <- function(study_pop,
                                    obs_files,
                                    codelists_aurum,
                                    study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "respiratory_disease") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Time-constrained codes: within 5 years on or before index date
  codes_with_time_constraints <- codelist %>%
    filter(!is.na(TimeConstraintyears) & TimeConstraintyears != "")
  
  flag_with_time_constraints <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_time_constraints$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    left_join(codes_with_time_constraints %>% select(medcodeid, TimeConstraintyears), by = "medcodeid") %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate & as.numeric(indexdate - event_date) <= 1826.25) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_time_constraints = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,                by = c("patid", "indexdate")) %>%
    left_join(flag_with_time_constraints, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple                = replace_na(flag_simple,                0L),
      flag_with_time_constraints = replace_na(flag_with_time_constraints, 0L),
      respiratory_disease        = if_else(
        flag_simple == 1L | flag_with_time_constraints == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, respiratory_disease)
  
}



#---- SKIN ULCER----

# Flags patients with a skin ulcer code on or before index date.
# Two code types:
#   Simple codes : no other instructions, any time on or before index
#   Value codes  : flag if value > 0, per eFI2 instructions

alg_skin_ulcer <- function(study_pop,
                           obs_files,
                           codelists_aurum,
                           study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "skin_ulcer") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Value codes: flag if value > 0, per eFI2 instructions
  codes_with_instruction <- codelist %>%
    filter(!is.na(Otherinstructions) & Otherinstructions != "")
  
  flag_with_instruction <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid, value) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_instruction$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    filter(value > 0) %>%   # per eFI2 instructions
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_instruction = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,           by = c("patid", "indexdate")) %>%
    left_join(flag_with_instruction, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple           = replace_na(flag_simple,           0L),
      flag_with_instruction = replace_na(flag_with_instruction, 0L),
      skin_ulcer            = if_else(
        flag_simple == 1L | flag_with_instruction == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, skin_ulcer)
  
}


#---- SMOKING STATUS----

# Assigns smoking status based on the most recent code recorded on or before
# index date. If a current smoker and ex-smoker code are recorded on the same
# day, current smoker takes precedence.
# Output: smoking = "current", "ex", or "none/missing"
#
alg_smoking <- function(study_pop,
                        obs_files,
                        codelists_aurum,
                        study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "smoking") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Classify codes as current or ex smoker; exclude non-smoker codes
  codelist_clean <- codelist %>%
    mutate(
      category  = case_when(
        grepl("current", Otherinstructions, ignore.case = TRUE) ~ "current",
        grepl("ex",      Otherinstructions, ignore.case = TRUE) ~ "ex",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(category))
  
  # Pull all relevant events
  events_raw <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    filter(medcodeid %in% codelist_clean$medcodeid) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    left_join(codelist_clean %>% select(medcodeid, category), by = "medcodeid") %>%
    collect()
  
  # Get most recent code on or before index date; current beats ex on same date
  smoking_final <- events_raw %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(obsdate <= indexdate) %>%
    mutate(category_rank = if_else(category == "current", 2L, 1L)) %>%
    group_by(patid, indexdate) %>%
    slice_max(order_by = data.frame(obsdate, category_rank), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(patid, indexdate, smoking = category)
  
  # Join back to study population; assign "none/missing" where no codes found
  study_pop %>%
    select(patid, indexdate) %>%
    left_join(smoking_final, by = c("patid", "indexdate")) %>%
    mutate(smoking = if_else(is.na(smoking), "none/missing", smoking))
  
}


#---- SOCIAL VULNERABILITY----

# Flags patients with a social vulnerability code on or before index date.
# Two code types:
#   Simple codes     : no other instructions, any time on or before index
#   Time-constrained : within 5 years on or before index (per eFI2 instructions)

alg_social_vulnerability <- function(study_pop,
                                     obs_files,
                                     codelists_aurum,
                                     study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "social_vulnerability") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Time-constrained codes: within 5 years on or before index date
  codes_with_time_constraints <- codelist %>%
    filter(!is.na(TimeConstraintyears) & TimeConstraintyears != "")
  
  flag_with_time_constraints <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_time_constraints$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    left_join(codes_with_time_constraints %>% select(medcodeid, TimeConstraintyears), by = "medcodeid") %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate & as.numeric(indexdate - event_date) <= 1826.25) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_time_constraints = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,                by = c("patid", "indexdate")) %>%
    left_join(flag_with_time_constraints, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple                = replace_na(flag_simple,                0L),
      flag_with_time_constraints = replace_na(flag_with_time_constraints, 0L),
      social_vulnerability       = if_else(
        flag_simple == 1L | flag_with_time_constraints == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, social_vulnerability)
  
}




#---- WEIGHT LOSS----

# Flags patients with a weight loss code on or before index date.
# Two code types:
#   Simple codes     : no other instructions, any time on or before index
#   Time-constrained : within 5 years on or before index (per eFI2 instructions)

alg_weight_loss <- function(study_pop,
                            obs_files,
                            codelists_aurum,
                            study_end) {
  
  codelist <- codelists_aurum %>%
    filter(name == "weight_loss") %>%
    pull(full) %>%
    .[[1]]
  
  
  
  # Simple codes: no other instructions, on or before index date
  codes_simple <- codelist %>%
    filter(is.na(Otherinstructions) & is.na(AgeLimit) & is.na(TimeConstraintyears)) %>%
    pull(medcodeid) %>%
    as.character()
  
  flag_simple <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_simple) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate <= study_end) %>%
    collect() %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_simple = 1L)
  
  # Time-constrained codes: within 5 years on or before index date
  codes_with_time_constraints <- codelist %>%
    filter(!is.na(TimeConstraintyears) & TimeConstraintyears != "")
  
  flag_with_time_constraints <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% codes_with_time_constraints$medcodeid) %>%
    
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect() %>%
    left_join(codes_with_time_constraints %>% select(medcodeid, TimeConstraintyears), by = "medcodeid") %>%
    select(patid, event_date = obsdate) %>%
    left_join(study_pop %>% select(patid, indexdate), by = "patid") %>%
    filter(event_date <= indexdate & as.numeric(indexdate - event_date) <= 1826.25) %>%
    distinct(patid, indexdate) %>%
    mutate(flag_with_time_constraints = 1L)
  
  # Combine: flag if either criterion met
  study_pop %>%
    left_join(flag_simple,                by = c("patid", "indexdate")) %>%
    left_join(flag_with_time_constraints, by = c("patid", "indexdate")) %>%
    mutate(
      flag_simple                = replace_na(flag_simple,                0L),
      flag_with_time_constraints = replace_na(flag_with_time_constraints, 0L),
      weight_loss                = if_else(
        flag_simple == 1L | flag_with_time_constraints == 1L, 1L, 0L
      )
    ) %>%
    select(patid, indexdate, weight_loss)
  
}


#==== FINAL eFI2 SCORE (EXCLUDING DEMENTIA)====

# Assembles all eFI2 component variables and computes the final frailty score.
# Dementia is EXCLUDED as it is the primary exposure in this study.
#
# Score logic (verified against eFI2 instructions — see header sanity check):
#   - Binary variables     : contribute coefficient if == 1, else 0
#   - Missing BMI          : contributes a coefficient (per eFI2 instructions)
#   - Missing alcohol      : contributes a coefficient (per eFI2 instructions)
#   - Missing polypharmacy : recoded to "0-4" upstream, contributes 0
#   - Smoking              : only "current" contributes coefficient; ex and
#                            none/missing contribute 0
#   - Final score          : rowSums of all component scores
#
# All component weights taken directly from eFI2 algorithm instructions.

alg_efi2_excl_dementia <- function(activity_limitation,
                                   alcohol,
                                   atrial_fibrillation,
                                   cancer,
                                   cognitive_impairment,
                                   copd,
                                   dressing_problems,
                                   environment_problems,
                                   falls,
                                   fracture,
                                   fragility_fracture,
                                   heart_failure,
                                   housebound,
                                   hypotension_syncope,
                                   liver_problems,
                                   med_management,
                                   memory_concerns,
                                   mobility_problems,
                                   motor_neuron_disease,
                                   bmi,
                                   palliative_care,
                                   parkinsonism_tremor,
                                   peptic_ulcer,
                                   peripheral_vascular_disease,
                                   polypharmacy,
                                   requirement_for_care,
                                   respiratory_disease,
                                   seizures,
                                   self_harm,
                                   skin_ulcer,
                                   smoking,
                                   social_vulnerability,
                                   stroke,
                                   tia,
                                   weight_loss) {
  
  all_dfs <- list(
    activity_limitation,
    alcohol,
    atrial_fibrillation,
    cancer,
    cognitive_impairment,
    copd,
    dressing_problems,
    environment_problems,
    falls,
    fracture,
    fragility_fracture,
    heart_failure,
    housebound,
    hypotension_syncope,
    liver_problems,
    med_management,
    memory_concerns,
    mobility_problems,
    motor_neuron_disease,
    bmi,
    palliative_care,
    parkinsonism_tremor,
    peptic_ulcer,
    peripheral_vascular_disease,
    polypharmacy,
    requirement_for_care,
    respiratory_disease,
    seizures,
    self_harm,
    skin_ulcer,
    smoking,
    social_vulnerability,
    stroke,
    tia,
    weight_loss
  )
  
  reduce(all_dfs, left_join, by = c("patid", "indexdate")) %>%
    
    mutate(
      # If cognitive impairment is flagged, zero out memory concerns to avoid
      # double-counting (per eFI2 instructions)
      memory_concerns = if_else(cognitive_impairment == 1L, 0L, memory_concerns),
      activity_limitation_score         = if_else(activity_limitation         == 1, 0.019, 0, missing = 0),
      atrial_fibrillation_score         = if_else(atrial_fibrillation         == 1, 0.016, 0, missing = 0),
      cancer_score                      = if_else(cancer                      == 1, 0.031, 0, missing = 0),
      cognitive_impairment_score        = if_else(cognitive_impairment        == 1, 0.014, 0, missing = 0),
      copd_score                        = if_else(copd                        == 1, 0.015, 0, missing = 0),
      dressing_problems_score           = if_else(dressing_problems           == 1, 0.006, 0, missing = 0),
      environment_problems_score        = if_else(environ_problems        == 1, 0.015, 0, missing = 0),
      falls_score                       = if_else(falls                  == 1, 0.079, 0, missing = 0),
      fracture_score                    = if_else(fractures                    == 1, 0.010, 0, missing = 0),
      fragility_fracture_score          = if_else(fragility_fracture     == 1, 0.022, 0, missing = 0),
      heart_failure_score               = if_else(heart_failure               == 1, 0.014, 0, missing = 0),
      housebound_score                  = if_else(housebound                  == 1, 0.042, 0, missing = 0),
      hypotension_score                 = if_else(hypotension_syncope         == 1, 0.024, 0, missing = 0),
      liver_problems_score              = if_else(liver_problems              == 1, 0.030, 0, missing = 0),
      medication_management_score       = if_else(med_management              == 1, 0.041, 0, missing = 0),
      memory_concerns_score             = if_else(memory_concerns             == 1, 0.015, 0, missing = 0),
      mobility_problems_score           = if_else(mobility_problems           == 1, 0.060, 0, missing = 0),
      motor_neuron_disease_score        = if_else(motor_neuron_disease        == 1, 0.045, 0, missing = 0),
      palliative_care_score             = if_else(palliative_care             == 1, 0.065, 0, missing = 0),
      parkinsonism_tremor_score         = if_else(parkinsonism_tremor         == 1, 0.004, 0, missing = 0),
      peptic_ulcer_score                = if_else(peptic_ulcer_disease        == 1, 0.006, 0, missing = 0),
      peripheral_vascular_disease_score = if_else(peripheral_vascular_disease == 1, 0.003, 0, missing = 0),
      requirement_for_care_score        = if_else(requirement_for_care        == 1, 0.027, 0, missing = 0),
      respiratory_disease_score         = if_else(respiratory_disease         == 1, 0.001, 0, missing = 0),
      seizures_score                    = if_else(seizures                    == 1, 0.003, 0, missing = 0),
      self_harm_score                   = if_else(self_harm                   == 1, 0.001, 0, missing = 0),
      skin_ulcer_score                  = if_else(skin_ulcer                  == 1, 0.028, 0, missing = 0),
      social_vulnerability_score        = if_else(social_vulnerability        == 1, 0.030, 0, missing = 0),
      stroke_score                      = if_else(stroke                      == 1, 0.014, 0, missing = 0),
      tia_score                         = if_else(tia                         == 1, 0.003, 0, missing = 0),
      weight_loss_score                 = if_else(weight_loss                 == 1, 0.025, 0, missing = 0),
      
      # --- Smoking score ---
      # Ex-smokers and none/missing contribute 0.
      smoking_score = if_else(smoking == "current", 0.013, 0, missing = 0),
      
      # --- Categorical component scores (per eFI2 instructions) ---
      
      # Alcohol: harmful and previous_harmful contribute coefficients;
      # missing also contributes a coefficient (per eFI2 instructions);
      # zero, lower, and higher contribute 0.
      # category labels assigned in alg_alcohol.
      alcohol_score = case_when(
        alcohol == "harmful"          ~ 0.029,
        alcohol == "previous_harmful" ~ 0.173,
        is.na(alcohol)                ~ 0.017,
        TRUE                          ~ 0
      ),
      
      # BMI: missing and underweight contribute coefficients (per eFI2 instructions);
      # all other categories contribute 0.
      bmi_score = case_when(
        is.na(bmi_classification)           ~ 0.032,
        bmi_classification == "Underweight" ~ 0.056,
        TRUE                                ~ 0
      ),
      
      # Polypharmacy: 5-9 chapters -> 0.041; 10+ chapters -> 0.064;
      # 0-4 (including missing, recoded upstream) -> 0.
      polypharmacy_score = case_when(
        polypharmacy == "5-9" ~ 0.041,
        polypharmacy == "10+" ~ 0.064,
        TRUE                  ~ 0
      )
    ) %>%
    # Sum all component scores to produce the final frailty score
    mutate(efi2_score_excl_dementia = rowSums(across(ends_with("_score")), na.rm = TRUE)) %>%
    select(patid, indexdate, efi2_score_excl_dementia)
  
}
