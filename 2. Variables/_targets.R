
#Monica Iyer
#Care pathways following infection CPRD project
#Data management April 2026

#Uses the study_pop dataframe in 'Create' to add outcome and covariate variables to the study population
#The outputs:
#study_pop_with_variables - one row per included patid/indexdate with all covariates
#outcomes - multiple rows per included patid/index date with the outcome events

# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

##----setup----
# Load packages required to define the pipeline:
library(targets)
library(arrow)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)
library(haven)
library(readxl)
library(purrr)
library(tarchetypes) 

library(crew)


tar_option_set(
  packages = c(
    "arrow",
    "dplyr",
    "data.table",
    "tidyverse",
    "stringr",
    "haven",
    "readxl",
    "tarchetypes",
    "purrr",
    "crew"
  )
)

# PARALLELISATION SETUP


# Safely determine how many workers to use.
# Uses physical cores only (logical = FALSE avoids hyperthreads),
# leaves one core free for the main process, and caps at 8.
n_workers <- max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
message("targets will use ", n_workers, " parallel worker(s).")

controller <- crew_controller_local(
  name    = "my_controller",
  workers = n_workers
)

tar_option_set(
  controller = controller,
  storage    = "worker",
  retrieval  = "worker",
  packages   = c(
    "arrow",
    "dplyr",
    "data.table",
    "tidyverse",
    "stringr",
    "haven",
    "readxl",
    "tarchetypes",
    "purrr",
    "collapse"
  )
)


#All functions called in this pipeine are located in functions.R
tar_source("./functions.R")
tar_source("./frailty_functions.R")

#global paths folder not provided - contains only file paths and no analytic code
tar_source("./global paths.R")



list(
    # ---- Codelists -----------------------------
    
    #read in all medcode coddelists
    
    tar_target(
      codelists_aurum,
      tribble(~name,                         ~label,                               ~codevar,     ~extract_from,
              
              # --- Exposure ---
              "dementia",                    "Dementia",                           "medcodeid",  "observation",
              
              # --- Outcomes ---

              
              # --- Living situation / residence status ---
              "living_alone",                "Living alone",                       "medcodeid",  "observation",
              "care_home",                   "Care home residence",                "medcodeid",  "observation",
              
              
              # --- Other covariates ---
              "delirium",                    "Delirium",                           "medcodeid",  "observation",
              "ihd",                         "Ischaemic heart disease",            "medcodeid",  "observation",
              "depression",                  "Depression",                         "medcodeid",  "observation",
              "ckd",                         "Chronic kidney disease",             "medcodeid",  "observation",
              "hypertension",                "Hypertension",                       "medcodeid",  "observation",
              "diabetes",                "Diabetes",                       "medcodeid",  "observation",
              
              # --- eFI2 frailty components: dedicated algorithm functions ---
              # (alg_* functions in frailty_functions.R handle these directly)
              "activity_limitation",         "Activity limitation",                "medcodeid",  "observation",
              "alcohol",                     "Alcohol",                            "medcodeid",  "observation",
              "atrial_fibrillation",         "Atrial fibrillation",                "medcodeid",  "observation",
              "bmi",                         "BMI",                                "medcodeid",  "observation",
              "cognitive_impairment",        "Cognitive impairment",               "medcodeid",  "observation",
              "copd",                        "COPD",                               "medcodeid",  "observation",
              "falls_for_efi2",              "Falls (eFI2 component)",             "medcodeid",  "observation",
              "housebound",                  "Housebound",                         "medcodeid",  "observation",
              "hypotension_syncope",         "Hypotension / syncope",              "medcodeid",  "observation",
              "liver_problems",              "Liver problems",                     "medcodeid",  "observation",
              "peripheral_vascular_disease", "Peripheral vascular disease",        "medcodeid",  "observation",
              "respiratory_disease",         "Respiratory disease",                "medcodeid",  "observation",
              "skin_ulcer",                  "Skin ulcer",                         "medcodeid",  "observation",
              "smoking",                     "Smoking",                            "medcodeid",  "observation",
              "social_vulnerability",        "Social vulnerability",               "medcodeid",  "observation",
              "weight_loss",                 "Weight loss",                        "medcodeid",  "observation",
              
              # --- eFI2 frailty components: routed through create_var_pre_index_medcode ---
              # (any time on or before index date)
              "cancer",                      "Cancer",                             "medcodeid",  "observation",
              "dressing_problems",           "Dressing and grooming problems",     "medcodeid",  "observation",
              "environ_problems",            "Environmental problems",             "medcodeid",  "observation",
              "heart_failure",               "Heart failure",                      "medcodeid",  "observation",
              "med_management",              "Medication management problems",     "medcodeid",  "observation",
              "mobility_problems",           "Mobility problems",                  "medcodeid",  "observation",
              "motor_neuron_disease",        "Motor neuron disease",               "medcodeid",  "observation",
              "palliative_care",             "Palliative care",                    "medcodeid",  "observation",
              "parkinsonism_tremor",         "Parkinsonism and tremor",            "medcodeid",  "observation",
              "requirement_for_care",        "Requirement for care",               "medcodeid",  "observation",
              "seizures",                    "Seizures",                           "medcodeid",  "observation",
              "stroke",                      "Stroke",                             "medcodeid",  "observation",
              "tia",                         "Transient ischaemic attack",         "medcodeid",  "observation",
              
              # --- eFI2 frailty components: routed through create_var_pre_index_medcode_5 ---
              # (within 5 years on or before index date; output named {name}_5)
              "fractures",                   "Fractures",                          "medcodeid",  "observation",
              "memory_concerns",             "Memory concerns",                    "medcodeid",  "observation",
              "peptic_ulcer_disease",        "Peptic ulcer disease",               "medcodeid",  "observation",
              "self_harm",                   "Self harm",                          "medcodeid",  "observation",
              
              # --- eFI2 fragility fracture component (separate codelist from outcome) ---
              "fragility_fractures_for_efi2","Fragility fractures (eFI2 component)","medcodeid", "observation"
      ) %>%
        mutate(
          # Path: prodcode codelists get _prodcode.txt, all others get _medcode.txt
          path  = paste0(codelists_path, name, "_", ifelse(codevar == "prodcodeid", "prodcode", "medcode"), ".txt"),
          full  = map(path, ~read_tsv(.x, col_types = cols(.default = "c"))),
          codes = map2(full, codevar, ~as.character(.x[[.y]])),
          # Extra columns used by frailty functions to split code types;
          # NA for non-frailty codelists (Option A: single unified mutate block)
          other_instructions    = map(full, ~if ("Otherinstructions"   %in% names(.x)) unique(na.omit(.x$Otherinstructions))   else NA_character_),
          time_constraint_years = map(full, ~if ("TimeConstraintyears" %in% names(.x)) unique(na.omit(.x$TimeConstraintyears)) else NA_character_),
          age_limit             = map(full, ~if ("AgeLimit"            %in% names(.x)) unique(na.omit(.x$AgeLimit))            else NA_character_)
        )
    ),
    
    # Meta-list of ICD-10 codelists (for HES linkage).
    tar_target(
      codelists_icd,
      tribble(~name,                 ~label,                ~codevar, ~extract_from,
              # Exposure
              "dementia",            "Dementia",            "icd",    "observation",
              "delirium",            "Delirium",            "icd",    "observation",
              
      ) %>%
        mutate(
          path  = paste0(codelists_path, name, "_icd10.txt"),
          full  = map(path, ~read_tsv(.x, col_types = cols(.default = "c"))),
          codes = map2(full, codevar, ~as.character(.x[[.y]]))
        )
    ),
    
    # Ethnicity codelist (TSV with additional columns eth5, eth16 needed by alg_ethnicity)
    tar_target(
      ethnicity_codelist,
      paste0(codelists_path, "ethnicity_medcode.txt"),
      format = "file"
    ),
    
    # Polypharmacy codelist (eFI2 algorithm)
    tar_target(
      polypharmacy_codelist,
      paste0(codelists_path, "polypharmacy_prodcodes.txt"),
      format = "file"
    ),
    
    # Aurum product codebrowser (March 2026) — used by alg_polypharmacy
    tar_target(
      aurum_product_codebrowser,
      aurum_product_codebrowser_path,
      format = "file"
    ),
    
    
    # ---- Raw data files ------------------------------------------------
    
    # Aurum observation parquet files 
    tar_target(
      obs_files,
      paste0(aurum_data_path, "observation/"),
      format = "file"
    ),
    
    # Aurum drug issue parquet files 
    tar_target(
      drug_files,
      paste0(aurum_data_path, "drug issue/"),
      format = "file"
    ),
    
    # HES hospitalisation file
    tar_target(
      hes_hosp,
      paste0(hes_data_path, "hes_hosp.parquet"),
      format = "file"
    ),
    
    # HES diagnosis file
    tar_target(
      hes_diag_epi,
      paste0(hes_data_path, "hes_diagnosis_epi.parquet"),
      format = "file"
    ),
    
    # Patient-level IMD (Index of Multiple Deprivation), 2019
    tar_target(
      imd_patient_data,
      read_tsv(
        paste0(imd_data_path, "patient_2019_imd_25_005905.txt"),
        col_types = cols(patid = col_character(), pracid = col_character())
      )
    ),
    
    # Practice-level IMD (used when patient-level IMD is missing)
    tar_target(
      imd_practice_data,
      read_tsv(
        paste0(imd_data_path, "practice_imd_25_005905.txt"),
        col_types = cols(pracid = col_character())
      )
    ),

    

    # ---- Study population ---------
  
    # Path to the previously-generated observation-level study population parquet
    tar_target(
      study_pop_file,
      study_pop_path,
      format = "file"
    ),
    
    # Read study population into memory
    tar_target(
      study_pop,
      open_dataset(study_pop_file) %>% collect()
    ),
    
    # ---- Study dates ----
    # Defined once here so they're consistent across all targets
    tar_target(study_start, as.Date("2010-01-01")),
    tar_target(study_end,   as.Date("2025-03-31")),
  
  #-----death data-------
  tar_target(
    name = hes_death, 
    command = open_dataset(paste0(hes_data_path, "hes_ons.parquet")) |> 
      
      # fix data formatting
      mutate(reg_date_of_death  = as.Date(reg_date_of_death , format = "%Y-%m-%d")) |> 
      
      collect ()  
    
  ), 
  
  # join death data to study pop: 
  # first remove duplicate data - multiple rows per person, sometimes with same death date, other times with different dates, keep first record when dates differ: 
  tar_target(name = hes_death_clean, 
             command = hes_death  |> 
               group_by(patid) |> 
               slice_min(order_by = reg_date_of_death, with_ties = FALSE) |> 
               ungroup(), 
             format = "parquet"
  ), 
  
  tar_target(name = ons_death_date, #add to the end
             command = study_pop |>  
               left_join(hes_death_clean, by = "patid"), 
             format = "parquet"
  ), 
  

  
#-------------------create study variables----------------------

## ---- Deprivation ----

# Combine patient- and practice-level IMD; use patient-level where available
tar_target(
  imd_all_data,
  imd_patient_data %>%
    left_join(imd_practice_data, by = "pracid") %>%
    rename(
      prac_imd = e_imd_5,
      pat_imd  = e2019_imd_5
    )
),


tar_target(
  study_pop_imd,
  study_pop %>%
    left_join(imd_all_data, by = "patid") %>%
    select(patid, indexdate, pat_imd, prac_imd) %>%
    mutate(imd = case_when(
      !is.na(pat_imd)                            ~ pat_imd,
      is.na(pat_imd) & !is.na(prac_imd)         ~ prac_imd,
      is.na(pat_imd) & is.na(prac_imd)          ~ NA
    )) %>%
    select(patid, indexdate, imd)
),


#make a list of variables defined by medcode only on or before index (to use in function)
tar_target(
  var_names_medcode_pre_index,
  c("ihd", "diabetes", "ckd", "hypertension", "depression"),
  iteration = "vector"
),
  
#----comorbidities on or before index----
  #produces a list of study_pops, each with one extra flag column for each variable in var_names_medcode_pre_index
  tar_target(
    study_pop_flagged_list,
    create_var_pre_index_medcode(
      study_pop       = study_pop,
      obs_files           = obs_files,
      codelists_aurum = codelists_aurum,
      var_name   = var_names_medcode_pre_index
    ),
    pattern = map(var_names_medcode_pre_index),
    iteration="list"

  ),
  #---- dementia variable----
  tar_target(
    study_pop_dementia,
    create_var_pre_index_medcode_hes(
      study_pop,
      obs_files,
      hes_diag_epi,
      codelists_aurum,
      codelists_icd,
      "dementia"
    )
    
  ),


  #---- ethnicity variable----
   tar_target(
    study_pop_ethnicity,
    alg_ethnicity(study_pop, 
                     obs_files,
                     ethnicity_codelist)
  
),

#----sex and age----
tar_target(
  sex_age,
  alg_sex_age(study_pop)
),

#---- primary care consultation in the year prior to index date----
tar_target(
  gp_consultation_prior_year,
  alg_gp_consultation_prior_year(
    obs_files = obs_files,
    study_pop   = study_pop,
    study_end = study_end
  )
),

#---- Delirium code on index date----
tar_target(
  delirium_index,
  alg_delirium_index(
    hes_diag_epi     = hes_diag_epi,
    study_pop       = study_pop,
    codelists_icd = codelists_icd
  )
),


# Calendar period (2-year bands from 2010 to 2025)
tar_target(
  calendar_period,
  study_pop %>%
    select(patid, indexdate) %>%
    mutate(
      index_year      = as.integer(format(indexdate, "%Y")),
      calendar_period = case_when(
        index_year %in% c(2010, 2011) ~ "2010-2011",
        index_year %in% c(2012, 2013) ~ "2012-2013",
        index_year %in% c(2014, 2015) ~ "2014-2015",
        index_year %in% c(2016, 2017) ~ "2016-2017",
        index_year %in% c(2018, 2019) ~ "2018-2019",
        index_year %in% c(2020, 2021) ~ "2020-2021",
        index_year %in% c(2022, 2023) ~ "2022-2023",
        index_year %in% c(2024, 2025) ~ "2024-2025",
        TRUE                          ~ NA_character_
      )
    ) %>%
    select(patid, indexdate, calendar_period)
),

# ---- Frailty (eFI2) — components with dedicated alg_* functions ----

tar_target(
  frailty_alg_names,
  c(
    "activity_limitation",
    "alcohol",
    "atrial_fibrillation",
    "bmi",
    "cognitive_impairment",
    "copd",
    "falls",
    "housebound",
    "hypotension_syncope",
    "liver_problems",
    "peripheral_vascular_disease",
    "respiratory_disease",
    "skin_ulcer",
    "smoking",
    "social_vulnerability",
    "weight_loss"
  ),
  iteration = "vector"
),

# Runs one branch per component on separate workers.
# Output: a named list; each element is (patid, indexdate, {component}).
tar_target(
  frailty_alg_components,
  {
    fn <- match.fun(paste0("alg_", frailty_alg_names))
    fn(
      study_pop       = study_pop,
      obs_files     = obs_files,
      codelists_aurum = codelists_aurum,
      study_end       = study_end
    )
  },
  pattern   = map(frailty_alg_names),
  iteration = "list"
),

# Convenience aliases — downstream targets reference these by name exactly
tar_target(frailty_activity_limitation,       frailty_alg_components[[1]]),
tar_target(frailty_alcohol,                   frailty_alg_components[[2]]),
tar_target(frailty_atrial_fibrillation,       frailty_alg_components[[3]]),
tar_target(frailty_bmi,                       frailty_alg_components[[4]]),
tar_target(frailty_cognitive_impairment,      frailty_alg_components[[5]]),
tar_target(frailty_copd,                      frailty_alg_components[[6]]),
tar_target(frailty_falls,                frailty_alg_components[[7]]),
tar_target(frailty_housebound,                frailty_alg_components[[8]]),
tar_target(frailty_hypotension_syncope,       frailty_alg_components[[9]]),
tar_target(frailty_liver_problems,            frailty_alg_components[[10]]),
tar_target(frailty_peripheral_vascular_disease, frailty_alg_components[[11]]),
tar_target(frailty_respiratory_disease,       frailty_alg_components[[12]]),
tar_target(frailty_skin_ulcer,                frailty_alg_components[[13]]),
tar_target(frailty_smoking,                frailty_alg_components[[14]]),
tar_target(frailty_social_vulnerability,      frailty_alg_components[[15]]),
tar_target(frailty_weight_loss,               frailty_alg_components[[16]]),



# alg_polypharmacy 
# kept as a standalone target.
tar_target(
  frailty_polypharmacy,
  alg_polypharmacy(
    study_pop                 = study_pop,
    drug_files             = drug_files,
    polypharmacy_codelist     = polypharmacy_codelist,
    study_end = study_end
  )
),

# ---- Frailty (eFI2) — components via create_var_pre_index_medcode ---------
#

tar_target(
  frailty_medcode_pre_names,
  c(
    "cancer",
    "dressing_problems",
    "environ_problems",
    "heart_failure",
    "med_management",
    "mobility_problems",
    "motor_neuron_disease",
    "palliative_care",
    "parkinsonism_tremor",
    "requirement_for_care",
    "seizures",
    "stroke",
    "tia",
    "fragility_fractures_for_efi2"
  ),
  iteration = "vector"
),

tar_target(
  frailty_medcode_pre_components,
  create_var_pre_index_medcode(
    study_pop       = study_pop,
    obs_files     = obs_files,
    codelists_aurum = codelists_aurum,
    var_name        = frailty_medcode_pre_names
  ),
  pattern   = map(frailty_medcode_pre_names),
  iteration = "list"
),

tar_target(frailty_cancer,                frailty_medcode_pre_components[[1]]),
tar_target(frailty_dressing_problems,     frailty_medcode_pre_components[[2]]),
tar_target(frailty_environ_problems,      frailty_medcode_pre_components[[3]]),
tar_target(frailty_heart_failure,         frailty_medcode_pre_components[[4]]),
tar_target(frailty_med_management,        frailty_medcode_pre_components[[5]]),
tar_target(frailty_mobility_problems,     frailty_medcode_pre_components[[6]]),
tar_target(frailty_motor_neuron_disease,  frailty_medcode_pre_components[[7]]),
tar_target(frailty_palliative_care,       frailty_medcode_pre_components[[8]]),
tar_target(frailty_parkinsonism_tremor,   frailty_medcode_pre_components[[9]]),
tar_target(frailty_requirement_for_care,  frailty_medcode_pre_components[[10]]),
tar_target(frailty_seizures,              frailty_medcode_pre_components[[11]]),
tar_target(frailty_stroke,               frailty_medcode_pre_components[[12]]),
tar_target(frailty_tia,                  frailty_medcode_pre_components[[13]]),

# Rename required: fragility_fractures_for_efi2 -> fragility_fractures_efi2
tar_target(
  frailty_fragility_fracture,
  frailty_medcode_pre_components[[14]] %>%
    rename(fragility_fracture = fragility_fractures_for_efi2)
),


# ---- Frailty (eFI2) — components via create_var_pre_index_medcode_5 -------

tar_target(
  frailty_medcode_pre5_names,
  c("fractures", "memory_concerns", "peptic_ulcer_disease", "self_harm"),
  iteration = "vector"
),

tar_target(
  frailty_medcode_pre5_components,
  create_var_pre_index_medcode_5(
    study_pop       = study_pop,
    obs_files     = obs_files,
    codelists_aurum = codelists_aurum,
    var_name        = frailty_medcode_pre5_names
  ),
  pattern   = map(frailty_medcode_pre5_names),
  iteration = "list"
),

tar_target(frailty_fractures,                frailty_medcode_pre5_components[[1]]),
tar_target(frailty_memory_concerns,     frailty_medcode_pre5_components[[2]]),
tar_target(frailty_peptic_ulcer_disease,      frailty_medcode_pre5_components[[3]]),
tar_target(frailty_self_harm,         frailty_medcode_pre5_components[[4]]),


# ---- Frailty scores ------------------------------------------------

## -- Final eFI2 score (excluding dementia) --
# Assembles all components and computes the weighted sum frailty score.
# Uses re-scaled coefficients excluding dementia
# See alg_efi2_excl_dementia in frailty_functions.R for score logic and weights.
tar_target(
  frailty_score_excl_dementia,
  alg_efi2_excl_dementia(
    activity_limitation         = frailty_activity_limitation,
    alcohol                     = frailty_alcohol,
    atrial_fibrillation         = frailty_atrial_fibrillation,
    cancer                      = frailty_cancer,
    cognitive_impairment        = frailty_cognitive_impairment,
    copd                        = frailty_copd,
    dressing_problems           = frailty_dressing_problems,
    environment_problems        = frailty_environ_problems,
    falls                       = frailty_falls,
    fracture                    = frailty_fractures,
    fragility_fracture          = frailty_fragility_fracture,
    heart_failure               = frailty_heart_failure,
    housebound                  = frailty_housebound,
    hypotension_syncope         = frailty_hypotension_syncope,
    liver_problems              = frailty_liver_problems,
    med_management              = frailty_med_management,
    memory_concerns             = frailty_memory_concerns,
    mobility_problems           = frailty_mobility_problems,
    motor_neuron_disease        = frailty_motor_neuron_disease,
    bmi                         = frailty_bmi,
    palliative_care             = frailty_palliative_care,
    parkinsonism_tremor         = frailty_parkinsonism_tremor,
    peptic_ulcer                = frailty_peptic_ulcer_disease,
    peripheral_vascular_disease = frailty_peripheral_vascular_disease,
    polypharmacy                = frailty_polypharmacy,
    requirement_for_care        = frailty_requirement_for_care,
    respiratory_disease         = frailty_respiratory_disease,
    seizures                    = frailty_seizures,
    self_harm                   = frailty_self_harm,
    skin_ulcer                  = frailty_skin_ulcer,
    smoking                     = frailty_smoking,
    social_vulnerability        = frailty_social_vulnerability,
    stroke                      = frailty_stroke,
    tia                         = frailty_tia,
    weight_loss                 = frailty_weight_loss
  )
),

#----living alone or care home----

tar_target(
  living_alone_care_home,
  alg_living_alone_care_home(obs_files,
                             study_pop,
                             codelists_aurum)
),

#----outcomes----

tar_target(
  outcomes,
  alg_outcomes(study_pop,
               hes_hosp) %>%
    left_join(ons_death_date %>% select(patid, indexdate, reg_date_of_death.x), by=c("patid","indexdate")) %>%
    rename(reg_date_of_death=reg_date_of_death.x)
)
,

  #----study pop with variables----
  #add all variables together in one dataframe 
tar_target(
  study_pop_with_variables,
  reduce(
    list(
      study_pop %>% select(patid, pracid, indexdate, observation_end, cprd_ddate, lcd, regstartdate, regenddate),
      sex_age,
      study_pop_ethnicity,
      study_pop_imd,
      study_pop_dementia,
      reduce(study_pop_flagged_list, left_join, by = c("patid", "indexdate")),
      living_alone_care_home,
      gp_consultation_prior_year,
      delirium_index,
      calendar_period,
      frailty_activity_limitation,
      frailty_alcohol,
      frailty_atrial_fibrillation,
      frailty_cancer,
      frailty_cognitive_impairment,
      frailty_copd,
      frailty_dressing_problems,
      frailty_environ_problems,
      frailty_falls,
      frailty_fractures,
      frailty_fragility_fracture,
      frailty_heart_failure,
      frailty_housebound,
      frailty_hypotension_syncope,
      frailty_liver_problems,
      frailty_med_management,
      frailty_memory_concerns,
      frailty_mobility_problems,
      frailty_motor_neuron_disease,
      frailty_bmi,
      frailty_palliative_care,
      frailty_parkinsonism_tremor,
      frailty_peptic_ulcer_disease,
      frailty_peripheral_vascular_disease,
      frailty_polypharmacy,
      frailty_requirement_for_care,
      frailty_respiratory_disease,
      frailty_seizures,
      frailty_self_harm,
      frailty_skin_ulcer,
      frailty_smoking,
      frailty_social_vulnerability,
      frailty_stroke,
      frailty_tia,
      frailty_weight_loss,
      frailty_score_excl_dementia

    ),
    left_join,
    by = c("patid", "indexdate")
  ) %>%
    left_join(ons_death_date %>% select(patid, indexdate, reg_date_of_death.x), by=c("patid","indexdate")) %>%
    rename(reg_date_of_death=reg_date_of_death.x)
) 
    
)




