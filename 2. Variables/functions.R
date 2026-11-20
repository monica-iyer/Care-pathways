
##Monica Iyer
#April 22, 2026
#Care pathways study
#Functions to create variables
#These functions are referenced in _targets


#---------Create medcode variable before index--------------------
#this function creates a binary variable for events in aurum occurring on or before index date

create_var_pre_index_medcode <- function(study_pop,       # study population dataframe
                                         obs_files,        # path(s) to Aurum observation parquet files
                                         codelists_aurum,  # codelist meta-list (from _targets.R)
                                         var_name         # variable name string (must match codelist name)
) {
  
  # Extract medcodeids for this variable from the codelist meta-list
  medcodes <- codelists_aurum %>%
    filter(name == var_name) %>%
    pull(codes) %>%
    .[[1]]      # character vector of medcodeids
  

  # Filter Aurum observation data to relevant medcodes, study population, and study period
  # with codelist vectors, which are always character
  events <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% medcodes) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%  
    inner_join(study_pop %>% select(patid, indexdate), by="patid") %>%
    collect() %>%
    filter(obsdate <=indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(!!var_name := 1L)
  
  
  # Join flags back to study population; patients without a qualifying event get 0
  study_pop %>%
    left_join(events, by = c("patid", "indexdate")) %>%
    mutate(!!var_name := tidyr::replace_na(.data[[var_name]], 0L)) %>%
    select(patid, indexdate, !!var_name)
  
}

#---------create variable using medcode on or before index, looking just 5 years before index---

create_var_pre_index_medcode_5 <- function(study_pop,       # study population dataframe
                                           obs_files,        # path(s) to Aurum observation parquet files
                                           codelists_aurum,  # codelist meta-list
                                           var_name          # variable name string
) {
  
  # Extract medcodeids for this variable
  medcodes <- codelists_aurum %>%
    filter(name == var_name) %>%
    pull(codes) %>%
    .[[1]]
  
  events <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% medcodes) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%  
    inner_join(study_pop %>% select(patid, indexdate), by="patid") %>%
    collect() %>%
    filter(obsdate <= indexdate & indexdate - obsdate <= 1826.25) %>% #flag patients with an event within 5 years before index
    distinct(patid, indexdate) %>%
    mutate(!!var_name := 1L)

  
  # Join flags back to study population; patients without a qualifying event get 0
  study_pop %>%
    left_join(events, by = c("patid", "indexdate")) %>%
    mutate(!!var_name := tidyr::replace_na(.data[[var_name]], 0L)) %>%
    select(patid, indexdate, !!var_name)
  
}
#---------Create medcode or HES variable before index--------------------

#this function creates a binary variable for events in aurum or HES occurring on or before index date

create_var_pre_index_medcode_hes <- function(study_pop,        # study population dataframe
                                             obs_files,         # path(s) to Aurum observation parquet files
                                             hes_diag_epi,          # path(s) to HES diagnosis parquet files
                                             codelists_aurum,   # medcode codelist meta-list
                                             codelists_icd,     # ICD-10 codelist meta-list
                                             var_name           # variable name string

                                             ) {
  
  # Extract medcodeids for this variable
  medcodes <- codelists_aurum %>%
    filter(name == var_name) %>%
    pull(codes) %>%
    .[[1]]
  
  # Extract ICD-10 codes for this variable
  icdcodes <- codelists_icd %>%
    filter(name == var_name) %>%
    pull(codes) %>%
    .[[1]]
  
  patids_in_study <- unique(study_pop$patid)
  
  # Aurum: filter to relevant medcodes and study population
  aurum_events <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    filter(medcodeid %in% medcodes) %>%
    filter(patid %in% patids_in_study) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%    # Aurum date format
    collect() %>%
    select(patid, event_date = obsdate)
  
  # HES: filter to relevant ICD-10 codes and study population
  hes_events <- open_dataset(hes_diag_epi) %>%
    select(patid, epistart, ICD) %>%
    filter(ICD %in% icdcodes) %>%
    filter(patid %in% patids_in_study) %>%
    mutate(epistart = as.Date(epistart, "%Y/%m/%d")) %>%  # HES date format
    collect() %>%
    select(patid, event_date = epistart)
  
  # Stack Aurum and HES events
  all_events <- bind_rows(aurum_events, hes_events)
  
  # Flag events that occurred on or before the patient's index date
  event_before_index <- all_events %>%
    left_join(study_pop %>% select(patid, indexdate),
              by = "patid") %>%
    filter(event_date <= indexdate) %>%
    distinct(patid, indexdate) %>%
    mutate(!!var_name := 1L)
  
  # Join flags back to study population; patients without a qualifying event get 0
  study_pop %>%
    left_join(event_before_index, by = c("patid", "indexdate")) %>%
    mutate(!!var_name := tidyr::replace_na(.data[[var_name]], 0L)) %>%
    select(patid, indexdate, !!var_name)
  
}




#-------Ethnicity; algorithm previously published-------

alg_ethnicity <- function(study_pop,
                          obs_files,         
                          ethnicity_codelist   # path to ethnicity codelist file (TSV)
) {
  
  ethnicity_codelist <- read_tsv(ethnicity_codelist,
                                 col_types = cols(medcodeid = col_character()))
  ethcodes <- as.character(ethnicity_codelist$medcodeid)
  
  ethn_recs <-       open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    mutate(medcodeid = as.character(medcodeid)) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(medcodeid %in% ethcodes) %>%
    collect() %>%
    left_join(ethnicity_codelist[c("medcodeid", "eth5", "eth16")], by = "medcodeid")
  
  df <- as.data.table(ethn_recs)
  
  #tagging patients with the same ethnicity entry within the same year
  df<- df[order(patid)]
  df <- df[, sysyear := year(obsdate)]
  df[, duplicates := .N, by = c("patid", "sysyear", "medcodeid")]
  #keeping the first ethnicity entry per year if there is a duplicate
  df[duplicates > 1, first_entry := lapply(.SD, min), by = c("patid", "sysyear"),
     .SDcols = "obsdate"]
  df[duplicates >1, keep_row := obsdate == first_entry]
  df[duplicates == 1, keep_row := TRUE]
  df <- df[keep_row == TRUE]
  
  df[, duplicates := NULL]
  df[, first_entry := NULL]
  df[, keep_row := NULL]
  
  
  #count the observations of ethnicity by patient
  df <- df[, by = c("patid", "eth5"), eth5_count := .N]
  
  
  
  
  ###adding up ethnicity - 5 categories and 16 categories
  df[eth5 == "0. White", eth5whitecount := .N, by = "patid"]
  df[is.na(eth5whitecount), eth5whitecount := 0]
  df[eth5 == "1. South Asian", eth5sacount := .N, by = "patid"]
  df[is.na(eth5sacount), eth5sacount := 0]
  df[eth5 == "2. Black", eth5blackcount := .N, by = "patid"]
  df[is.na(eth5blackcount), eth5blackcount := 0]
  df[eth5 == "3. Other", eth5othercount := .N, by = "patid"]
  df[is.na(eth5othercount), eth5othercount := 0]
  df[eth5 == "4. Mixed", eth5mixedcount := .N, by = "patid"]
  df[is.na(eth5mixedcount), eth5mixedcount := 0]
  df[eth5 == "5. Not Stated", notstated5count := .N, by = "patid"]
  df[is.na(notstated5count), notstated5count := 0]
  
  
  df[eth16 =="1. British", british16count := .N, by = "patid"]
  df[is.na(british16count), british16count := 0]
  df[eth16 =="2. Irish", irish16count := .N, by = "patid"]
  df[is.na(irish16count), irish16count := 0]
  df[eth16 =="3. Other White", otherwhite16count := .N, by = "patid"]
  df[is.na(otherwhite16count), otherwhite16count := 0]
  df[eth16 =="4. White and Black Caribbean", whitecarib16count := .N, by = "patid"]
  df[is.na(whitecarib16count), whitecarib16count := 0]
  df[eth16 =="5. White and Black African", whiteaf16count := .N, by = "patid"]
  df[is.na(whiteaf16count), whiteaf16count := 0]
  df[eth16 =="6. White and Asian", whiteasian16count := .N, by = "patid"]
  df[is.na(whiteasian16count), whiteasian16count := 0]
  df[eth16 =="7. Other Mixed", othermixed16count := .N, by = "patid"]
  df[is.na(othermixed16count), othermixed16count := 0]
  df[eth16 =="8. Indian", indian16count := .N, by = "patid"]
  df[is.na(indian16count), indian16count := 0]
  df[eth16 =="9. Pakistani", pak16count := .N, by = "patid"]
  df[is.na(pak16count), pak16count := 0]
  df[eth16 =="10. Bangladeshi", bangla16count := .N, by = "patid"]
  df[is.na(bangla16count), bangla16count := 0]
  df[eth16 =="11. Other Asian", otherasian16count := .N, by = "patid"]
  df[is.na(otherasian16count), otherasian16count := 0]
  df[eth16 =="12. Caribbean", carib16count := .N, by = "patid"]
  df[is.na(carib16count), carib16count := 0]
  df[eth16 =="13. African", african16count := .N, by = "patid"]
  df[is.na(african16count), african16count := 0]
  df[eth16 =="14. Other Black", otherblack16count := .N, by = "patid"]
  df[is.na(otherblack16count), otherblack16count := 0]
  df[eth16 =="15. Chinese", chinese16count := .N, by = "patid"]
  df[is.na(chinese16count), chinese16count := 0]
  df[eth16 =="16. Other ethnic group", other16count := .N, by = "patid"]
  df[is.na(other16count), other16count := 0]
  df[eth16 =="17. Not Stated", notstated16count := .N, by = "patid"]
  df[is.na(notstated16count), notstated16count := 0]
  
  
  
  #flagging those patients who only have non stated ethnicity
  df<- df[,  nonstatedonly := (notstated5count > 0 &
                                 eth5whitecount == 0 &
                                 eth5sacount == 0 &
                                 eth5blackcount == 0 &
                                 eth5othercount == 0 &
                                 eth5mixedcount == 0)]
  
  
  ###finding the most common ethnicity, excluding ethnicity not stated
  # Get Row wise max in R
  df <- df[,eth5max_f :=  names(.SD)[max.col(.SD, ties.method = "first")], .SDcols = 8:12] #the eth5 columns
  df <- df[,eth5max_l :=  names(.SD)[max.col(.SD, ties.method = "last")], .SDcols = 8:12]
  df <- df[,eth16max_f :=  names(.SD)[max.col(.SD, ties.method = "first")], .SDcols = 14:29]
  df <- df[,eth16max_l :=  names(.SD)[max.col(.SD, ties.method = "last")], .SDcols = 14:29]
  
  #get the rows with equally common ethnicities, excluding non stated
  df<- df[, equaleth5 := eth5max_f != eth5max_l & notstated5count == 0 ]
  df<- df[, equaleth16 := eth16max_f != eth16max_l & notstated16count  ==0]
  
  
  
  ###---FINAL STEP - assigning ethnicity to patid
  #case 1: ethnicity is not stated - assign not stated
  eth <- df[nonstatedonly==T, ethnicity_5 := 5]
  eth <- df[nonstatedonly==T, ethnicity_16 := 17]
  
  
  
  ####case 2: assign the most common ethnicity
  eth[equaleth5 == F & nonstatedonly == F, ethnicity_5c := eth5max_f]
  eth[equaleth16 == F & nonstatedonly == F, ethnicity_16c := eth16max_f]
  
  #intermediate step, translating the character back into ethnicity numbers
  eth[ethnicity_5c == "eth5whitecount", ethnicity_5 := 0]
  eth[ethnicity_5c == "eth5sacount", ethnicity_5 := 1]
  eth[ethnicity_5c == "eth5blackcount", ethnicity_5 := 2]
  eth[ethnicity_5c == "eth5othercount", ethnicity_5 := 3]
  eth[ethnicity_5c == "eth5mixedcount", ethnicity_5 := 4]
  
  eth[ethnicity_16c ==  "british16count", ethnicity_16 := 1]
  eth[ethnicity_16c ==  "irish16count", ethnicity_16 := 2]
  eth[ethnicity_16c ==  "otherwhite16count", ethnicity_16 := 3]
  eth[ethnicity_16c ==  "whitecarib16count", ethnicity_16 := 4]
  eth[ethnicity_16c ==  "whiteaf16count", ethnicity_16 := 5]
  eth[ethnicity_16c ==  "whiteasian16count", ethnicity_16 := 6]
  eth[ethnicity_16c ==  "othermixed16count", ethnicity_16 := 7]
  eth[ethnicity_16c ==  "indian16count", ethnicity_16 := 8]
  eth[ethnicity_16c ==  "pak16count", ethnicity_16 := 9]
  eth[ethnicity_16c ==  "bangla16count", ethnicity_16 := 10]
  eth[ethnicity_16c ==  "otherasian16count", ethnicity_16 := 11]
  eth[ethnicity_16c ==  "carib16count", ethnicity_16 := 12]
  eth[ethnicity_16c ==  "african16count", ethnicity_16 := 13]
  eth[ethnicity_16c ==  "otherblack16count", ethnicity_16 := 14]
  eth[ethnicity_16c ==  "chinese16count", ethnicity_16 := 15]
  eth[ethnicity_16c ==  "other16count", ethnicity_16 := 16]
  
  
  
  ###case 3: if equal common
  #finding the latest date of ethnicity record and using this record as eth5/16 code
  eth[equaleth5 ==T, latestethdate := lapply(.SD, max), by = "patid", .SDcols = "obsdate"]
  eth[equaleth5 == T, eth_pick := obsdate == latestethdate]
  eth[equaleth5 == T & eth_pick == T, ethnicity_5 := eth5]
  
  eth[equaleth16 ==T, latestethdate := lapply(.SD, max), by = "patid", .SDcols = "obsdate"]
  eth[equaleth16 == T, eth_pick := obsdate == latestethdate]
  eth[equaleth16 == T & eth_pick == T, ethnicity_16 := eth16]
  
  
  ###--- creating a clean ethnicity file
  eth <- unique(eth[, list(patid, ethnicity_5, ethnicity_16)])
  
  
  ###case 4: different ethnicities assigned on the same day
  # assign to "other" ethnicitiy group
  eth[duplicated(patid) | duplicated(patid, fromLast=TRUE), ethnicity_16 := 16]
  eth[duplicated(patid) | duplicated(patid, fromLast=TRUE), ethnicity_5 := 3]
  
  eth <- unique(eth[, list(patid, ethnicity_5, ethnicity_16)])
  # Join back to study population and rename numeric to character
  study_pop %>%
    left_join(eth, by = "patid") %>%
    select(patid, indexdate, ethnicity_5, ethnicity_16) %>%
    mutate(
      ethnicity_5 = case_when(
        ethnicity_5 == 0 ~ "0. White",
        ethnicity_5 == 1 ~ "1. South Asian",
        ethnicity_5 == 2 ~ "2. Black",
        ethnicity_5 == 3 ~ "3. Other",
        ethnicity_5 == 4 ~ "4. Mixed",
        ethnicity_5 == 5 ~ "5. Not Stated"
      ),
      ethnicity_16 = case_when(
        ethnicity_16 == 1  ~ "1. British (White)",
        ethnicity_16 == 2  ~ "2. Irish (White)",
        ethnicity_16 == 3  ~ "3. Other White",
        ethnicity_16 == 4  ~ "4. White and Black Caribbean",
        ethnicity_16 == 5  ~ "5. White and Black African",
        ethnicity_16 == 6  ~ "6. White and Asian",
        ethnicity_16 == 7  ~ "7. Other Mixed",
        ethnicity_16 == 8  ~ "8. Indian",
        ethnicity_16 == 9  ~ "9. Pakistani",
        ethnicity_16 == 10 ~ "10. Bangladeshi",
        ethnicity_16 == 11 ~ "11. Other Asian",
        ethnicity_16 == 12 ~ "12. Black Caribbean",
        ethnicity_16 == 13 ~ "13. Black African",
        ethnicity_16 == 14 ~ "14. Other Black",
        ethnicity_16 == 15 ~ "15. Chinese",
        ethnicity_16 == 16 ~ "16. Other ethnic group",
        ethnicity_16 == 17 ~ "17. Not Stated"
      )
    )
  
}

#----Outcomes-----------------------------------------
#create hospitalisations and spells

alg_outcomes <- function(study_pop,
                                hes_hosp) {
  
#get hospitalisations for each patid that are on or after index
  open_dataset(hes_hosp) %>%
    select(patid, admidate, discharged, admisorc, disdest) %>%
    mutate(admidate=as.Date(admidate, "%Y/%m/%d"),
           discharged=as.Date(discharged, "%Y/%m/%d")) %>%
    inner_join(study_pop %>% select(patid, indexdate), by="patid") %>%
    filter(admidate >= indexdate) %>%
    filter(!is.na(discharged) &
             discharged>=admidate) %>% #only include admissions where we have a discharge date and it's after admidate
    arrange(patid, indexdate, admidate) %>%
    group_by(patid, indexdate) %>%
    collect() %>%
    mutate(
      # Check if this admission is a continuation (gap from prev discharge <= 2 days)
      prev_discharged = lag(discharged),
      new_spell = is.na(prev_discharged) | as.numeric(admidate - prev_discharged) > 2, #spell = hospitalisations with 2 day allowable gap
      spell_n = cumsum(new_spell),
      #Use NHS definitions for disdest== care home
      #Per protocol, any discharge that's not to care home or death -> home
      disdest=case_when(disdest>=54 & disdest<=65 ~ "Care home",
                        disdest==69 ~ "Care home",
                        disdest==79 ~ "Died",
                        disdest>=85 & disdest<=86~ "Care home",
                        TRUE ~ "Home"),
      admisorc=case_when(admisorc>=54 & admisorc<=65 ~ "Care home",
                         admisorc==69 ~ "Care home",
                         admisorc==79 ~ "Died",
                         admisorc>=85 & admisorc<=86~ "Care home",
                         TRUE ~ "Home")) %>%
    group_by(patid, indexdate, spell_n) %>% #spell_n= which spell number per cohort entry
    summarise(
      spell_start = min(admidate), #for each spell, take the min admidate as spell start date and last discharge date
      spell_end   = max(discharged),
      disdest=last(disdest),#discharge destination is the last disdest in the spell
      admisorc=first(admisorc), #admission source is first admisorc in the spell
       .groups = "drop"
    ) %>%
    # keep only admissions within a year of index date
    filter(as.numeric(spell_start-indexdate)<=365.25) %>%
    mutate(failed_discharge=ifelse(as.numeric(spell_start -indexdate)>2 & as.numeric(spell_start-indexdate)<=30, 1, 0),
           readmission=ifelse(as.numeric(spell_start-indexdate)>2 & as.numeric(spell_start-indexdate)<365.25, 1, 0),
           los=as.numeric(spell_end-spell_start)) %>%
    mutate(spell_end = if_else(spell_end < spell_start, as.Date(NA), spell_end)) %>%
    select(patid, indexdate, spell_n, spell_start, spell_end, admisorc, disdest, los, failed_discharge, readmission)
  

  
}

#---- Living alone or care home----

alg_living_alone_care_home <- function(obs_files,
                                       study_pop,
                                       codelists_aurum) {
  
  patids_in_study <- unique(study_pop$patid)
  
  # Extract each codelist separately
  living_alone_codes <- codelists_aurum %>%
    filter(name == "living_alone") %>%
    pull(codes) %>%
    .[[1]]
  
  care_home_codes <- codelists_aurum %>%
    filter(name == "care_home") %>%
    pull(codes) %>%
    .[[1]]
  
  all_codes <- c(living_alone_codes, care_home_codes)
  
  # Read and filter obs data
  obs_data <- open_dataset(obs_files) %>%
    select(patid, obsdate, medcodeid) %>%
    filter(medcodeid %in% all_codes) %>%
    filter(patid %in% patids_in_study) %>%
    mutate(obsdate = as.Date(obsdate, "%d/%m/%Y")) %>%
    collect()
  
  # Join to study_pop and apply index date filter
  obs_data <- study_pop %>%
    select(patid, indexdate) %>%
    left_join(obs_data, by = "patid") %>%
    filter(obsdate < indexdate) %>%                     # any time before index date
    mutate(
      code_type = case_when(
        medcodeid %in% living_alone_codes ~ "living_alone",
        medcodeid %in% care_home_codes    ~ "care_home"
      )
    ) %>%
    arrange(patid, indexdate, obsdate, desc(code_type == "care_home")) %>%  # care_home = TRUE sorts first (TRUE > FALSE)
    group_by(patid, indexdate) %>%
    slice_max(obsdate, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Join back to study_pop to preserve all patid/indexdate combos
  study_pop %>%
    select(patid, indexdate) %>%
    left_join(obs_data %>% select(patid, indexdate, code_type), 
              by = c("patid", "indexdate")) %>%
    mutate(
      living_alone = as.integer(code_type == "living_alone"),
      care_home    = as.integer(code_type == "care_home")
    ) %>%
    mutate(across(c(living_alone, care_home), ~replace_na(., 0))) %>%
    select(patid, indexdate, living_alone, care_home)
}


#----PRIMARY CARE CONSULTATION IN YEAR PRIOR TO INDEX----

# Flags whether a patient had at least one primary care consultation (any Aurum observation record) in the 365 days before their index date.
# Output: gp_consultation_prior_year = 1 if at least one record, 0 otherwise



alg_gp_consultation_prior_year <- function(obs_files, 
                                           study_pop,   # study start date (Date)
                                           study_end     # study end date (Date)
) {
  
  
  
  temp <- open_dataset(obs_files) |> 
    select(patid, obsdate) |> 
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |> 
    inner_join(study_pop[c("patid","indexdate")], by="patid") |>
    mutate(indexdate = as.Date(indexdate, "%d/%m/%Y")) %>%
    filter(obsdate<indexdate) |> 
    filter(as.integer(indexdate)-as.integer(obsdate)<365.25) |> 
    distinct(patid, indexdate) %>%
    mutate(gp_consultation_prior_year = 1L) %>%
    collect()
  
  study_pop %>%
    left_join(temp, by = c("patid", "indexdate")) %>%
    mutate(gp_consultation_prior_year = replace_na(gp_consultation_prior_year, 0L)) %>%
    select(patid, indexdate, gp_consultation_prior_year)
}


#----DELIRIUM ON INDEX DATE----

# Flags whether a patient had a delirium code recorded within seven days of index
# Output: delirium_index = 1 if any delirium code on index date, 0 otherwise


alg_delirium_index <- function(hes_diag_epi,      
                               study_pop,         # study population dataframe
                               codelists_icd          # study end date (Date)
) {
  
  icdcodes <- codelists_icd %>%
    filter(name == "delirium") %>%
    pull(codes) %>%
    .[[1]]
  
  temp <- open_dataset(hes_diag_epi) %>%
    select(patid, epistart, ICD) %>%
    filter(ICD %in% icdcodes) %>%
    mutate(epistart = as.Date(epistart, "%Y/%m/%d")) %>%  
    inner_join(study_pop %>% select(patid, indexdate), by="patid") %>%
    collect() %>%
    filter(epistart>=indexdate & as.numeric(epistart-indexdate)<8) %>%
    distinct(patid, indexdate) %>%
    mutate(delirium_index=1L) %>%
    collect()
  
  study_pop %>%
    left_join(temp, by = c("patid", "indexdate")) %>%
    mutate(delirium_index = replace_na(delirium_index, 0L)) %>%
    select(patid, indexdate, delirium_index)
  
}


#----sex and age----
alg_sex_age <- function(study_pop) {
  
  result <- study_pop %>%
    select(patid, indexdate, gender, yob) %>%
    mutate(
      sex = case_when(
        gender == 1L ~ "Male",
        gender == 2L ~ "Female",
        gender == 3L ~ "Indeterminate",
        TRUE         ~ NA_character_
      ),
      dob_mid_year = as.Date(paste0(yob, "-07-01")),
      age          = as.numeric(floor(
        as.numeric(indexdate - dob_mid_year) / 365.25
      ))
    ) %>%
    select(patid, indexdate, sex, age)
  
  unmapped <- sum(is.na(result$sex) & !is.na(study_pop$gender))
  if (unmapped > 0) {
    warning("alg_sex_age: ", unmapped,
            " record(s) had unrecognised gender codes and were set to NA")
  }
  
  result
  
}


