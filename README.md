## Data management and analysis scripts for manuscript

Title: A multi-state model of care pathways following hospitalisation for infection in people with and without dementia in England: a cohort study

Authors: Monica Iyer, Kwabena Asare, Georgia R Gore-Langton, Anne Suffel, John Tazare, Louis Tunnicliffe, Alex Lewin, Jatinder S Minhas, Daniel Stow, Philip Williams, Kathryn E Mansfield,\* Charlotte Warren-Gash\*

\*Contributed equally

Abstract:

Background: Infections may lead to worse outcomes among people with dementia, and care pathways after infection-related hospitalisation are poorly understood, inhibiting service development.

Methods: We conducted a cohort study using primary care records from Clinical Practice Research Datalink Aurum, and linked hospital admissions and death data in England. We included people aged ≥65 years hospitalised for infection (2010-2025), followed from index admission until earliest of death, data end, or one year. The exposure was prior dementia diagnosis from primary or secondary care. Markov multi-state models described transitions between seven states (including discharge to usual residence, discharge to new care home, readmission, and death). Cox proportional hazards models estimated transition-specific hazard ratios (HRs), adjusting for key potential baseline confounders.

Results: We included 859,701 infection-related hospitalisations (680,552 unique participants; 52.2% women, median age 79 years). Compared to those without dementia, individuals with dementia had slightly lower hazards of death in hospital (HR: 0.96, 95% confidence interval [CI]: 0.94-0.98), but higher hazards of discharge to new care home (1.78 [1.72-1.84]). Among those discharged to usual residence, individuals with dementia had lower hazards of readmission (0.85 [0.84-0.86]), but higher hazards of death (1.82 [1.77-1.87]).

Conclusions: After infection-related hospitalisation, adults with dementia were more likely to be discharged to a care home and to die post-discharge, but less likely to be readmitted or die during the index admission, than those without dementia. For people with dementia, infection-related admissions may signal a need for better anticipatory care, discharge planning, and longer-term support.

## Repository scope

This repository includes data management code to create variables and all analysis code used in the study. 

Raw data are not provided per CPRD regulations. There are placeholders for file paths to data.

1.  The Create folder contains the script for creating the study eligibility criteria for the main analysis and a sensitivity analysis. It has placeholders where raw CPRD and linked data should be to create the study population.
2.  The Variables folder contains code for creating the study variables using the targets package. The \_targets file runs using the scripts functions.R and frailty_functions.R. An additional script named global_paths with path names of where the data and codelists are kept is required but not provided for privacy.
3.  The Analysis folder contains scripts for the descriptive analysis, main multi-state model analysis, and sensitivity analyses. These are not in targets format.

## Protocol

A link to the protocol can be found: <https://www.cprd.com/approved-studies/care-pathways-following-hospitalisation-infection-people-and-without-dementia>

## Codelists
Codelists used in this project can be downloaded from: <https://datacompass.lshtm.ac.uk/id/eprint/5417/>
