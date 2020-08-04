library(tidyverse)
library(data.table)

#### LOAD DATA

# load rare disease data
rare_diseases = read.csv("~/Documents/rare_diseases_participant_dise_2019-09-26_16-44-06.csv")
# SM using this dataset, I'm using slightly earlier one
# disease_label = fread("~/Documents/v7_data/rare_diseases_participant_dise_2019-10-26_19-07-08.tsv") %>% 
disease_label = fread("~/Documents/rare_diseases_participant_dise_2019-09-26_16-44-06.csv") %>% 
  select(participant_id, normalised_specific_disease) 

# SM using this dataset, I'm using slightly earlier one
# participants = read.delim("~/Documents/v7_data/participant_2019-10-26_18-32-36.tsv") %>% 
participants = fread("~/Documents/participant_2019-10-24_16-53-25.csv") %>% 
  select(participant_id, programme, year_of_birth, participant_phenotypic_sex, participant_type)


# --- load potential control data, controlled for having more than 5 years of follow up
rare_disease_EHR_controlled = read.csv("~/re_gecip/health_records/PheRS/participants_admidate_pre2014.txt") %>% 
  select(participant_id) %>% 
  distinct()

# --- read in all rare participants with EHR data
participantsWithEHR = fread("~/re_gecip/health_records/PheRS/extract_ICD_codes/hes_apc_clean_icds_v7.txt") %>% 
  select(participant_id) %>% 
  distinct()




####################################


df <- data.frame(unique_norm_spec_dis)

templist <- list()

for (i in unique_norm_spec_dis){
  print(i)
  countnum <- tally(unique((disease_label[(disease_label$normalised_specific_disease == i)])))
  print(countnum)
  
  templist <- c(templist, countnum)
  
}

#add tmeplist to df
df[, 'num_recruited_cases'] <- unlist(templist)

tally(df[df$num_recruited_cases < 10, ])
## 44 diseases with recruited cases < 10


##### count how many participants per disease have EHR  AND  have admission data of over 5 years #####

create_case_sets <- function(disease_name){
  ### GET DISEASE OF INTEREST AND CASES
  ### UN-hardcode this
  disease_interest = disease_label %>% 
    # filter(normalised_specific_disease == "Familial Hypercholesterolaemia") %>% 
    filter(normalised_specific_disease == disease_name) %>% 
    select(participant_id, normalised_specific_disease)
  
  #print(disease_interest)
  
  # cases of the disease of interest, with added info on age and gender
  cases = disease_interest %>% 
    inner_join(participants) %>% 
    mutate(age = 2019 - year_of_birth,
           binaryGender = ifelse(participant_phenotypic_sex == "Male", 0, 1)) ### TO DO - CHECK
  
  #print(cases)
  
  # check population demographics data e.g. mean age, gender and number of the cases for the disease of interest
  demographics_cases = cases %>% 
    summarise(Avg_age = mean(age),
              PercFem = mean(binaryGender),
              Count = n())
  
  # out of all participants with EHR disease for rare diseases and admidate>5 years, create a table with age, gender and disease label
  population = rare_disease_EHR_controlled %>% 
    inner_join(participants) %>% 
    mutate(age = 2019 - year_of_birth,
           binaryGender = ifelse(participant_phenotypic_sex == "Male", 0, 1)) %>% ### TO DO - CHECK
    inner_join(disease_label)
  
  #print(population)
  
  # remove duplicate participants in disease of interest cases, and find which have EHR data
  cases.cur = cases %>%
    filter(!duplicated(participant_id)) %>%
    mutate(Pheno = "Case") %>%
    filter(participant_id %in% rare_disease_EHR_controlled$participant_id)

  #print(cases.cur)  
  selected_case_num <- length(unique(cases.cur$participant_id))
  print(selected_case_num)
  
  # any(duplicated(cases.cur$participant_id))
  return(selected_case_num)
  
}

num_selected_cases <- list()

for (i in unique_norm_spec_dis){
  selected_cases <- create_case_sets(i)
  num_selected_cases <- c(num_selected_cases, selected_cases)
  #break
}

# add num_selected_cases list to df
df[, 'num_selected_cases'] <- unlist(num_selected_cases)

tally(df[df$num_selected_cases < 10, ])
## 44 diseases with recruited cases < 10

df %>%
  write.csv(paste0("Num_Cases_Normalised_Specific_Disease.csv"))
