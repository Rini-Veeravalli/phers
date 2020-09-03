# Residual PheRS #


# Load dependencies -----------------------------------

library(tidyverse)
library(data.table)



# Load data and mappings ------------------------------

map_disease_omim = read.csv("~/re_gecip/health_records/PheRS/code/omimdb.csv")

map_phe_icd10 = read_delim("~/re_gecip/health_records/PheRS/ref_data/icd10_phe_mapping_cleaned.txt", col_types = "cc", delim = "\t") %>%
  distinct() %>%
  select("ICD10" = 1, PHECODE) %>%
  mutate(PhecodeCHR = as.character(as.numeric(PHECODE))) 
### PhecodeCHR preserves original phecode format i.e. removes .0


icds = read.delim("~/re_gecip/health_records/PheRS/extract_ICD_codes/hes_apc_clean_icds_v7.txt", na.strings = "NA", stringsAsFactors = F)

disease_label = fread("~/Documents/rare_diseases_participant_dise_2019-09-26_16-44-06.csv") %>% 
  select(participant_id, starts_with("normalised")) 

participants = fread("~/Documents/participant_2019-10-24_16-53-25.csv") %>% 
  select(participant_id, programme, year_of_birth, participant_phenotypic_sex, participant_type)

icds_admidate_pre2014 = read.csv("participants_admidate_pre2014_updated.txt") %>% 
  select(participant_id) %>% 
  distinct()



# Select population -----------------------------------

# ICD10 data of Total participants with EHR data and admission date of > 5 years

population_icds = icds %>%
  select(participant_id, ICD10_CLEAN1) %>%
  filter(participant_id %in% icds_admidate_pre2014$participant_id) %>%
  left_join(map_phe_icd10, by = c("ICD10_CLEAN1" = "ICD10")) 

# add record length (unique number of years of hes records)

record_length = fread("~/re_gecip/health_records/PheRS/extract_ICD_codes/hes_apc_icds_v7_long_format.csv") %>%
  select(participant_id, admidate) %>%
  filter(participant_id %in% unique(population_icds$participant_id)) %>%
  distinct() %>%
  mutate(admidate = format(as.Date(admidate, format="%Y-%m-%d"), "%Y")) %>%
  distinct() %>%
  group_by(participant_id) %>%
  summarise(record_len = n_distinct(admidate))

# Participant data of Total participants with EHR data and admission date of > 5 years and disease data

population_participants = population_icds %>%
  select(participant_id) %>%
  distinct() %>%
  left_join(participants, by = "participant_id") %>% 
  mutate(age = 2019 - year_of_birth,
         binaryGender = ifelse(participant_phenotypic_sex == "Male", 1, 0)) %>% 
  left_join(disease_label) %>%
  filter(!is.na(normalised_specific_disease)) %>%
  left_join(record_length)



# Create functions ------------------------------------

# create_filename()           # get disease name to use for naming files
# calculate_residual_phers()  # calculate rphers and create table of rphers results
# evaluate_residual_phers()   # evaluate rphers per participant across diseases


create_filename <- function(disease_index) {
  disease_filename <- gsub(" ", "_", map_disease_omim$NormalisedSpecificDisease[disease_index])
  return(disease_filename)
}

calculate_residual_phers <- function(disease_filename) {
  # calculate z-score to compare PheRS across different diseases 
  # (z score - represent rphers by its deviation from the mean by sd)
  # by how many standard deviations does the rPheRS differ from the mean?
  
  # PRS in population_participants
  PRS = read.delim(paste0("phers_test/", disease_filename, "_phers.txt"))
  
  # same as population_participants, but without joining for rare disease data or recordlen
  rphers_participants = population_icds %>%
    dplyr::select(participant_id) %>%
    distinct() %>%
    left_join(participants, by = "participant_id") %>% 
    mutate(age = 2019 - year_of_birth,
           binaryGender = ifelse(participant_phenotypic_sex == "Male", 1, 0))
  
  
  phers_participants <- PRS[PRS$participant_id %in% rphers_participants$participant_id, ]
  participantPheRS = phers_participants %>%
    inner_join(rphers_participants, by = "participant_id") %>%
    dplyr::select(participant_id, age, binaryGender, participant_type, PRS) %>%
    distinct()
  participantPheRS$rPheRS <- 0
  
  # fit linear model 
  linearmodel <- lm(PRS ~ age + binaryGender, data = participantPheRS)
  #print(linearmodel)
  #print(summary(linearmodel))
  
  coeffs = coefficients(linearmodel); coeffs

  # for each participant in PRS, get phers
  # then calculate EPheRS
  for (row in 1:nrow(participantPheRS)) {

    i_age = participantPheRS$age[row]
    i_gender = participantPheRS$binaryGender[row]
    Ephers = coeffs[1] + coeffs[2]*i_age + coeffs[3]*i_gender

    # calculate rPheRS = PheRS - EPheRS
    rPheRS = (participantPheRS$PRS[row] - Ephers)
 
    # add rPheRS column to sampHers
    participantPheRS$rPheRS[row] <- rPheRS
    
  }
  
  # z-score of rPheRS
  participantPheRS$zrPheRS <- scale(participantPheRS$rPheRS)
 
  # add disease name column to sampPheRS
  participantPheRS$disease <- disease_filename
  
  participantPheRS %>%
    write.csv(paste0("residual_test/", disease_filename,"_participant_residuals.csv"))
  
  return(participantPheRS)
}


evaluate_residual_phers <- function() {

  list_of_files <- list.files(path="residual_test", pattern = ".*.csv", full.names = TRUE) 
  residual_part_phers <- rbindlist(sapply(list_of_files, fread, simplify = F),
                                   use.names = T, idcol = "FileName", fill = T)
  
  # for each participant, get zrPheRS for all diseases (descending order) 
  
  for (unique_participant in unique(residual_part_phers$participant_id)){
    zrPheRS_profile = (residual_part_phers[residual_part_phers$participant_id == unique_participant,]) 
    zrPheRS_profile =  zrPheRS_profile[order(zrPheRS, decreasing = TRUE),]
    print(zrPheRS_profile)
    hist(zrPheRS_profile$zrPheRS)
    pid = unique_participant 
    break
  }

  print(population_participants[population_participants$participant_id == pid,])
  
  print(disease_label[disease_label$participant_id == pid,])

}



# main script

for (row in 1:nrow(map_disease_omim)) {
  tryCatch({
    disease_filename = create_filename(map_disease_omim$NormalisedSpecificDisease[row])
    participantPheRS = calculate_residual_phers(disease_filename)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
# for each participant, across all diseases
evaluate_residual_phers()
