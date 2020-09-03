# code snippet to investage cases with PheRS = 0 for 3 chosen diseases
# can be run in R 3.5.0 after running phers.R and loading in and creating the appropriate datasets and variables


# number of unique admissions in HES APC
N_admissions = fread("~/re_gecip/health_records/PheRS/extract_ICD_codes/hes_apc_icds_v7_long_format.csv") %>%
  dplyr::select(participant_id, admidate) %>%
  filter(participant_id %in% unique(population_icds$participant_id)) %>%
  distinct() %>%
  mutate(admidate = as.Date(admidate, format="%Y-%m-%d")) %>%
  distinct() %>%
  group_by(participant_id) %>%
  summarise(N_admissions = n_distinct(admidate))

# number of unique ICD10 codes in HES APC
N_diagnoses = icds %>%
  dplyr::select(participant_id, ICD10_CLEAN1) %>%
  filter(participant_id %in% unique(population_icds$participant_id)) %>%
  distinct() %>%
  group_by(participant_id) %>%
  summarise(N_diagnoses = n_distinct(ICD10_CLEAN1))


for (disease_index in c(11, 31, 39)) {

  disease_name <- map_disease_omim$NormalisedSpecificDisease[disease_index]
  disease_filename <- gsub(" ", "_", disease_name)
    
  selected_cases = population_participants %>%
    filter(normalised_specific_disease == disease_name) %>%
    filter(!duplicated(participant_id)) %>%
    mutate(Pheno = "Case")
  
  PRS = read.delim(paste0("phers_test/", disease_filename, "_phers.txt"))
  case_phers = selected_cases %>% 
    inner_join(PRS, by = "participant_id") %>%
    left_join(N_admissions) %>%
    left_join(N_diagnoses)
  
  print(disease_name)
  print(paste0("no. of cases: ", length(unique(selected_cases$participant_id))))
  print(paste0("N cases with phers > 0: ", length(unique(case_phers[case_phers$PRS > 0, ]$participant_id))))
  print(paste0("mean age: ", mean(case_phers[case_phers$PRS > 0, ]$age)))
  print(paste0("median record length: ", median(case_phers[case_phers$PRS > 0, ]$record_len)))
  print(paste0("IQR record length: ", IQR(case_phers[case_phers$PRS > 0, ]$record_len)))
  print(paste0("median no. of hes apc admissions: ", median(case_phers[case_phers$PRS > 0, ]$N_admissions)))
  print(paste0("IQR no. of hes apc admissions: ", IQR(case_phers[case_phers$PRS > 0, ]$N_admissions)))
  print(paste0("median no. of diagnosis codes in hes apc: ", median(case_phers[case_phers$PRS > 0, ]$N_diagnoses)))
  print(paste0("IQR no. of diagnosis codes in hes apc: ", IQR(case_phers[case_phers$PRS > 0, ]$N_diagnoses)))
  
  print(paste0("N cases with phers = 0: ", length(unique(case_phers[case_phers$PRS == 0, ]$participant_id))))
  print(paste0("mean age: ", mean(case_phers[case_phers$PRS == 0, ]$age)))
  print(paste0("median record length: ", median(case_phers[case_phers$PRS == 0, ]$record_len)))
  print(paste0("IQR record length: ", IQR(case_phers[case_phers$PRS == 0, ]$record_len)))
  print(paste0("median no. of hes apc admissions: ", median(case_phers[case_phers$PRS == 0, ]$N_admissions)))
  print(paste0("IQR no. of hes apc admissions: ", IQR(case_phers[case_phers$PRS == 0, ]$N_admissions)))
  print(paste0("median no. of diagnosis codes in hes apc: ", median(case_phers[case_phers$PRS == 0, ]$N_diagnoses)))
  print(paste0("IQR no. of diagnosis codes in hes apc: ", IQR(case_phers[case_phers$PRS == 0, ]$N_diagnoses)))

}
