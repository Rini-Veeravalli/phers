# Script to carry out association tests between PheRS and Case Status


#-------- Load libraries

library(data.table)
install.packages("devtools")
devtools::install_github("tidyverse/dplyr")
library(dplyr)
install.packages("caret")
library(caret)


#-------- Load data

icds = read.delim("~/re_gecip/health_records/PheRS/extract_ICD_codes/hes_apc_clean_icds_v7.txt", na.strings = "NA", stringsAsFactors = F)

disease_label = fread("~/Documents/rare_diseases_participant_dise_2019-09-26_16-44-06.csv") %>% 
  select(participant_id, starts_with("normalised")) 

participants = fread("~/Documents/participant_2019-10-24_16-53-25.csv") %>% 
  select(participant_id, programme, year_of_birth, participant_phenotypic_sex, participant_type)

icds_admidate_pre2014 = read.csv("participants_admidate_pre2014_updated.txt") %>% 
  select(participant_id) %>% 
  distinct()

map_phe_icd10 = read_delim("~/re_gecip/health_records/PheRS/ref_data/icd10_phe_mapping_cleaned.txt", col_types = "cc", delim = "\t") %>%
  distinct() %>%
  select("ICD10" = 1, PHECODE) %>%
  mutate(PhecodeCHR = as.character(as.numeric(PHECODE))) 

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

# Number of total participants with EHR data and admission date of > 5 years

ALLNUM = unique(population_icds$participant_id) %>%
  length()


#-------- Create functions

create_filename <- function(disease_index) {
  disease_filename <- gsub(" ", "_", map_disease_omim$NormalisedSpecificDisease[disease_index])
  return(disease_filename)
}

test_association <- function(disease_filename, disease_index, ccratio) {
  print(disease_filename)
  disease_name = map_disease_omim$NormalisedSpecificDisease[disease_index]
  
  # --- combine case/control samples with their phers
  sam = read.csv(paste0("samples_test/", disease_filename,"_samples_ccratio_", ccratio,".csv"))
  PRS = read.delim(paste0("phers_test/", disease_filename, "_phers.txt"))
  sampPheRS = sam %>% 
    inner_join(PRS, by = "participant_id") 

  # --- Non-parametric comparison between case and control phers
  compare_means <- wilcox.test(PRS ~ Pheno, data = sampPheRS)
  zscore <- qnorm(compare_means$p.value/2)

  # Create and save summary of test results
  summary_wilcox <- data.frame(Disease = disease_name, WilcoxRS_Test_p_value = compare_means$p.value, WilcoxRS_Test_z_score = zscore)
  write.csv(summary_wilcox, paste0("analysis_test/", disease_filename, "_wilcox.txt"), row.names = FALSE)
  
  
  # --- Investigate relationship

  # Logistic reg - does PheRS affect case status?
  
  selected_cases = population_participants %>%
    filter(normalised_specific_disease == map_disease_omim$NormalisedSpecificDisease[disease_index]) %>%
    filter(!duplicated(participant_id)) %>%
    mutate(Pheno = "Case")
    
  # --- remove participants with diagnosis from possible controls before matching
  possible_controls = population_participants %>% 
    # remove all participants_ids of participants in selected_cases
    filter(!participant_id %in% selected_cases$participant_id) %>%
    # remove all partcipant_ids with diseases in selected_cases
    filter(!normalised_disease_group %in% selected_cases$normalised_disease_group) %>% 
    filter(!normalised_specific_disease %in% selected_cases$normalised_specific_disease) %>% 
    distinct() %>%
    filter(!duplicated(participant_id)) %>% 
    mutate(Pheno = "Control")
    
  popPheRS = rbind(selected_cases[, c("participant_id", "age", "binaryGender", "record_len", "Pheno")],
                   possible_controls[, c("participant_id", "age", "binaryGender", "record_len", "Pheno")]) %>% 
    mutate(Pheno = ifelse(Pheno == "Case", 1, 0) %>%
             as.factor) %>% 
    inner_join(PRS, by = "participant_id")
 
  require(MASS)
  popPheRS.fit <- glm(Pheno ~ PRS, data = popPheRS, family = binomial(link = "logit"))
  popPheRS_confounding.fit <- glm(Pheno ~ PRS + age + binaryGender + record_len, data = popPheRS, family = binomial(link = "logit"))
  
  
  summary_logreg_prs <- data.frame(Disease    = disease_name,
                                   p_value    = coef(summary(popPheRS.fit))["PRS", 4],
                                   OR         = exp(coef(popPheRS.fit)["PRS"]),
                                   OR_CI_2_5  = exp(confint(popPheRS.fit)[2, 1]),
                                   OR_CI_97_5 = exp(confint(popPheRS.fit)[2, 2]))

  summary_logreg_confounders <- data.frame(Disease                 = disease_name,
                                           PRS_p_value             = coef(summary(popPheRS_confounding.fit))["PRS", 4],
                                           age_p_value             = coef(summary(popPheRS_confounding.fit))["age", 4],
                                           binaryGender_p_value    = coef(summary(popPheRS_confounding.fit))["binaryGender", 4],
                                           record_len_p_value      = coef(summary(popPheRS_confounding.fit))["record_len", 4],
                                           PRS_OR                  = exp(coef(popPheRS_confounding.fit)["PRS"]),
                                           age_OR                  = exp(coef(popPheRS_confounding.fit)["age"]),
                                           binaryGender_OR         = exp(coef(popPheRS_confounding.fit)["binaryGender"]),
                                           record_len_OR           = exp(coef(popPheRS_confounding.fit)["record_len"]),
                                           PRS_OR_CI_2_5           = exp(confint(popPheRS_confounding.fit)["PRS", 1]),
                                           age_OR_CI_2_5           = exp(confint(popPheRS_confounding.fit)["age", 1]),
                                           binaryGender_OR_CI_2_5  = exp(confint(popPheRS_confounding.fit)["binaryGender", 1]),
                                           record_len_OR_CI_2_5    = exp(confint(popPheRS_confounding.fit)["record_len", 1]),
                                           PRS_OR_CI_97_5          = exp(confint(popPheRS_confounding.fit)["PRS", 2]),
                                           age_OR_CI_97_5          = exp(confint(popPheRS_confounding.fit)["age", 2]),
                                           binaryGender_OR_CI_97_5 = exp(confint(popPheRS_confounding.fit)["binaryGender", 2]),
                                           record_len_OR_CI_97_5   = exp(confint(popPheRS_confounding.fit)["record_len", 2]))

  
  write.csv(summary_logreg_prs, paste0("analysis_test/", disease_filename, "_association_phers.txt"), row.names = FALSE)
  write.csv(summary_logreg_confounders, paste0("analysis_test/", disease_filename, "_association_confounders.txt"), row.names = FALSE)
 
  return(popPheRS)
}

test_prediction <- function(disease_filename, disease_index, threshold) {
  
  # --- Investigate the predictive ability of PheR - does PheRS predict case status?
  
  # paper: high-scoring controls = PheRS is > 3rd quartile PheRS for cases
  # here: high-scorers = PheRS is in the 99th quantile
  
  disease_filename = create_filename(map_disease_omim$NormalisedSpecificDisease[disease_index])
  PRS = read.delim(paste0("phers_test/", disease_filename, "_phers.txt"))
  disease_name = map_disease_omim$NormalisedSpecificDisease[disease_index]
  
  recruited_cases = population_icds %>%
    inner_join(disease_label) %>%
    filter(normalised_specific_disease == disease_name) %>%
    filter(!duplicated(participant_id)) %>%
    mutate(Pheno = "Case")
  
  possible_controls = population_icds %>% 
    inner_join(disease_label) %>%
    # remove all participants_ids of participants in recruited_cases
    filter(!participant_id %in% recruited_cases$participant_id) %>%
    # remove all partcipant_ids with diseases in recruited_cases
    # filter(!normalised_disease_group %in% recruited_cases$normalised_disease_group) %>% 
    # filter(!normalised_specific_disease %in% recruited_cases$normalised_specific_disease) %>% 
    distinct() %>%
    filter(!duplicated(participant_id)) %>% 
    mutate(Pheno = "Control")
  
  allPheRS <- rbind(recruited_cases[, c("participant_id", "normalised_specific_disease", "Pheno")],
                   possible_controls[, c("participant_id", "normalised_specific_disease", "Pheno")]) %>% 
    mutate(Pheno = ifelse(Pheno == "Case", 1, 0) %>%
             as.factor) %>% 
    inner_join(PRS, by = "participant_id") 
  
 
  case_threshold <- quantile(PRS$PRS, threshold)

  # get info of participants in the 99th quantile
  # check if they are cases or controls
  summary(allPheRS[allPheRS$PRS >= case_threshold, ]$Pheno)
  
  allPheRS <- allPheRS %>%
    mutate(predicted = ifelse(PRS >= case_threshold, 1, 0) %>%
             as.factor)


  confmat <- confusionMatrix(data = allPheRS$predicted, reference = allPheRS$Pheno, positive = "1")

  summary_pred_stats <- data.frame(disease                 = disease_name,
                                   prevalence              = confmat$byClass[['Prevalence']],
                                   quantile_threshold      = quantile_threshold,
                                   case_threshold          = case_threshold,
                                   sensitivity             = confmat$byClass[['Sensitivity']],
                                   specificity             = confmat$byClass[['Specificity']],
                                   posPredValue            = confmat$byClass[['Pos Pred Value']],
                                   negPredValue            = confmat$byClass[['Neg Pred Value']],
                                   accuracy                = confmat$overall[['Accuracy']],
                                   precision               = confmat$byClass[['Precision']],
                                   recall                  = confmat$byClass[['Recall']],
                                   f1                      = confmat$byClass[['F1']]
                                   )

  write.csv(summary_pred_stats, paste0("analysis_test/", disease_filename, "_threshold_", quantile_threshold, "_predictive_summary.txt"), row.names = FALSE)
  # hist(allPheRS[allPheRS$PRS >= case_threshold, ]$PRS)
}

create_all_disease_summaries <- function() {
  list_of_files <- list.files(path="analysis_test", pattern = "wilcox*.txt", full.names = TRUE) 
  summary_wilcox <- rbindlist(sapply(list_of_files, fread, simplify = F),
                              use.names = T, idcol = "FileName", fill = T)
  write.csv(summary_wilcox, paste0("analysis_test/all_disease_wilcox_summary.csv"), row.names = FALSE)
  
  list_of_files <- list.files(path="analysis_test", pattern = "association_phers*.txt", full.names = TRUE) 
  summary_association_phers <- rbindlist(sapply(list_of_files, fread, simplify = F),
                                         use.names = T, idcol = "FileName", fill = T)
  write.csv(summary_association_phers, paste0("analysis_test/all_disease_association_phers_summary.csv"), row.names = FALSE)
  
  list_of_files <- list.files(path="analysis_test", pattern = "association_confounders*.txt", full.names = TRUE) 
  summary_association_confounders <- rbindlist(sapply(list_of_files, fread, simplify = F),
                                               use.names = T, idcol = "FileName", fill = T)
  write.csv(summary_association_confounders, paste0("analysis_test/all_disease_association_confounders_summary.csv"), row.names = FALSE)
  
  list_of_files <- list.files(path="analysis_test", pattern = "predictive_summary*.txt", full.names = TRUE) 
  summary_association_confounders <- rbindlist(sapply(list_of_files, fread, simplify = F),
                                               use.names = T, idcol = "FileName", fill = T)
  write.csv(summary_association_confounders, paste0("analysis_test/all_disease_predictive_stats_summary.csv"), row.names = FALSE)
}




# main script

for (row in 1:nrow(map_disease_omim)) { #c(11, 31, 39)){ # 3 test diseases BS, DCM, NF1
  tryCatch({ 
    ccratio = 4
    case_threshold = 0.99
    disease_filename = create_filename(map_disease_omim$NormalisedSpecificDisease[row])
    sampPheRS = test_association(disease_filename, row, ccratio)
    test_prediction(disease_filename, row, case_threshold)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
create_all_disease_summaries()
