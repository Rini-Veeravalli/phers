# part II of variant analysis
# uses participant genetic information per extracted variant obtained in part I to perform genetic variant analysis
# e.g. association tests between PheRS and carrier status, case status and carrier status
# hardcoded for a specific disease, would need adjustment for batch processing of all diseases available 


#-------- Load dependencies

library(tidyverse)
library(readr)
library(dplyr)
library(stringr)
library(data.table)
install.packages("lifecycle")


#-------- Load data

samples <- read.csv("Brugada_syndrome_samples_ccratio_4.csv", row.names = 1)
sampleGT <- read.csv("samplesGT_Brugada.csv")
phers <- read.csv("Brugada_syndrome_phers.txt", header = TRUE, sep = '\t')

participants_genes <- read.csv("genome_file_paths_and_types_2020-05-01_00-19-10.csv")
participant_gene_id <- subset(participants_genes, select = c('Participant.Id', 'Platekey'))
participants <- fread("participant_2019-10-24_16-53-25.csv") %>%
    select(participant_id, year_of_birth, participant_phenotypic_sex, participant_type)
icds <- read.delim("hes_apc_clean_icds_v7.txt", na.strings = "NA", stringsAsFactors = F)
map_phe_icd10 <- read_delim("icd10_phe_mapping_cleaned.txt", col_types = "cc", delim = "\t") %>%
    distinct() %>%
    select("ICD10" = 1, PHECODE) %>%
    mutate(PhecodeCHR = as.character(as.numeric(PHECODE)))
population_icds <- icds %>%
    select(participant_id, ICD10_CLEAN1) %>%
    filter(participant_id %in% icds_admidate_pre2014$participant_id) %>%
    left_join(map_phe_icd10, by = c("ICD10_CLEAN1" = "ICD10"))
participants_age_sex <- population_icds %>%
    dplyr::select(participant_id) %>%
    distinct() %>%
    left_join(participants, by = "participant_id") %>%
    mutate(age = 2019 - year_of_birth,
           binaryGender = ifelse(participant_phenotypic_sex == "Male", 1, 0))


#-------- Find carriers and noncarriers

# find carriers with '0/1' in any columns (variants) in sampleGT 
carriers <- sampleGT %>%
    filter_all(any_vars(grepl('0/1', .)))

# find non-carriers who must not have '0/1' in all of the columns (variants)
noncarriers <- sampleGT %>%
    filter_all(all_vars(!grepl('0/1', .)))
    
# get participant genetic id of carriers and noncarriers
find_carriers <- substr(as.character(carriers[1]), 1, nchar(as.character(carriers[1])) - 4)
find_noncarriers <- substr(as.character(noncarriers[1]), 1, nchar(as.character(noncarriers[1])) - 4)

# match genetic ids to participant ids
carrier_participant_id <- unique(participant_gene_id[participant_gene_id$Platekey %in% unique(find_carriers), 'Participant.Id'])
noncarrier_participant_id <- unique(participant_gene_id[participant_gene_id$Platekey %in% unique(find_noncarriers), 'Participant.Id'])  


#-------- Association between case status and carrier status

#####---- Check case status of carriers and noncarriers

# check if carriers are cases or matched controls in samples
case_status_carrier <- samples[samples$participant_id %in% carrier_participant_id, ]

# check if noncarriers are cases or matched controls in samples
case_status_noncarrier <- samples[samples$participant_id %in% noncarrier_participant_id, ]

#####---- Check carrier status of cases and matched controls in samples

samples_carriers <- carrier_participant_id[carrier_participant_id %in% samples$participant_id]
carrier_status_samples <- samples[samples$participant_id %in% samples_carriers, ]
samples_noncarriers <- carrier_participant_id[noncarrier_participant_id %in% samples$participant_id]
noncarrier_status_samples <- samples[samples$participant_id %in% samples_noncarriers, ]

# add isCarrier binary column
carrier_status_samples$isCarrier = c('1')
noncarrier_status_samples$isCarrier = c('0')

# merge carrier_status of cases and controls in samples for logreg
carrier_status_samples <- rbind(carrier_status_samples, noncarrier_status_samples)
carrier_status_samples$isCase <- lapply(carrier_status_samples$isCase, as.numeric)
carrier_status_samples$isCarrier <- as.numeric(carrier_status_samples$isCarrier)   

# add age, gender to phers and carrier status dataset           
logregdata <- carrier_status_samples %>%
    inner_join(participants_age_sex, by = "participant_id") %>%
    dplyr::select(participant_id, age, binaryGender, participant_type, isCarrier, isCase) %>%
    distinct()

# Logistic Regression model 1
carrier_case_logreg <- glm(isCase ~ isCarrier, data = logregdata, family = binomial(link = "logit"))

# Logistic Regression model 2 with covariates
carrier_case_logreg_confounders <- glm(isCase ~ isCarrier + age + binaryGender, data = logregdata, family = binomial(link = "logit"))

# summary statistics, odds ratios and 95% CI
summary(carrier_case_logreg)
exp(cbind(OR = coef(carrier_case_logreg), confint(carrier_case_logreg)))

summary(carrier_case_logreg_confounders)
exp(cbind(OR = coef(carrier_case_logreg_confounders), confint(carrier_case_logreg_confounders)))


#-------- Association between PheRS and carrier status

carrier_phers <- phers
carrier_phers$isCarrier <- c('NA')
carrier_phers[carrier_phers$participant_id %in% carrier_participant_id, "isCarrier"] <- c('1')
carrier_phers[carrier_phers$participant_id %in% noncarrier_participant_id, "isCarrier"] <- c('0')
carrier_phers <- carrier_phers[(carrier_phers$isCarrier != "NA), ]
carrier_phers$isCarrier <- as.numeric(carrier_phers$isCarrier)

# add age, gender to phers and carrier status dataset           
logregdata <- carrier_phers %>%
    inner_join(participants_age_sex, by = "participant_id") %>%
    dplyr::select(participant_id, age, binaryGender, participant_type, PRS, isCarrier) %>%
    distinct()

# Logistic Regression model 1
carrier_prs_logreg <- glm(isCarrier ~ PRS, data = logregdata, family = binomial(link = "logit"))

# Logistic Regression model 2 with covariates
carrier_prs_logreg_confounders <- glm(isCarrier ~ PRS + age + binaryGender, data = logregdata, family = binomial(link = "logit"))

# summary statistics, odds ratios and 95% CI
summary(carrier_prs_logreg)
exp(cbind(OR = coef(carrier_prs_logreg), confint(carriers_prs_logreg)))

summary(carrier_prs_logreg_confounders)
exp(cbind(OR = coef(carrier_prs_logreg_confounders), confint(carrier_prs_logreg_confounders)))
