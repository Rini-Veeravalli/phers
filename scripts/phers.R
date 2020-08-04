# This script calculates the PheRS for rare disease participants in GeL
# and visualises the PheRS distribution of cases and controls for X normalised specific diseases.
# Input: disease list and associated OMIM IDs. 
# Output: summary of OMIM/HPO/Phecodes/ICD10 codes mapped to each disease, 
#         summary of case/control cohort data for each disease, 
#         visualisations of PheRS between case and controls for each disease. 
# Code developed by Stefanie Mueller and Rini Veeravalli.


# Load dependencies ---------------------------------

library(tidyverse)
library(data.table)
library(MatchIt)
library(pacman)
library(knitr)
library(tableone)
library(captioner)
library(ggpubr)
library(ggplot2)
# if(!require(devtools)) install.packages(devtools)
#  devtools::install_github("kassambara/ggpubr")
library(ggstatsplot)
library(optmatch)


# Load data and mappings ------------------------------

### TO DO change filepaths in final version

map_disease_omim = read.csv("omimdb.csv")

map_omim_hpo = read.delim("phenotype_annotation_new.tab", stringsAsFactors = F) %>%
  filter(DB == "OMIM") %>%
  select(DB_Object_ID, HPO_ID, Evidence_Code)

map_hpo_phe = read_csv("hpo_phecode_mapping.csv", col_types = "iccccccc")
# changing column names by replacing spaces with ..
names(map_hpo_phe)[2] = "HPO.term_id"
names(map_hpo_phe)[7] = "suppress."

map_phe_icd10 = read_delim("icd10_phe_mapping_cleaned.txt", col_types = "cc", delim = "\t") %>%
  distinct() %>%
  select("ICD10" = 1, PHECODE) %>%
  mutate(PhecodeCHR = as.character(as.numeric(PHECODE))) 
### PhecodeCHR preserves original phecode format i.e. removes .0


icds = read.delim("hes_apc_clean_icds_v7.txt", na.strings = "NA", stringsAsFactors = F)

disease_label = fread("rare_diseases_participant_dise_2019-09-26_16-44-06.csv") %>% 
  select(participant_id, starts_with("normalised")) 

participants = fread("participant_2019-10-24_16-53-25.csv") %>% 
  select(participant_id, programme, year_of_birth, participant_phenotypic_sex, participant_type)

icds_admidate_pre2014 = read.csv("participants_admidate_pre2014_updated.txt") %>% 
  select(participant_id) %>% 
  distinct()


# ICD10 data of Total participants with EHR data and admission date of > 5 years

population_icds = icds %>%
  select(participant_id, ICD10_CLEAN1) %>%
  filter(participant_id %in% icds_admidate_pre2014$participant_id) %>%
  left_join(map_phe_icd10, by = c("ICD10_CLEAN1" = "ICD10")) 

# Participant data of Total participants with EHR data and admission date of > 5 years and disease data

population_participants = population_icds %>%
  select(participant_id) %>%
  distinct() %>%
  left_join(participants, by = "participant_id") %>% 
  mutate(age = 2019 - year_of_birth,
         binaryGender = ifelse(participant_phenotypic_sex == "Male", 1, 0)) %>% 
  left_join(disease_label) %>%
  filter(!is.na(normalised_specific_disease))

# Number of Total participants with EHR data and admission date of > 5 years

ALLNUM = unique(population_icds$participant_id) %>%
  length()


# Create functions ------------------------------------

# create_filename()           # get disease name to use for naming files
# calculate_phers()           # calculate phers and create table of phers results, and create summary text, returns disease dataset
# create_case_control_sets()  # create case and control datasets
# compare_case_controls()     # create boxplot and histogram to compare cases and controls


create_filename <- function(disease_index) {
  disease_filename <- gsub(" ", "_", map_disease_omim$NormalisedSpecificDisease[disease_index])
  return(disease_filename)
}


calculate_phers <- function(disease_filename, disease_index) {
  disease = list()
  disease$OMIM = unique(unlist(strsplit(unlist(as.character(map_disease_omim$OMIM_ID[disease_index])), split = ",")))  
  disease$hpo = unique(map_omim_hpo[map_omim_hpo$DB_Object_ID %in% disease$OMIM, "HPO_ID"])
  disease$hpo_clean = str_extract(disease$hpo, "[0-9]*$") %>%
    as.integer()
  # exclude HPO terms: with NA phecode, not suppress 
  disease$phe = map_hpo_phe %>%
    filter(HPO.term_id %in% disease$hpo_clean) %>%
    filter(!is.na(phecode)) %>%
    filter(suppress. < 1) %>%
    select(2, 3, 4, 5) %>%
    distinct()
  disease$icds = map_phe_icd10 %>%
    filter(PhecodeCHR %in% disease$phe$phecode) %>%
    distinct()
  # phe and icds are not unique, instead tables with distinct rows
  
  # Calculate frequency of disease phecodes
  # Filter icds of population to disease icds and get associated disease phecodes

  disease$weight_phe = population_icds %>%
    filter(ICD10_CLEAN1 %in% disease$icds$ICD10) %>% 
    filter(PhecodeCHR %in% disease$icds$PhecodeCHR) %>%
    select(participant_id, PhecodeCHR) %>%
    distinct() %>%
    count(PhecodeCHR) %>%
    mutate(weight = log(ALLNUM / n))
  

  # Calculate PRS for all participants with at least one of ICDs of interest
  PRS_one = population_icds %>%
    inner_join(disease$weight_phe, by = "PhecodeCHR") %>%
    group_by(participant_id) %>%
    filter(!duplicated(PhecodeCHR)) %>%
    summarise(PRS = sum(weight))
  
  # Fill up missing PRS, owning to no ICD hit, with zeros
  PRS_two = population_icds %>%
    select(participant_id) %>%
    filter(!participant_id %in% PRS_one$participant_id) %>%
    distinct() %>%
    mutate(PRS = 0)
  
  # Bind both databases and write out
  PRS = rbind(PRS_one, PRS_two) %>%
    arrange(participant_id) %>%
    mutate(OMIM = paste0(disease$OMIM, collapse = ","))
  
  # Save PheRS results
  write.table(PRS, paste0("phers/", disease_filename, "_phers.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  
  # Create and save summary of uniqe disease-related codes
  outSummary <- setNames(data.frame(matrix(ncol = 13, nrow = 0)), 
                         c("Disease", "number_OMIM", "number_HPO_terms", "number_Phecodes",
                           "number_observed_Phecodes", "number_ICD", "max_PRS_theo", "max_PRS_obs", "percentage_zero"))
  outSummary = data.frame(Disease = map_disease_omim$NormalisedSpecificDisease[disease_index],
                          number_OMIM = length(unique(disease$OMIM)),
                          number_HPO_terms = length(unique(disease$hpo_clean)),
                          number_Phecodes = length(unique(disease$phe$phecode)),
                          number_observed_Phecodes = length(unique(disease$phe$phecode) %in% icds$PHECODE),
                          number_ICD = length(unique(disease$icds$ICD10)),
                          max_PRS_theo = sum(disease$weight_phe$weight),
                          max_PRS_obs = max(PRS$PRS),
                          percentage_zero = (nrow(PRS_two) / ALLNUM)
  )
  write.csv(outSummary, paste0("summary/", disease_filename, "_summary.txt"), row.names = FALSE, append = TRUE)

  # save disease list as R object for future evaluation of which codes mapped to disease
  save(disease, file = paste0("summary/", disease_filename, "_codes.RData"))
  
  return(disease)
} 


create_case_control_sets <- function(disease_filename, disease_index, ccratio) {
  disease_interest = disease_label %>%
    filter(normalised_specific_disease == map_disease_omim$NormalisedSpecificDisease[disease_index]) %>% 
    select(participant_id, normalised_specific_disease)
  
  selected_cases = disease_interest %>% 
    inner_join(population_participants) %>% 
    filter(!duplicated(participant_id)) %>%
    mutate(Pheno = "Case")

  any(duplicated(selected_cases$participant_id))
  
  demographics_cases = selected_cases %>% 
    summarise(Avg_age = mean(age),
              PercFem = mean(binaryGender),
              Count = n())

  
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
  
  together = rbind(selected_cases[, c("participant_id", "age", "binaryGender", "Pheno")],
                   possible_controls[, c("participant_id", "age", "binaryGender", "Pheno")]) %>% 
    mutate(isCase = ifelse(Pheno == "Case", T, F)) %>% 
    mutate(factorGender = ifelse(binaryGender == 1, "male", "female") %>%
             as.factor)
  
  
  # --- match controls
  set.seed(1)
  options("optmatch_max_problem_size" = Inf)
  match <- matchit(isCase ~ age + factorGender, data = together, method = "optimal", ratio = ccratio)
  
  # Save descriptive summary of case/control samples
  capture.output(summary(match), file = paste0("samples/", disease_filename,"_samples_summary_ccratio_", ccratio,".txt"))
  
  df.match <- match.data(match)[1:ncol(together)]
  
  # --- write out case/control samples
  df.match %>% 
    write.csv(paste0("samples/", disease_filename,"_samples_ccratio_", ccratio,".csv"))
  
}


compare_case_controls <- function(disease_filename, disease, disease_index, ccratio) {
  
  # --- combine case/control samples with their phers
  sam = read.csv(paste0("samples/", disease_filename,"_samples_ccratio_", ccratio,".csv"))
  PRS = read.delim(paste0("phers/", disease_filename, "_phers.txt"))
  sampPheRS = sam %>% 
    inner_join(PRS, by = "participant_id") 
  
 
  # --- Non-parametric comparison between phers and case status
  wilcox.test(PRS ~ Pheno, data = sampPheRS)

 
  
  # --- Plot
  # set color levels
  col <- c("grey30", "#045a8d")
  
  
  sampPheRS = sampPheRS %>% 
    mutate(PhenotypeStatus = plyr::revalue(Pheno,
                                        c("Case" = paste0('Case \n','n=', nrow(sampPheRS[sampPheRS$Pheno == "Case", ])),
                                          "Control" = paste0('Control', '\n', 'n=', nrow(sampPheRS[sampPheRS$Pheno == "Control", ])))))
  
  mean.data = sampPheRS %>% 
    group_by(PhenotypeStatus) %>% 
    summarise(mean = mean(PRS, na.rm = T),
              label = paste("list(~italic(mu) ==", round(mean, digits = 3), ")"))

  boxplot = sampPheRS %>% 
    ggplot(aes(y = PRS, x = PhenotypeStatus, col = PhenotypeStatus)) + 
    theme_bw() +
    geom_point(position = position_jitter(width = 0.3), alpha = 0.3, size = 3) +
    geom_boxplot(width = 0.8, alpha = 0.2, fill = "white",
                 outlier.shape = NA, color = "grey20", size = 1.5,
                 position = ggplot2::position_dodge(width = NULL)) +
    stat_summary(fun.y = mean, geom = "point", size = 4, color = "red") +
    ggrepel::geom_label_repel(data = mean.data, aes(x = PhenotypeStatus, y = mean, label = label),
                              color = "black", parse = TRUE, nudge_x = -0.3) + 
    scale_colour_manual(values = col, guide = F) +
    labs(y = "Phenotype Risk Score", x = "Phenotype Status",
         title = paste0("Phenotype / Disease Label: ", map_disease_omim$NormalisedSpecificDisease[disease_index]),
         subtitle = paste0("OMIM ID: ", paste(disease$OMIM, collapse = ","))) +
    theme(title        = element_text(size = rel(0.5), face = "bold"),
          axis.title.y = element_text(vjust = 1.5),
          axis.title.x = element_text(vjust = -1.5),
          axis.text    = element_text(size = rel(1.25)),
          legend.text  = element_text(size = rel(1.3)),
          strip.text   = element_text(size = rel(1.3)),
          plot.margin  = unit(c(1, 1, 1, 2), "cm"),
          panel.grid.major = element_line(colour = "grey80")) +
    stat_compare_means(comparisons = list(c(1, 2)), method = "wilcox.test",
                       label = "p.signif", bracket.size = 0.5, text.size = 2) + 
    stat_compare_means(label.y = max(sampPheRS$PRS) * 1.2) 

  
  boxplot %>% 
    ggsave(filename = paste0("boxplot_test/", disease_filename, "_ccratio_", ccratio, ".pdf"),
           width = 9, height = 6)
  
  # --- Histogram
  ggsave(
    (sampPheRS %>%
       ggplot(aes(x = PRS, fill = PhenotypeStatus, col = PhenotypeStatus)) + 
       theme_bw() +
       geom_histogram() +
       scale_colour_manual(values = col) +
       scale_fill_manual(values = col) +
       facet_grid(PhenotypeStatus~., scales = "free_y") +
       labs(x = paste0("PheRS of ", map_disease_omim$NormalisedSpecificDisease[disease_index])) +   
       theme(title        = element_text(size = rel(1.3), face = "bold"),
             axis.title.y = element_text(vjust = 1.5),
             axis.title.x = element_text(vjust = -1.5),
             axis.text    = element_text(size = rel(1.25)),
             legend.text  = element_text(size = rel(1.3)),
             strip.text   = element_text(size = rel(1.3)),
             plot.margin  = unit(c(1, 1, 1, 2), "cm"),
             panel.grid.major = element_line(colour="grey80"))),
    
    filename = paste0("hist_test/", disease_filename, "_ccratio_", ccratio, "_hist.pdf"),
    width = 9, height = 6
  )  
}



# Run main script -------------------------------------

# ccratio specifies ratio of controls to cases
ccratio = 4
for (row in 1:nrow(map_disease_omim)) {
  tryCatch({
    disease_filename = create_filename(map_disease_omim$NormalisedSpecificDisease[row])
    disease = calculate_phers(disease_filename, row)
    create_case_control_sets(disease_filename, row, ccratio)
    compare_case_controls(disease_filename, disease, row, ccratio)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}


