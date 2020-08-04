library(tidyverse)
library(icd)

icds = read.csv("/home/rveeravalli/re_gecip/health_records/PheRS/extract_ICD_codes/hes_apc_icds_v7_long_format.csv", stringsAsFactors = F, row.names = 1) 


tidy_icd = icds %>%
  select(-admidate) %>%
  # slice(1:100) %>% # for testing
  gather(DIAG, CODE, -participant_id) %>%
  select(-DIAG) %>% 
  distinct() %>% 
  filter(!is.na(CODE)) %>% 
  mutate(ICD10_CLEAN1 = str_extract(CODE, "[A-Z0-9]*"))

malformed = tidy_icd[!is_valid(as.icd10(tidy_icd$ICD10_CLEAN1),),] 
clean_icd = tidy_icd[is_valid(as.icd10(tidy_icd$ICD10_CLEAN1),),] 
# couldn't find function 'icd_is_valid.icd10' so used the above

# count individuals in clean_icd
clean_icd %>% 
  count(participant_id) %>% 
  nrow() #78957

#reformat ICD codes, include decimal divider
clean_icd = 
  clean_icd %>% 
  mutate(ICD10_decimal = short_to_decimal(clean_icd$ICD10_CLEAN1))
  # couldn't find function short_to_parts() so used the above
# check 'K20.X' 'R69.X6' 'R51.X' 'R31.X' 'N12.X' 'C61.X' 'O16.X' ?

### save clean icds
write.table(clean_icd, "hes_apc_clean_icds_v7.txt", row.names = F, col.names = T, sep="\t", quote = F)
