### data

#### Mapping files - description of scripts and data sources used to create the final mapping files used in PheRS analysis.

###### Disease to OMIM

`~/scripts/create_omimdb.R`   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;       Create Disease to OMIM ID Map  
`omimdb.csv`                  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                  Disease-OMIM Mapping created and used - Normalised Specific Diseases | Phenotype Series | OMIM IDs   

###### OMIM to HPO

`~/scripts/prs_omim.py`      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Create OMIM to HPO Map, using:   
`phenotype_annotations.tab`  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    Mapping between OMIM ID and HPO terms, obtained from HPO website.  
`mapping_omim_to_hpo.csv`    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    HPO-OMIM Mapping used - cleaned version of above .tab file, containing only HPO to OMIM map  

###### HPO to Phecode

`hpo_phecode_mapping.csv`         &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    HPO to Phecode mapping, obtained from Bastarache 2018 paper Supplementary Table 12  

###### ICD10 to Phecode

`~/scripts/map_icdphec.py`        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    Create ICD10 to Phecode Map using:  
`ICD10list.csv`                   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of all ICD10 codes from UK Biobank, WHO dataset (encoding 19)  
`Phecode_map_v1_2_icd10_beta.csv` &nbsp; ICD10 to Phecode mapping file obtained from PheWAS catalog  
`mapping_icd10_to_phecode.csv`    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ICD10-Phecode Mapping used.  


