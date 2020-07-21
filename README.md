# Phenotype Risk Scores

Phenotype Risk Scores (PheRS) in Genomics England (GeL). In order to evaluate whether the PheRS can detect undiagnosed Mendlian disorders, this project calculates the PheRS for rare diseases in GeL and analyses the results between cases and controls.  


## Structure

folder              | description
--------------------|---------------------------------------------------
scripts            | for analysis and creating datasets used
data               | data and mappings used and/or created
phers              | phers of all participants for diseases 
phers/description  | summary of codes mapped from each disease
samples            | case/controls identified for diseases
samples/summary    | summary of case/controls for diseases
diagrams           | boxplot and histograms to compare case/controls

## How to Use
Requires: GeL cohort datasets and ontology mappings in `/data`  
Input: adjudicated GeL rare diseases mapped to OMIM, specificed in `data/omimdb.csv`  
Run: `/scripts/phers.R` to calculate the PheRS for GeL rare disease cohort data in `/data`    
Output:  
-  the phenotype risk scores for each disease in `/phers`
-  descriptive information of the codes used for each disease in `/phers/description`
-  case and matched control samples for each disease in `/samples`
-  descriptive information about the case/control samples for each disease in `/samples/summary`
-  visualisations of case/control phers distributions for each disease in `/diagrams`


## Details

### scripts

`create_omimdb.R` 
generates adjudicated dataset of GeL rare diseases mapped to OMIM IDs, `/data/omimdb.csv`    


`phers.R`  
1. Load mappings and gel data  
1. Create functions to calculate PheRS:
    - create_filename
    - calculate_phers
    - create_case_control_sets
    - compare_case_controls
1. Run full analysis (functions above)


### data

`/data/omimdb.csv`     Normalised Specific Diseases | Phenotype Series | OMIM IDs  
