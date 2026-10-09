# IntOGen RELEASE 31/05/2023
# Headers of the cohorts files
# For further information about methods please visit the documentation website: https://intogen.readthedocs.io/en/latest/index.html
# To access the source code of intogen please visit: https://bitbucket.org/intogen/intogen-plus/src/master/


Cohorts.tsv


1   COHORT: Name of the cohort. 
2   CANCER: Acronym of the cancer type associated with the cohort. 
3   CANCER_NAME: Long name of the cancer type associated with the cohort. 
3   SOURCE: Source of the data (TCGA, PCAWG, HARTWIG, etc.)
4   PLATFORM: Whole-exome sequencing (WXS) or Whole-genome sequencing (WGS). 
6   REFERENCE: Pubmed ID of the publication. 
7   TYPE: Type of cohort {“Primary”, “Metastatic”, “Relapse”}
8   TREATED: Treatment status {“Treated”, “Untreated”}
9   AGE: Age status {“Adult”, “Pediatric”}
10  SAMPLES: Number of samples (before filtering by intOGen). 
11  MUTATIONS: Number of total of mutations in the cohort (before filtering by intOGen). 
12  COHORT_NICK:  Short name in the intOGen website. 
13  COHORT_NAME: Long name in the intOGen website.
14  DRIVERS: Number of drivers gene detected per cohort