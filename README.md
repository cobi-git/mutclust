# COVID-19 Cytokine Analysis
## Overview
This repository contains the code and data used in our study on the cytokine profiles of COVID-19 patients. Our research focuses on identifying significant cytokines that can serve as biomarkers for disease severity and progression, with an emphasis on understanding the immune response during the deterioration and recovery phases of severe COVID-19 cases.

## Repository Contents
- data_preprocessing/: Scripts for preprocessing cytokine and clinical data.
- models/: Implementation of the Random Forest model used for analysis.
- analysis/: Scripts for calculating SHAP values and performing statistical tests.
- visualization/: Code for generating the figures and plots included in the manuscript.
- docs/: Documentation and supplementary materials related to the study.
## Data
The dataset used in this study includes cytokine data from plasma and laboratory data (e.g., neutrophil and lymphocyte counts) from blood samples, alongside clinical data from patients' Electronic Health Records (EHRs). The cytokine profiles were measured by the Korea National Institute of Health (KNIH) using the Luminex MAGPIX system with a customized panel, following a standardized protocol.

The dataset of the COVID-19 cohort used in this study is available in the Clinical and Omics Data Archive (CODA) database under the accession number CODA_D23017. Due to patient privacy concerns and data protection regulations, access to this data requires appropriate authorization.

## Code Availability
The code used for data preprocessing, Random Forest model implementation, SHAP value calculation, and visualization is available in this repository. To reproduce our results, please follow the instructions provided in the respective directories.

## Getting Started
To get started with the code in this repository, clone the repository to your local machine and install the required dependencies listed in requirements.txt.

