# Mitochondrial Haplogroup Calleing
## HG calling using machine Random Forest Model.

We developed a custom method for mitochondrial haplogroup (HG) calling from GWAS array data, training a Random Forest classifier. Approximately 189 of the 231 positions were associated with haplogroup (HG) markers, and the probes defining the genotype in our array met the necessary quality standards. Our objective was to train a random forest classifier using these 189 positions (which serve as the variables for the random forest) derived from 61,134 HG-labeled sequences obtained from the Mitomap database (https://www.mitomap.org/foswiki/bin/view/Main/WebHome). The perfomance of our model was measured as the mean of Cohen’s kappa calculated for each cross-validation cycle (3-fold CV). Aditionally, external validation of this model was carried out using data from the third phase of 1,000 genome project (https://www.internationalgenome.org).<br><br>
<p align="center">
  <img src="https://github.com/Cabrera-alarcon/GENOXPHOS/blob/master/HG_Caller_and_analysis/Train_validation_workflow.png" width="1000" title="hover text">
</p>

<br><br>
Trained model is available as HG_caller.rds.<br><br>
An app powered by shiny (https://shiny.posit.co/) to do mitochondrial haplogroup calling from GWAS array data is available here as HG_caller_v1.0.tar.gz. This app runs in R version 4.4.1 and the following R-packages:
* caret
* dplyr
* irr
* ggplot2
* shiny
* shinyFiles
* shinythemes

## Analysis of HGs as independent risk factors

We evaluate the potential of haplogroups (HGs) as independent risk factors for the severity of SARS-CoV-2. Severe disease was defined as cases with a fatal outcome, admission to the ICU, or the need for mechanical ventilation (invasive or non-invasive). HGs assessed were those that were represented in at least 1% o the patients used in the study. In addition HGs H, V, and HV were grouped under the HV branch, given that they share same residue changes determined by top-level HG-markers (present in ≥80% of HGs).
Then, we assessed the explanatory significance of HGs for SARS-CoV-2 severity in two groups:<br><br>
1. The script Analysis_1_scourge.R study the role of mitochondrial HGs as risk facro of severity in the SCOURGE cohort (group with comorbidity data available). Here we analyze mitochondrial HGs as independent risk factors of SARS-CoV-2 severity from comorbidities, sex, age and genetic background.<br><br>
2. The Analysis_2_all_patients.R script explores HG as an independent risk factor in all patients, cosidering sex, age and genetic background (here comorbidity data were not available for all individuals).<br><br>
## Citation

If you use this program in your research, please cite the following publication:

Cabrera-Alarcon JL, Cruz R, Rosa-Moreno M, Latorre-Pellicer A, de Almeida SD; Scourge Cohort Group; Riancho JA, Rojas-Martinez A, Flores C, Lapunzina P, Sánchez-Cabo F, Carracedo Á, Enriquez JA.  
**Shaping current European mitochondrial haplogroup frequency in response to infection: the case of SARS-CoV-2 severity.**  
*Commun Biol.* 2025 Jan 9;8(1):33. doi: [10.1038/s42003-024-07314-y](https://doi.org/10.1038/s42003-024-07314-y). PMID: 39789223.
individuals).<br><br>
## Contact: 
jaenriquez@cnic.es or jlcabreraa@cnic.es 

