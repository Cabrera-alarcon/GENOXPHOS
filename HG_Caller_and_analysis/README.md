# Mitochondrial Haplogroup Calleing
## HG calling using machine Random Forest Model.

We developed a custom method for mitochondrial haplogroup (HG) calling from GWAS array data, training a Random Forest classifier. Approximately 189 of the 231 positions were associated with haplogroup (HG) markers, and the probes defining the genotype in our array met the necessary quality standards. Our objective was to train a random forest classifier using these 189 positions (which serve as the variables for the random forest) derived from 61,134 HG-labeled sequences obtained from the Mitomap database. The perfomance of our model was measured as the mean of Cohen’s kappa calculated for each cross-validation cycle (3-fold CV). Aditionally, external validation of this model was carried out using data from the third phase of 1,000 genome project.<br><br>
<p align="center">
  <img src="https://github.com/Cabrera-alarcon/GENOXPHOS/blob/master/HG_Caller_and_analysis/Train_validation_workflow.png" width="1000" title="hover text">
</p>

<br><br>
An app powered by shiny (https://shiny.posit.co/) to do mitochondrial haplogroup calling from GWAS array data will be available here soon.

## Analysis of HGs as independent risk factors

We evaluate the potential of haplogroups (HGs) as independent risk factors for the severity of SARS-CoV-2. Severe disease was defined as cases with a fatal outcome, admission to the ICU, or the need for mechanical ventilation (invasive or non-invasive). HGs assessed were those that were represented in at least 1% o the patients used in the study. In addition HGs H, V, and HV were grouped under the HV branch, given that they share same residue changes determined by top-level HG-markers (present in ≥80% of HGs).
We then assessed the explanatory significance of HGs for SARS-CoV-2 severity in two groups:<br><br>
1. The script Analysis_1_scourge.R study the role of mitochondrial HGs as risk facro of severity in the SCOURGE cohort (group with comorbidity data available). Here we analyze mitochondrial HGs as independent risk factors of SARS-CoV-2 severity from comorbidities, sex, age and genetic background.<br><br>
2. The Analysis_2_all_patients.R script explores HG as an independent risk factor in all patients, cosidering sex, age and genetic background (here comorbidity data were not available for all individuals).<br><br>
## Contact: 
jaenriquez@cnic.es or jlcabreraa@cnic.es 

