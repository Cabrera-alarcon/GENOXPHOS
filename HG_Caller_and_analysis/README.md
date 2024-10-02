# Mitochondrial Haplogroup Calleing
## HG calling using machine Random Forest Model.

We developed a custom method for mitochondrial haplogroup (HG) calling from GWAS array data, training a Random Forest classifier. Approximately 189 of the 231 positions were associated with haplogroup (HG) markers, and the probes defining the genotype in our array met the necessary quality standards. Our objective was to train a random forest classifier using these 189 positions (which serve as the variables for the random forest) derived from 61,134 HG-labeled sequences obtained from the Mitomap database. The perfomance of our model was measured as the mean of Cohen’s kappa calculated for each cross-validation cycle (3-fold CV). Aditionally, external validation of this model was carried out using data from the third phase of 1,000 genome project.

<p align="center">
  <img src="https://github.com/Cabrera-alarcon/VSRFM/blob/VSRFM/ROC_curves.png" width="1000" title="hover text">
</p>


Variant pathogenic prediction models VSRFM and VSRFM-s, the importance of splicing and allele frequency. Cabrera-Alarcon JL, Garcia-Martinez J. Available in bioRxiv: https://www.biorxiv.org/content/early/2018/10/03/430975?rss=1

## Contact: 
joseluiscabreraalarcon@gmail.com or xurde.garcia.martinez@gmail.com 

