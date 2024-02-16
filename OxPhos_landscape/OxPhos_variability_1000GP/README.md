# Landscape of respiratory complexes at an individual level.

The aim of this script is to assess genetic diversity within the Respiratory Complexes (RC) among humans, measured as the number of different RCs an individual is capable of assembling. To this end, we analyzed the genotypes of 2,504 individuals who participated in the third phase of the 1,000 Genomes Project (1). SNPs from the third phase were annotated with the VEP variant effect predictor (2) and nonsense variants affecting OxPhos genes were selected, considering onli aa changes that are going to be in the final structure (filterout those in the signal peptide). In addition, the script took into account stochiometric particularities of OxPhos complexes. 

## Installation

This script runs in R 4.1.2.

Install R with the following packages:

    vcfR
    ggplot2

## Run this script

This script takes two input files, also included in the folder:
1- A vcf.gz file including third phase missense variants from 1000 genome project, that affects to canonical transcripts.
2- A tsv file used to annotate variants, with gene symbol (intended to unify ATP-synthase nomenclature), startin aa included in the final structure according to Uniprot scifications (3).

