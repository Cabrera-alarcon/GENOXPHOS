# Landscape of respiratory complexes at an individual level.

The aim of this script is to assess genetic diversity within the Respiratory Complexes (RC) among humans, measured as the number of different RCs an individual is capable of assembling. To this end, we analyzed the genotypes of 2,504 individuals who participated in the third phase of the 1,000 Genomes Project (1). SNPs from the third phase were annotated with the VEP variant effect predictor (2) and nonsense variants affecting OxPhos genes were selected, considering onli aa changes that are going to be in the final structure (filterout those in the signal peptide). In addition, the script took into account stochiometric particularities of OxPhos complexes, as described in .......... 

## Installation

This script runs in R 4.1.2.

Install R with the following packages:

    vcfR
    ggplot2

## Run this script

This script takes two input files, also included in the folder:

1- A vcf.gz file including third phase missense variants from 1000 genome project, that affects to canonical transcripts.
2- A tsv file used to annotate variants, with gene symbol (intended to unify ATP-synthase nomenclature), startin aa included in the final structure according to Uniprot scifications (3).

## References

1. 1000 Genomes Project Consortium; Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, Abecasis GR. A global reference for human genetic variation. Nature. 2015 Oct 1;526(7571):68-74. doi: 10.1038/nature15393. PMID: 26432245; PMCID: PMC4750478.
2. McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Jun 6;17(1):122. doi: 10.1186/s13059-016-0974-4. PMID: 27268795; PMCID: PMC4893825.
3. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res. 2023 Jan 6;51(D1):D523-D531. doi: 10.1093/nar/gkac1052. PMID: 36408920; PMCID: PMC9825514.
