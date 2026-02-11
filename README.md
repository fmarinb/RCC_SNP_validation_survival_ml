#  Validation of ITPR2, DPF3, EPAS1 and PVT1 -associated SNPs as biomarkers for RCC in an independent case-control cohort.

This repository contains:

1. The R scripts used for the analysis of overall survival (OS) and 5-year survival from the date of diagnosis, consisting of:
1.1 The generation of multiple Cox proportional hazard models, including the SNPs studied quantified by qPCR and clinical and sociodemographic covariates (age group and sex of the patient). Various models were adjusted according to the relationship between predictors (additive model, interaction model, model with stratification of metastasis status at diagnosis).
1.2 Generation of Kaplan-Meier curves as a univariate strategy for evaluating the prognostic value of the SNPs of interest, following an interaction of covariates that best fit once the Cox models generated in 1.1 had been evaluated.

2. The corresponding R scripts with the machine learning models used for the selection of the genes studied according to their expression levels quantified by RT-PCR on paired samples of FFPE tumour and adjacent healthy tissue. The corresponding training, testing and evaluation functions were generated using metrics such as Accuracy, Sensitivity, Specificity, Precision, Recall and F1-score.

## Abstract
Introduction: Renal cell carcinoma (RCC) is a heterogeneous malignancy influenced by genetic and environmental factors. Previous genome-wide association studies (GWAS) have identified risk single nucleotide polymorphisms (SNPs) associated with RCC susceptibility, particularly within genes such as ITPR2, DPF3, EPAS1, PVT1, and MYC. These SNPs are in regions implicated in key cellular processes like calcium signaling, chromatin remodeling, hypoxia response and oncogenesis. These pathways are highly relevant to RCC pathogenesis, although the functional significance of these genetic variations in sporadic RCC remains insufficiently characterized. 

Methods: This study analyzed five GWAS-identified SNPs—rs1049380 and rs10771279 (ITPR2), rs4903064 (DPF3), rs7579899 (EPAS1), and rs35252396 (PVT1/MYC)—in a Spanish case-control cohort comprising 168 RCC patients and 259 healthy controls. Genotyping was performed from buccal swabs, and gene expression levels were assessed in 33 paired formalin-fixed paraffin-embedded (FFPE) tumor and adjacent normal kidney tissue samples. Associations between SNPs, overall survival, and expression of quantitative trait loci (eQTLs) were evaluated in relation to RCC risk and RCC progression in the case of survival curves.

Results: The C/C genotype of ITPR2 rs10771279 was nominally associated with a protective effect (OR:0.41), with higher ITPR2 expression observed in healthy tissues than in RCC. The C/C genotype of DPF3 rs4903064 was nominally correlated with increased RCC risk (OR: 2.21) and higher DPF3 expression, potentially linked to hypoxia-inducible pathways. Similarly, EPAS1 rs7579899 A/A genotype was nominally associated with RCC risk (OR:1.78) While PVT1/MYC rs35252396 did not show susceptibility relation, both genes showed upregulated expression in RCC tissue.
In survival analyses, the G allele of rs1049380 (ITPR2) was significantly associated with reduced  5-year survival in metastic and non-metastatic patients. Additionally, the AC genotype of rs35252396 showed nominal associations with  highest risk in 5-year survival models.

Conclusion:  This study provides independent evidence supporting the biological relevance of GWAS-identified loci in RCC. While several variants showed nominal associations with disease risk, ITPR2 rs1049380 emerged as a variant of potential prognostic relevance for five-year overall survival. Overall, these findings highlight the differential contribution of genetic variants to RCC susceptibility and progression and should be considered hypothesis-generating, warranting validation in larger, independent cohorts.

## Data avaliability

The datasets supporting this study (RCC_genotype_clin_survival.xlsx, RCC_genotype_expression_HvT.xlsx, and RCC_genotype_clin_PT_HC.xlsx), containing genotyping, clinical, sociodemographic, survival, and RT-qPCR gene expression data (RQ, 2⁻ΔΔCt), have been deposited in Zenodo (DOI: 10.5281/zenodo.18600366) and are available under restricted access upon reasonable request.

