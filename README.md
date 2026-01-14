# A predicted cancer dependency map for paralog pairs #

This repository contains the scripts and Jupyter notebooks for the research paper titled:
**"A predicted cancer dependency map for paralog pairs"**, currently available in XXX

The repository includes:
- `notebooks/01_preprocessing/` - Preprocess the input files from CRISPR screens and DepMap portal
- `notebooks/02_feature_calculation/` - Compute network based features for paralog pairs
- `notebooks/03_feature_annotation/`- Annotate and integrate the calculated features into training and test datasets
- `notebooks/04_model_evaluation/` - Train and evaluate the context-specific random forest classifier
- `notebooks/05_visualization/` - Generate all main and supplementary figures reported in the manuscript
- `scripts/` - Python scripts to implement the prediction pipeline for annotating paralog pairs and generating prediction scores

## Overview of the Notebooks

1. Data Preprocessing

All the input data were processed in these notebooks

| Notebook                       | Description                                                                                  |
|:-------------------------------|:---------------------------------------------------------------------------------------------|
| preprocess_ito.ipynb           | Processing of LFC and FDR datasets from the CRISPR screen by Ito et al.|
| preprocess_klingbeil.ipynb     | Processing of LFC and FDR datasets from the CRISPR screen by Klingbeil et al. |
| preprocess_parrish.ipynb       | Processing of LFC dataset from the CRISPR screen by Parrish et al. |
| preprocess_DepMap_22Q4.ipynb   | Processing of the transcriptomics, gene effect scores, copy number and somatic mutation profiles from the DepMap22Q4 release|   

2. Feature Calculation and Annotation

All the features calculated and mapped using the notebooks below

| Notebook                                                            | Description                                                                             |
|:--------------------------------------------------------------------|:----------------------------------------------------------------------------------------|
| 02_feature_calculation/01_ranked_essentiality.ipynb                 | Computation of ranked and normalized gene essentiality features                         |
| 02_feature_calculation/02a_IdentifySharedInteractors_BioGRID.ipynb  | Identification of shared protein–protein interactors (BioGRID)                          |
| 02_feature_calculation/02b_ParalogPPIEssentiality_BioGRID.ipynb     | Calculation of the weighted average essentiality of PPI (BioGRID)                       |
| 02_feature_calculation/03a_IdentifySharedInteractors.ipynb          | Identification of shared protein–protein interactors (STRING)                           |
| 02_feature_calculation/03b_ParalogPPIEssentiality.ipynb             | Calculation of the weighted average essentiality of PPI (STRING)                        |
| 02_feature_calculation/03c_ParalogPPIExpression.ipynb               | Calculation of the weighted average expression of PPI (STRING)                          |
| 02_feature_calculation/04_GO_expression.ipynb                       | Calculation of the average gene expression across annotated gene ontology terms         |
| 02_feature_calculation/05_GO_ranked_essentiality.ipynb              | Calculation of the average ranked gene essentiality across annotated gene ontology terms|

| Notebook                                                     | Description                                                    |
|:-------------------------------------------------------------|:---------------------------------------------------------------|
| 03_feature_annotation/01_annotate_DepMap.ipynb               | Annotation of preprocessed genomics features                   |
| 03_feature_annotation/02_annotate_networkfeatures.ipynb      | Annotation of network features                                 |
| 03_feature_annotation/03_annotate_ranked_essentiality.ipynb  | Annotation of ranked and normalized gene essentiality features |
| 03_feature_annotation/04_annotate_scores.ipynb               | Integration of all features and prediction score from Context-Agnostic Classifier |
| 03_feature_annotation/05_annotate_gemini.ipynb               | Integration of GEMINI score and labeling of the training and test dataset |

3. Model Evaluation

| Notebook                                                  | Description                                                                    |
|:----------------------------------------------------------|:-------------------------------------------------------------------------------|
| 04_model_evaluation/01_preprocess_training_dataset.ipynb  | Identification and handling of missing values in the training and test dataset |
| 04_model_evaluation/02_cross_validation.ipynb             | Cross-validation of context-specific random forest classifier                  |

The `scripts/` directory contains a standalone pipeline of the annotation and prediction, generate prediction scores for the provided data

4. Building Breast Cancer Specific Prediction Map

These notebooks construct and analyze a breast cancer cell line specific paralog dependency map

| Notebook                                                                 | Description                                                                         |
|:-------------------------------------------------------------------------|:------------------------------------------------------------------------------------|
| 05_breast_cancer_prediction_map/01_breast_cancer_data_generation.ipynb   | Expand paralog pairs for breast cancer cell lines and annotate them by key mutations|
| 05_breast_cancer_prediction_map/02_breast_cancer_data_annotation.ipynb   | Annotation of breast cancer cell line specific paralog pair dataset with all features |

5. Data Visualization

These notebooks reproduce all figures reported in the manuscript

| Notebook                                | Description                                                                                                         |
|:----------------------------------------|:--------------------------------------------------------------------------------------------------------------------|
| fig2_genomics_feature_analysis.ipynb    | Visualize the predictive performance of genomics-related individual features (Fig. 2)                               |
| fig3_network_feature_analysis.ipynb     | Visualize the predictive performance of network-based individual features (Fig. 3)                                  |
| fig4_cross_validation_visuals.ipynb     | Visualize ROC and PR curves from cross-validation (Fig. 4, Supp. Fig. 2)                                            |
| fig5_evaluate_classifier.ipynb          | Visualize the performance of the context-specific classifier on independent dataset (Fig. 5, Supp. Fig. 4)          |
| fig6_breast_cancer_map.ipynb            | Visualize the distribution of the prediction scores for selected gene pairs among breast cancer cell lines (Fig. 6) |
| sup_fig1_auc_scatter_features.ipynb     | Visualize the consistency of feature performance across independent screens (Supp. Fig. 1)                          |
| sup_fig3_reduced_cross_validation.ipynb | Visualize ROC and PR curves of the cross-validation when the training set is reduced (Supp. Fig. 3)                 |

## Data Sources

Visit [data_sources.md](/data_sources.md) to access the details and download links for all external datasets used in this study

## Environment

The code was developed and tested using a Conda environment with **Python 3.11**

To recreate the environment:
```bash
conda env create -f environment.yml
conda activate paralog_sl_env
```

Notebook execution order follows the directory numbering (01_ → 05_).

### Citation

If you use this code, please cite the paper:  
**"A predicted cancer dependency map for paralog pairs"**  
*(In preparation)*.