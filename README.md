# A predicted cancer dependency map for paralog pairs #

This repository contains the scripts and Jupyter notebooks for the research paper titled:
**"A predicted cancer dependency map for paralog pairs"**, currently available in XXX

The repository includes:
- `notebooks/01_preprocessing/` - Preprocess the input files 
- `notebooks/02_feature_calculation/` - Compute network based features
- `notebooks/03_feature_annotation/`- Map calculated features to training and test datasets
- `notebooks/04_model_evaluation/` - Train and evaluate the context-specific Random Forest classifier
- `notebooks/05_visualization/` - Visualize the performance of individual features and the classifier
- `scripts/` - Python scripts to calculate the prediction score

## Overview of the Notebooks

### Data Preprocessing

All the input data were processed in these notebooks

| Notebook                       | Description                                                                                  |
|:-------------------------------|:---------------------------------------------------------------------------------------------|
| preprocess_ito.ipynb           | Process LFC and FDR datasets in the CRISPR screen by Ito et al.|
| preprocess_klingbeil.ipynb     | Process LFC and FDR datasets in the CRISPR screen by Klingbeil et al. |
| preprocess_parrish.ipynb       | Process LFC dataset in the CRISPR screen by Parrish et al. |
| preprocess_DepMap_22Q4.ipynb   | Process the transcriptomics, gene effect scores, etc datasets from DepMap22Q4 release|   

### Feature Calculation and Annotation

All the features calculated and mapped using the notebooks below

| Notebook                                                 | Description                                                                 |
|:---------------------------------------------------------|:----------------------------------------------------------------------------|
| 02_feature_calculation/01_ranked_essentiality.ipynb      | Compute ranked gene essentiality feature     |
| 02_feature_calculation/02a_IdentifySharedInteractors_BioGRID.ipynb  | Calculate shared protein interactors of the paralog pairs (BioGRID)    |
| 02_feature_calculation/02b_ParalogPPIEssentiality_BioGRID.ipynb  | Calculate the weighted average essentiality of PPI (BioGRID)              |
| 02_feature_calculation/03a_IdentifySharedInteractors.ipynb        | Calculate shared protein interactors of the paralog pairs (STRING)       |
| 02_feature_calculation/03b_ParalogPPIEssentiality.ipynb     | Calculate the weighted average essentiality of PPI (STRING)                    |
| 02_feature_calculation/03c_ParalogPPIExpression.ipynb     | Calculate the weighted average expression of PPI (STRING)                        |
| 02_feature_calculation/04_GO_expression.ipynb             | Calculate the average gene expression for the annotated gene ontology terms      |
| 02_feature_calculation/05_GO_ranked_essentiality.ipynb     | Calculate the average gene essentiality for the annotated gene ontology terms   |

| Notebook                              | Description                                                                 |
|:--------------------------------------|:----------------------------------------------------------------------------|
| 03_feature_annotation/01_annotate_DepMap.ipynb             | Map preprocessed genomics features to the desired format for downstream analysis. |
| 03_feature_annotation/02_annotate_networkfeatures.ipynb          | Map calculated network features to the desired format for downstream analysis.      |
| 03_feature_annotation/03_annotate_ranked_essentiality.ipynb | Map ranked and normalized gene essentiality features to the desired format for downstream analysis.|
| 03_feature_annotation/04_annotate_scores.ipynb  | Map all features and prediction score from Context-Agnostic Classifier in the desired format for downstream analysis.|
| 03_feature_annotation/05_annotate_gemini.ipynb             | Map GEMINI score and label the training and test dataset in the desired format for downstream analysis.      |

### Model Evaluation

| Notebook                              | Description                                                                 |
|:--------------------------------------|:----------------------------------------------------------------------------|
| 04_model_evaluation/01_preprocess_training_dataset.ipynb            | Identify missing values and fill/drop them in the training and test dataset |
| 04_model_evaluation/02_cross_validation.ipynb            | Run the cross-validation of context-specific random forest classifier on the training dataset |

### Data Analysis

| Notebook                                     | Description                                                                                   |
|:---------------------------------------------|:----------------------------------------------------------------------------------------------|
| fig2_genomics_feature_analysis.ipynb         | Visualize the predictive performance of genomics-related individual features (Fig. 2)         |
| fig3_network_feature_analysis.ipynb          | Visualize the predictive performance of network-based individual features (Fig. 3)            |
| fig4_cross_validation_visuals.ipynb          | Visualize ROC and PR curves of the cross-validation (Fig. 4)                                  |
| fig5_evaluate_classifier.ipynb               | Visualize the performance of the context-specific classifier on independent dataset (Fig. 5)  |
| fig6_breast_cancer_map.ipynb                 | Visualize the distribution of the prediction scores for selected gene pairs (Fig. 6)          |

## Data Sources

Visit [data_sources.md](/data_sources.md) to access the input data.

## Environment

The code was developed and tested using a Conda environment with **Python 3.11**

To recreate the environment:
```bash
conda env create -f environment.yml
conda activate paralog_sl_env
```

### Citation

If you use this code, please cite the paper:  
**"A predicted cancer dependency map for paralog pairs"**  
*(In preparation)*.