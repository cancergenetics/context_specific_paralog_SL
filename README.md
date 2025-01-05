# Predicting context-specific synthetic lethality between paralogs in cancer #

This repository contains the codebase and Jupyter notebooks for the research paper titled:
**"Predicting context-specific synthetic lethality between paralogs in cancer"**, currently available in 

The repository includes:
- Preprocessing of input datasets 
- Calculation and mapping of features
- Training of the Random Forest classifier
- Visualization of the performance of individual features and the classifier

## Preprocessing of Input Datasets

All the input data were processed in these notebooks

| Notebook                       | Description                                                                                  |
|:-------------------------------|:---------------------------------------------------------------------------------------------|
| preprocess_ito.ipynb           | Preprocess data for Cell Line 1, including filtering and normalization.                      |
| preprocess_parrish.ipynb       | Preprocess data for Cell Line 2, handling missing values and scaling.                        |
| preprocess_dede.ipynb          | Combine preprocessed datasets from individual cell lines for further feature calculation.    |
| preprocess_chymera.ipynb       | Combine preprocessed datasets from individual cell lines for further feature calculation.    |
| preprocess_thompson.ipynb      | Combine preprocessed datasets from individual cell lines for further feature calculation.    |

## Calculation and Mapping of Features

All the features calculated and mapped using the notebooks below

| Notebook                          | Description                                                                 |
|:----------------------------------|:----------------------------------------------------------------------------|
| map_depmap_features.ipynb         | Map calculated features to the desired format for downstream analysis.      |
| calculate_ppi.ipynb               | Compute features (e.g., gene expression, mutation status) for Cell Line 1. |
| calculate_go_terms.ipynb          | Compute features (e.g., gene expression, mutation status) for Cell Line 1. |

## Training of the Random Forest Classifier

| Notebook                     | Description                                                 |
|:-----------------------------|:------------------------------------------------------------|
| build_random_forest.ipynb    | Train and evaluate the Random Forest classifier on datasets. |
| independent_predictions.ipynb    | Train and evaluate the Random Forest classifier on datasets. |

## Data Analysis Notebooks

| Notebook                                    | Figures        | Description           |
|:--------------------------------------------|:---------------|:----------------------|
|1_proteomic_replicate_correlation_analysis   | Figure 1       | Visualize the performance of individual features |
|2_compare_replicate-cor_with_mRNA-protein-   | Figure 2       | Visualize the performance of individual features |
|3_robust_rank_aggregation                    | Figure 3       | Visualize the performance metrics of the context-specific classifier (e.g., ROC, PR curves) |
|4_aggregate_ranks_comparison                 | Figure 4       | Compare the performance of  |

### Citation

If you use this code, please cite the paper:  
**"Predicting context-specific synthetic lethality between paralogs in cancer"**  
*(In preparation)*.

### Data Sources

Visit [data_sources.md](/data_sources.md) to access the input data.
