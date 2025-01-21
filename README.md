# Predicting context-specific synthetic lethality between paralogs in cancer #

This repository contains the codebase and Jupyter notebooks for the research paper titled:
**"Predicting context-specific synthetic lethality between paralogs in cancer"**, currently available in 

The repository includes:
- Preprocessing of input datasets 
- Calculation and mapping of features
- Training of the Random Forest classifier
- Visualization of the performance of individual features and the classifier

## Data Processing Notebooks

### Preprocessing of input datasets

All the input data were processed in these notebooks

| Notebook                       | Description                                                                                  |
|:-------------------------------|:---------------------------------------------------------------------------------------------|
| preprocess_ito.ipynb           | Process gene pairs screened by Ito et al.|
| preprocess_klingbeil.ipynb     | Process gene pairs screened by Klingbeil et al. |  
| preprocess_parrish.ipynb       | Process gene pairs screened by Parrish et al. |
| preprocess_dede.ipynb          | Process gene pairs screened by Dede et al. |
| preprocess_chymera.ipynb       | Process gene pairs screened by Gonatopoulos-Pournatzis et al. |
| preprocess_thompson.ipynb      | Process gene pairs screened by Thompson et al. |

### Calculation and mapping of features

All the features calculated and mapped using the notebooks below

| Notebook                          | Description                                                                 |
|:----------------------------------|:----------------------------------------------------------------------------|
| map_depmap_data.ipynb             | Map calculated features to the desired format for downstream analysis.      |
| map_depmap_data.ipynb             | Map calculated features to the desired format for downstream analysis.      |
| map_depmap_data.ipynb             | Map calculated features to the desired format for downstream analysis.      |
| map_depmap_data.ipynb             | Map calculated features to the desired format for downstream analysis.      |
| map_depmap_data.ipynb             | Map calculated features to the desired format for downstream analysis.      |
| map_ranked_essentiality.ipynb     | Compute features (e.g., gene expression, mutation status) for Cell Line 1. |
| map_ppi_data.ipynb                | Compute features (e.g., gene expression, mutation status) for Cell Line 1. |
| map_go_data.ipynb                 | Compute features (e.g., gene expression, mutation status) for Cell Line 1. |
| map_sequence_identity.ipynb       | Compute features (e.g., gene expression, mutation status) for Cell Line 1. |
| map_gemini_score.ipynb            | Compute features (e.g., gene expression, mutation status) for Cell Line 1. |

## Data Analysis Notebooks

| Notebook                                    | Figures        | Description           |
|:--------------------------------------------|:---------------|:----------------------|
|fig1_genomics_feature_analysis.ipynb         | Figure 1       | Visualize the performance of individual features |
|fig2_network_feature_analysis.ipynb          | Figure 2       | Visualize the performance of individual features |
|fig3_train_classifier.ipynb                  | Figure 3       | Visualize the performance metrics of the context-specific classifier (e.g., ROC, PR curves) |
|fig4_evaluate_classifier.ipynb               | Figure 4       | Visualize the performance of the context-specific classifier on independent datasets |

### Citation

If you use this code, please cite the paper:  
**"Predicting context-specific synthetic lethality between paralogs in cancer"**  
*(In preparation)*.

### Data Sources

Visit [data_sources.md](/data_sources.md) to access the input data.
