# Blueprint Document for `ml4time2event` Codebase

## Overview
The `ml4time2event` repository is a comprehensive R-based framework for survival analysis and competing risks modeling. It includes utilities for data preparation, model training, evaluation, and visualization.

## Repository Structure

### Root Directory
- **DESCRIPTION**: Metadata about the R package.
- **LICENSE**: Licensing information.
- **ml4time2event.Rproj**: RStudio project file.
- **NAMESPACE**: Exported functions and imports.
- **README.md**: Overview and usage instructions.

### `data/`
- **bmt_competing_risks_prepared.rda**: Preprocessed dataset for competing risks analysis.
- **lung_survival_prepared.rda**: Preprocessed dataset for survival analysis.

### `data-raw/`
- **prepare_vignette_datasets.R**: Script to prepare datasets for vignettes.

### `inst/`
- **examples/**: Example scripts.
- **extdata/**: External data files.

### `man/`
Contains `.Rd` files documenting functions, datasets, and utilities.

### `R/`
Contains R scripts implementing core functionality. Key components include:
- Data preparation utilities.
- Model definitions (e.g., Cox, Fine-Gray, BART, etc.).
- Ensemble methods.
- Evaluation metrics (e.g., Brier Score, Concordance Index).
- Visualization tools (e.g., survival and CIF curves).

### `tests/`
Unit tests for validating functionality.

### `vignettes/`
R Markdown documents for detailed examples and tutorials.

## Key Functionalities

### Data Preparation
- **`bake_data_recipe`**: Applies preprocessing steps to datasets.
- **`minimal_data_recipe`**: Minimal preprocessing pipeline.
- **`CollapseRareCategoriestoOtherData`**: Collapses rare categories in categorical variables.

### Models
- **Competing Risks Models**:
  - `CRModel_Cox`
  - `CRModel_FineGray`
  - `CRModel_RF`
  - `CRModel_xgboost`
- **Survival Models**:
  - `SurvModel_Cox`
  - `SurvModel_RF`
  - `SurvModel_xgboost`

### Evaluation
- **`BrierScore`**: Calculates Brier Score for survival models.
- **`integratedBrier`**: Computes integrated Brier Score.
- **`integratedC`**: Calculates integrated Concordance Index.

### Visualization
- **`plot_survival_curves`**: Plots Kaplan-Meier survival curves.
- **`plot_cif_curves`**: Plots cumulative incidence functions.
- **`plot_model_performance`**: Visualizes model performance metrics.

## Development Workflow
1. **Data Preparation**: Use scripts in `data-raw/` to preprocess datasets.
2. **Model Training**: Leverage functions in `R/` to train models.
3. **Evaluation**: Use evaluation metrics to assess model performance.
4. **Visualization**: Generate plots for analysis and reporting.

## Dependencies
- R (installed at `/usr/local/bin/R`)
- Required R packages (specified in `DESCRIPTION`):
  - `survival`
  - `xgboost`
  - `randomForest`
  - `ggplot2`

## Testing
- Unit tests are located in `tests/`.
- Run tests using `devtools::test()`.

## Documentation
- Function documentation is in `man/`.
- Tutorials and examples are in `vignettes/`.

## Contribution Guidelines
1. Follow the structure and naming conventions.
2. Document all new functions in `.Rd` format.
3. Add unit tests for new features.

## Future Enhancements
- Expand model support (e.g., neural networks).
- Add more evaluation metrics.
- Improve visualization capabilities.