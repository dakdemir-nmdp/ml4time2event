# ml4t2event

Machine learning for time to event analysis. This package provides a collection of tools for predicting time-to-event outcomes using various statistical and machine learning methods.

## Installation

To install the package from GitHub, run the following commands in R:

```r
# install.packages("devtools") # Uncomment if devtools is not installed
library(devtools)
install_github("dakdemir-nmdp/ml4t2event")
```

This command will install ml4t2event along with its required dependencies as specified in the DESCRIPTION file. The package depends on R (>= 3.5.0) and utilizes packages such as magrittr, party, prodlim, survival, labelled, and xgboost.

## R Dependencies and Components

The package is designed to work with several R packages:

- **magrittr**: Provides the pipe operator `%>%` for cleaner code.
- **party**: Implements recursive partitioning methods.
- **prodlim**: Offers functions for survival analysis.
- **survival**: Essential for creating survival objects and related analyses.
- **labelled**: Helps manage data labels and dictionaries.
- **xgboost**: Used for gradient boosting methods.

Additionally, the package suggests other packages (energy, gtsummary, data.table, openxlsx, pacman, checkmate, dplyr, stringr, missRanger, randomForestSRC, glmnet, mgcv, partykit, rpart, gbm, mlr3, mlr3learners, pec, pracma, fastcmprsk, haven, BART, forcats, pseudo, recipes, rsample, pbapply) to extend its functionality.

## DESCRIPTION and NAMESPACE

The DESCRIPTION file outlines the package metadata including:
- Package title: Machine learning for time to event analysis
- Version: 0.1.0
- Author and Maintainer: Deniz Akdemir <dakdemir@nmdp.org>
- License: GPL (>= 2)

The NAMESPACE file exposes package functions and imports essential functions from several R packages to facilitate a smooth workflow.

## Usage

Load the package in your R session with:

```r
library(ml4t2event)
```

Then, you can explore various functions for:
- Building and evaluating time-to-event models.
- Creating survival objects.

## Contributing

Contributions, suggestions, and bug reports are welcome. Please use the GitHub repository to submit issues or pull requests.

## License

This project is licensed under the GPL (>= 2) License.

## Acknowledgements

We thank all contributors and users for their support and feedback.

Enjoy using ml4t2event!
