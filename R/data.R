## Prepared datasets for vignettes ------------------------------------------------

#' Prepared Lung Cancer Survival Data
#'
#' A cleaned and imputed version of the NCCTG lung cancer study used in the
#' survival-analysis vignette. The dataset has missing values imputed with
#' `missRanger`, categorical variables encoded as factors, and includes a few
#' engineered predictors that illustrate feature creation.
#'
#' @format A data frame with 228 rows and 11 variables:
#' \describe{
#'   \item{time}{Follow-up time in days.}
#'   \item{status}{Event indicator (1 = death, 0 = censored).}
#'   \item{age}{Age in years at registration.}
#'   \item{sex}{Factor indicating sex (`Male`, `Female`).}
#'   \item{ph.ecog}{Eastern Cooperative Oncology Group performance score (factor).}
#'   \item{ph.karno}{Karnofsky performance score assessed by physician.}
#'   \item{pat.karno}{Karnofsky performance score assessed by patient.}
#'   \item{meal.cal}{Calories consumed per day.}
#'   \item{wt.loss}{Weight loss in the last six months (kg).}
#'   \item{age_group}{Factor with age groups (`Young`, `Middle`, `Elderly`).}
#'   \item{performance_good}{Indicator (1/0) for physician Karnofsky score >= 80.}
#' }
#' @details Derived from the NCCTG lung cancer dataset included with the package.
#'   Missing values were imputed using `ImputeMissinRecordsData()` and engineered
#'   features were added for demonstration purposes.
#' @source North Central Cancer Treatment Group (NCCTG) clinical trial.
"lung_survival_prepared"


#' Prepared Bone Marrow Transplant Competing Risks Data
#'
#' A preprocessed version of the bone marrow transplant dataset used in the
#' competing-risks vignette. The data include factor encodings of the categorical
#' predictors and retain the original competing-risks event coding.
#'
#' @format A data frame with 177 rows and 7 variables:
#' \describe{
#'   \item{sex}{Factor indicating donor sex (`F`, `M`).}
#'   \item{d}{Factor indicating diagnosis (`ALL`, `AML`).}
#'   \item{phase}{Factor for disease phase at transplant (`CR1`, `CR2`, `CR3`, `Relapse`).}
#'   \item{age}{Age (years) at transplant.}
#'   \item{status}{Event indicator (0 = censored, 1 = relapse, 2 = treatment-related mortality).}
#'   \item{source}{Factor for stem-cell source (`BM+PB`, `PB`).}
#'   \item{ftime}{Follow-up time in months.}
#' }
#' @details Derived from the bone marrow transplant competing risks dataset
#'   distributed with the package and prepared for direct use with
#'   `RunCRModels()`.
#' @source Keiding et al. (1996) bone marrow transplant study.
"bmt_competing_risks_prepared"
