# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build/Test Commands
- Load package during development: `devtools::load_all()`
- Run all tests: `devtools::test()` or `R CMD check`
- Run single test: `testthat::test_file("tests/testthat/path/to/test_file.R")`
- Build package: `devtools::build()` or `R CMD build`
- Linting: `lintr::lint_package()`

## Code Style Guidelines
- **Naming Conventions**:
  - Function names: CamelCase (e.g., `SurvModel_Cox`, `PredictSurvModels`)
  - File names: snake_case (e.g., `surv_cox.R`, `data_cleaning.R`)
  - Variables: snake_case for internal variables
- **Documentation**: Use Roxygen2 style with `@title`, `@description`, `@param`, `@return`
- **Imports**: Use `@importFrom package function` in Roxygen blocks
- **Error Handling**: Validate inputs at function start
- **Assignments**: Use `<-` instead of `=` for variable assignment
- **Formatting**: Indent with 2 spaces, max line length ~80 characters
- **Functions**: Separate arguments on new lines for long argument lists
- **Return Values**: Return named lists for consistent function outputs

File organization follows functional grouping (survival models, competing risks models, utilities, data processing).