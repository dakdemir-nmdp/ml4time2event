# R Environment Setup Notes

This document outlines the R environment configuration for the ml4time2event project, used for testing and development.

## R Installation
- **Version**: 4.4.1 (2024-06-14) -- "Race for Your Life"
- **Path**: `/usr/local/bin/R`
- **Platform**: aarch64-apple-darwin20 (Apple Silicon macOS)

## Package Management
- **Tool**: renv (version 1.1.1 in the autoloader)
- **Activation**: Automatic when starting R in the project directory (`/Users/dakdemir/Library/CloudStorage/OneDrive-NMDP/Year2025/Github/ml4time2event`)
- **Library Location**: `renv/library/` (local to the project)
- **Snapshot Type**: Explicit (packages are explicitly listed)

## Usage Instructions
- Use `/usr/local/bin/R` or ensure `/usr/local/bin` is in your PATH.
- When working in VS Code, the R terminal should use this setup automatically.
- renv will isolate packages to avoid conflicts with system or other project libraries.
- To check renv status: Run `renv::status()` in an R session within the project directory.

## Additional Notes
- No specific R version is pinned in `renv/settings.json` (r.version: null).
- Bioconductor version is not specified.
- External libraries: None configured.
- PPM (Package Manager) is not explicitly enabled.