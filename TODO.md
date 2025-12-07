# ml4time2event Roadmap

## 🚀 Immediate Priorities: Release Readiness & Stability
- [x] **Package Health**: Ensure `devtools::check()` passes with 0 errors/warnings.
- [x] **Dependencies**: Move heavy dependencies (e.g., specific engines) to `Suggests` to ensure a lightweight default install. Implement `ml4t2e_install_extras()` helper if needed.
- [x] **Documentation**:
    - Finalize README with a clear "Quick Start".
    - [x] Build and deploy `pkgdown` site.
    - [x] Ensure all Vignettes (Survival, Conformal, SHAP) render correctly and reproducible.
- [x] **Testing**:
    - Add integration tests for full pipeline (save/load cycles).
    - Verify Conformal Prediction implementation on varied datasets.

## ✨ Core Feature Polish
- [x] **Pipeline Persistence**: Implement/Verify robust `ml4t2e_save()` / `ml4t2e_load()` that validates recipes and model state, ensuring predictions on new data work seamlessly.
- [x] **Explainability**: Ensure SHAP/Explainability tools are consistent with the latest pipeline API and handle new model types correctly.

## 🔜 Next Steps (High Impact)
- [x] **Benchmarking/Leaderboard**: Implement a unified way to view model performance comparisons (C-index, Brier score) from a fitted pipeline.
- [x] **Unified Fit Interface**: Refine `ml4t2e_fit` to better handle resource budgets (time, max models) and auto-ensembling defaults.
