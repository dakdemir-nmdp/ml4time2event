# AI Agent Instructions for ml4time2event

This directory contains specialized AI agent instructions for working inside **ml4time2event**, an R package for survival and competing-risk machine learning. Each agent is tuned for the workflows we rely on when preparing the package for CRAN- or GitHub-style releases (roxygen2 docs, `testthat`, pkgdown, renv, etc.).

## Available Agents

### Core Agents

#### [Codebase Finisher](finisher.agent.md)
- **Role**: Stabilize structure, roxygen docs, tests, and packaging so ml4time2event is ready for a minimal tagged release.
- **Expertise**: R package hygiene, `testthat` integration, roxygen2 documentation, pkgdown/publishing checklists.
- **Use when**: The repo needs coordinated cleanup (DESCRIPTION, NAMESPACE, tests, vignettes, CI).

#### [Test Pilot](test-pilot.agent.md)
- **Role**: Add targeted `testthat` coverage for survival/competing-risk helpers without refactoring the package source.
- **Expertise**: `testthat` 3e, fixtures via `withr`, mocking model registries/data helpers, deterministic seeds.
- **Use when**: Need runnable tests for a module or regression scenario.

#### [Statistical Methods Auditor](stats-auditor.agent.md)
- **Role**: Review modelling pipelines (task constructors, ensembling, evaluation metrics) for leakage, resampling mistakes, and survival-method assumptions.
- **Expertise**: Survival analysis, competing risks, resampling, metrics like C-index/Brier/SHAP, reproducibility.
- **Use when**: Validating methodological soundness or deciding what to test/fix before release.

### Specialized Agents

#### [Data Validator](data-validator.agent.md)
- **Role**: Verify task inputs (time/status/cause columns, factor levels, predictors) align with `ml4t2e_task_*()` requirements.
- **Expertise**: Tidy survival data checks, missingness audits, recipe compatibility, caching minimal fixtures.

#### [Docs Doctor](docs-doctor.agent.md)
- **Role**: Improve roxygen docs, README/vignettes, pkgdown articles, and runnable examples.
- **Expertise**: Technical writing for R packages, roxygen skeletons, `ml4t2e_fit()` examples, reproducible seeds.

#### [Survival/Competing-Risk Theory Expert](survival-theory.agent.md)
- **Role**: Verify theoretical correctness of survival/CIF estimators, metrics, and inference logic; compare against authoritative references.
- **Expertise**: Survival analysis theory, competing-risk methodology, derivations, replication against `survival`/`cmprsk`/`riskRegression`.

#### [Performance/Hardware Tuner](perf-tuner.agent.md)
- **Role**: Recommend safe performance knobs (parallelism, registry defaults, recipe caching) for CPU/GPU training workloads.
- **Expertise**: Benchmarks with `bench` or `tictoc`, parallel backends (`future`, `foreach`), tuning survival engines.

#### [Release Manager](release-manager.agent.md)
- **Role**: Handle DESCRIPTION/NAMESPACE/versioning, devtools/pkgdown workflows, CRAN checks, and GitHub Releases for ml4time2event.
- **Expertise**: usethis release process, `devtools::check()`, tarball verification, renv snapshot coordination.

#### [Security Auditor](security-auditor.agent.md)
- **Role**: Review dependency/licensing risks (CRAN/Bioc/GitHub packages), sandboxing for example data, and file I/O surfaces.
- **Expertise**: R dependency graphing, license compatibility, reproducible environment pinning.

## Agent Communication

Agents hand off via the `handoffs` metadata in each `.agent.md`. Typical flows:

- **Finisher → Test Pilot**: Source cleanup exposes areas that still need deterministic `testthat` coverage.
- **Test Pilot → Stats Auditor**: Tests surface methodological doubts needing a deeper audit.
- **Stats Auditor → Finisher**: Findings require code or configuration changes.
- **Any agent → Docs Doctor**: When documentation/examples lag behind the recent changes.

## File Format

Each agent file has:
- **YAML frontmatter** for metadata, tool permissions, and handoff options.
- **Markdown guidance** describing behavior, scope, and deliverables tailored to ml4time2event.

## Usage Guidelines

1. Pick the agent whose scope matches your task (tests, release checklist, docs, etc.).
2. Provide context when invoking an agent (module path, dataset description, release target).
3. Follow the recommended handoffs so specialists can work sequentially without losing context.
4. Keep these docs current when the R package architecture or tooling changes.

## Related Files

- `../copilot-instructions.md`: Copilot-wide context for the ml4time2event package.
- `../workflows/ci.yml`: GitHub Actions workflow the agents assume when asking for CI proof.

---

VS Code recognizes `.agent.md` files as agent definitions. Use the handoff buttons (e.g., Finisher → Test Pilot → Stats Auditor) to keep context flowing between specialists.
