# lsirm12pl 2.0.0 (2026-01-29)

## Major Changes

* **NEW**: Added Graded Response Model (GRM) for ordinal/Likert-scale data
  - `lsirmgrm()`: 1PL GRM with latent space
  - `lsirmgrm2pl()`: 2PL GRM with item discrimination parameters
  - Supports complete data, MAR, and MCAR missing data mechanisms
  - Implements threshold ordering constraint (β₁ > β₂ > ... > βₖ₋₁) for identifiability
  - Reference: De Carolis, Kang, & Jeon (2025), DOI: 10.1080/00273171.2025.2605678

* **NEW**: Adaptive MCMC algorithm
  - Automatic tuning of proposal standard deviations during burn-in
  - Parameter-specific target acceptance rates
  - Improved mixing and convergence for complex models
  - Available for all model types (1PL, 2PL, GRM)

* **IMPROVED**: Code refactoring and optimization
  - Reorganized C++ codebase for better maintainability
  - Enhanced numerical stability in likelihood calculations
  - Improved multi-chain diagnostics support

## Bug Fixes

* Fixed acceptance ratio calculation to exclude burn-in period when adaptive MCMC is enabled
* Fixed column name ordering in GRM beta matrix output (column names now match data storage order)
* Improved numerical stability in GRM log-likelihood calculation (eps: 1e-16 → 1e-10)
* Fixed multi-chain diagnostic issues for GRM models (proper beta matrix structure)

## Documentation

* Updated all function documentation with detailed descriptions of adaptive MCMC
* Added comprehensive examples for GRM models
* Clarified threshold ordering constraints in GRM documentation
* Improved documentation for missing data handling (MAR/MCAR)

---

# lsirm12pl 1.3.9

* Last stable release before major 2.0.0 update
* Supported 1PL and 2PL models for binary and continuous data
* Basic spike-and-slab priors for model selection
