# Agent Guidance for EMC2

This document provides guidance for AI agents interacting with the `EMC2` R package. It covers the high-level architecture, key workflows, and detailed information about specific components.

## EMC2 Package Overview

`EMC2` is an R package for fitting Bayesian hierarchical cognitive models of choice. It supports several models, including the Diffusion Decision Model (DDM), Linear Ballistic Accumulator (LBA), Racing Diffusion Model (RDM), and Lognormal Race (LNR).

The core workflow is divided into four main phases:

### 1. Model Specification and Prior Definition

This phase involves defining the model structure, specifying the experimental design, and setting up the priors.

*   **`make_emc`**: The central function to create an `emc` object. This object encapsulates the model, design, and data.
*   **`make_data`**: Prepares the raw data for use with `EMC2`. It links the data to the model's factor structure.
*   **`get_prior`**: Used to specify and inspect the prior distributions for the model parameters.
*   **`plot.emc.prior`**: Visualizes the prior distributions, which is crucial for understanding the implications of the chosen priors.

### 2. Model Fitting (Sampling)

This phase runs the MCMC sampler to draw samples from the posterior distribution of the model parameters.

*   **`run_emc`**: The main workhorse function that initiates the sampling process. It takes the `emc` object created in the first phase and runs the particle Metropolis MCMC sampler.

### 3. Assessing Convergence

After sampling, it's essential to check whether the MCMC chains have converged to a stable posterior distribution.

*   **`plot.emc`**: A versatile plotting function. When used on the output of `run_emc`, it can generate trace plots of the MCMC chains.
*   **`gd_summary`**: Provides the Gelman-Rubin diagnostic (`Rhat`) for assessing convergence across multiple chains.
*   **`ess_summary`**: Calculates the Effective Sample Size (ESS), indicating how many independent samples the chains are equivalent to.

### 4. Inference and Model Evaluation

Once convergence is confirmed, this phase involves analyzing the posterior distributions, evaluating model fit, and comparing models.

*   **`plot_chains` / `plot_density`**: Functions to visualize the posterior distributions of the parameters.
*   **`plot_fit` / `plot_predict`**: Generate posterior predictive plots to assess how well the model captures the observed data patterns.
*   **`get_BayesFactor`**: Computes Bayes factors for model comparison, allowing for the quantification of evidence in favor of one model over another.

## C++ Level Functions

The R functions in `EMC2` rely on underlying C++ code for performance-critical operations, especially the likelihood calculations. The C++ functions are exposed to R using `Rcpp`. Here is a brief overview of some key C++ components:

*   **Likelihood Functions**: For each cognitive model (LBA, DDM, etc.), there is a corresponding C++ function that calculates the log-likelihood for a given set of parameters and data. These are the core of the model fitting process. For example, `dlba_c` is the likelihood function for the LBA model.
*   **Particle Likelihood**: The package uses a particle Metropolis MCMC sampler. The `particle_ll` C++ function is a key part of this, calculating the likelihood for a single participant's data, which is then used within the MCMC sampler.
*   **Random Effects**: The hierarchical nature of the models is handled by C++ functions that manage the random effects (i.e., the individual differences). These functions deal with the multivariate normal distributions for the random effects and their covariance matrices.
*   **Solvers and Numerical Methods**: Some models, like the DDM, do not have a simple closed-form likelihood function. The package includes C++ code to solve the underlying stochastic differential equations (SDEs) or to use other numerical methods to approximate the likelihood.

## Sampling Architecture and Likelihood Calculation

The core of the package is its sampling engine, which uses a Particle Metropolis with Gibbs (PMwG) algorithm. This is a hybrid MCMC strategy well-suited for hierarchical models.

### 1. High-Level Sampling Process (R)

The fitting process is managed by the `run_emc` function in `R/fitting.R`, which orchestrates the sampler through several stages: "preburn", "burn", "adapt", and "sample". The core logic resides in `run_stage` (in `R/sampling.R`), which executes the main MCMC loop for a given number of iterations.

Each iteration of the loop consists of two main steps:

1.  **Gibbs Step (`gibbs_step`)**: The group-level (hyper) parameters (`theta_mu` and `theta_var`) are sampled directly from their full conditional distributions. This is possible because of the choice of conjugate priors and is highly efficient.

2.  **Particle Metropolis Step (`new_particle`)**: The individual-level parameters (`alpha`, one set for each subject) are updated using a Metropolis-Hastings step. Because the likelihood for a subject's parameters depends on the group-level parameters in a complex way, this step is more involved. This step is parallelized across subjects.

### 2. The Particle Metropolis Step in Detail

The `new_particle` function in `R/sampling.R` is the heart of the individual-level parameter update. For each subject, it generates a set of candidate parameter vectors (particles) and then uses their likelihoods to decide on the next sample.

*   **Mixture of Proposals**: To ensure robust and efficient exploration of the parameter space, new particles are drawn from a mixture of different proposal distributions. The mixture adapts across the sampling stages:
    *   **Group-level (Prior)**: Particles are drawn from the group-level distribution (`theta_mu`, `theta_var`). This encourages global exploration.
    *   **Previous Sample**: Particles are drawn from a distribution centered on the subject's previous parameter values. This encourages local exploration.
    *   **Empirical Covariance**: In later stages, proposals are also drawn from a distribution based on the empirical covariance of the subject's posterior samples so far (`chains_var`).
    *   **Efficient Proposals**: In the final "sample" stage, an "efficient" proposal distribution is added, which is based on a conditional approximation of the posterior (`eff_mu`, `eff_var`).

*   **Weighting and Resampling**: The log-likelihood for each particle is calculated by the C++ function `calc_ll`. These log-likelihoods are combined with the proposal and prior densities to calculate importance weights. A new parameter vector for the subject is then sampled from the particles, weighted by these importance weights.

### 3. Likelihood Calculation Engine (C++)

The most computationally intensive part is calculating the log-likelihood for a set of particles. This is handled in C++ for speed, primarily in `src/particle_ll.cpp`.

*   **Entry Point**: The `calc_ll` function is called from R. It receives a matrix of particles (`p_matrix`), where each row is a parameter vector for one particle.

*   **Main Loop**: The function iterates through each particle in `p_matrix`. For each particle, it performs two major steps:
    1.  **Parameter Expansion (`get_pars_matrix_rcpp`)**: A single particle contains a compact vector of parameters. This function expands that vector into a full matrix of parameters, one row for each unique trial condition. This involves applying pre-transformations (e.g., `exp` to enforce positivity), multiplying by the design matrix, and applying final transformations.
    2.  **Likelihood Dispatching**: The expanded parameter matrix is passed to a model-specific C++ function based on a `type` string (e.g., "DDM", "LBA").

*   **Race Model Engine (`c_log_likelihood_race_cens_trunc`)**: This is a highly generic and powerful function for all race models (LBA, RDM, LNR).
    *   **Function Pointers**: It uses C++ function pointers (`RacePdfFun`, `RaceCdfFun`) to call the appropriate probability density and cumulative distribution functions for the specific model (e.g., `lba_dfun_adapter`).
    *   **Batching**: It intelligently separates trials with standard, finite response times from trials that are censored or fall outside truncation bounds. The finite RTs are processed in a single, fast, vectorized operation.
    *   **Numerical Integration**: For censored or truncated trials, where the likelihood is an integral over a time interval, it uses the GNU Scientific Library (GSL) to perform numerical integration (`integrate_for_kth_winner_cpp`).
    *   **Caching**: To avoid re-computing the expensive normalization integral for truncated distributions, it caches the results in a `std::unordered_map` (`inv_Z_cache`).

### 4. Potential Optimizations

While the sampler is already highly optimized, several areas could be explored for further performance gains:

*   **Parameter Mapping (`get_pars_matrix`)**: This expansion process is executed for every particle at every iteration. If the same parameter vector is proposed multiple times (within or across subjects), the result of this expensive mapping could be cached.
*   **Truncation Cache Key**: The cache key for the truncation integral is currently a string generated by serializing the relevant parameters. This is robust but incurs overhead. Using a custom C++ struct as the key with a specialized hash function could offer a speed-up.
*   **GSL Integration Workspace**: The `integrate_for_kth_winner_cpp` function allocates and frees a GSL integration workspace for every call. For a given particle, this function might be called multiple times. Reusing a single workspace object per particle could reduce allocation overhead.
*   **Redundant Target Model Logic**: The likelihood functions for the specialized redundant target models (`c_log_likelihood_redundant_target_race` and its variants) are very complex with many conditional branches. While correct, they could be a target for refactoring to improve clarity and potentially identify shared computations.
*   **Gauss-Hermite Quadrature**: The `..._substitution` model uses a fixed 5-point Gauss-Hermite quadrature for integration. Making the number of nodes configurable would allow a user to trade off between accuracy and speed.

## Model Design Architecture

The process of specifying a model's design in EMC2 involves creating a blueprint that maps experimental factors to model parameters. This is primarily handled by functions in `R/design.R` and `R/make_data.R`.

### 1. The Design Blueprint (`design()`)

The user's journey begins with the `design()` function. This function does not build the final numerical matrices itself, but rather collects all the necessary specifications into a comprehensive `emc.design` object.

*   **Formula-Based Specification**: The core of the design is a list of R formulas (e.g., `list(v ~ S*E, a ~ E)`). Each formula defines a linear model for a specific cognitive model parameter (like `v` or `a`) based on experimental factors (`S`, `E`).
*   **Components**: The `design` object also stores information about experimental `factors`, response levels (`Rlevels`), custom `contrasts` for categorical variables, and `constants`.
*   **Output**: The result is a list object that acts as a self-contained blueprint for how to construct the design for any given dataset that conforms to the specified factors.

### 2. Realizing the Design (`design_model()`)

This internal function is the workhorse that translates the `emc.design` blueprint into the concrete data and matrices needed by the sampler. It is called by `make_emc()` before sampling begins.

1.  **Accumulator Expansion (`add_accumulators`)**: For race models (LBA, RDM, LNR), the data must be structured in a "one row per accumulator per trial" format. This function takes the original data and expands it, creating `n_accumulators` rows for each original trial. It adds an `lR` (latent response) column to identify each accumulator.

2.  **Design Matrix Creation (`make_dm`)**: For each parameter (e.g., `v`, `a`, `t0`), this function calls R's standard `model.matrix` utility. It uses the corresponding formula and the (expanded) data frame to generate a numerical design matrix. Each column in this matrix corresponds to a parameter to be estimated (e.g., `v`, `v_S2`, `v_E2`, `v_S2:E2`).

3.  **Data Compression (`compress_dadm`)**: This is a critical performance optimization. Since many trials in an experiment can be identical in terms of their design, calculating the likelihood for each one would be redundant. This function identifies unique trial configurations and creates a compressed data frame (`dadm`).
    *   **The "Hack"**: Uniqueness is determined by creating a string key for each row by `paste`-ing together all the values from all design matrices for that row, plus subject, RT, and response information. While effective, this is not maximally efficient.
    *   **The `expand` Attribute**: The function stores an integer vector, `expand`, as an attribute of the compressed data. This vector serves as an index to map the results from the unique (compressed) trials back to the full, original trial structure.

### 3. Group-Level Design (`group_design()`)

For hierarchical models, `group_design()` allows for specifying linear models for the group-level means. For example, it can model how the mean of a subject-level parameter (e.g., `v_S1`) varies as a function of a subject-level covariate (e.g., `age`). It also produces design matrices that are used during the Gibbs step of the sampler.

### 4. Potential Optimizations and Enhancements

The design architecture is powerful but has areas that could be enhanced for better performance, clarity, and robustness.

*   **Principled Data Compression**: The `paste`-based compression in `compress_dadm` could be replaced with more modern, efficient, and robust methods.
    *   **Suggestion**: Use `data.table` or `dplyr::group_by` followed by `.gid`. These tools are written in C and are exceptionally fast and memory-efficient at finding unique group indices, which is exactly what `compress_dadm` is doing.
    *   **Impact**: High. This would make the pre-computation phase significantly faster, especially for large datasets. It is a non-breaking change as long as the final `expand` attribute is identical.

*   **Use of Sparse Matrices**: Design matrices, especially after dummy-coding factors, are often very sparse (mostly zeros). Currently, they are stored as standard dense matrices.
    *   **Suggestion**: Modify `make_dm` to return sparse matrices from the `Matrix` package (e.g., `dgCMatrix`). The C++ backend (`RcppArmadillo`) has excellent support for sparse matrices (`sp_mat`), so the matrix multiplication in `c_map_p` would become much more efficient.
    *   **Impact**: High. This would dramatically reduce the memory footprint of the `emc` object and speed up the parameter expansion step in the C++ code. This is a key, non-breaking enhancement.

*   **Refactor Data Expansion**: The `add_accumulators` function physically replicates the data frame, which can cause a large memory spike. 
    *   **Suggestion**: A more advanced approach would be to avoid this replication. The C++ likelihood function could be re-written to understand the race model structure and implicitly calculate the contribution of each accumulator without needing the data to be in the expanded long format. This is a more significant, breaking change but would offer the largest memory savings.

*   **Improved Validation and User Experience**: The `design` function could provide more upfront validation.
    *   **Suggestion**: Add checks to `design()` to ensure that all factors used in formulas are present in the `factors` list, that contrast names match factor names, and that constants do not conflict with estimated parameters. Providing clearer, more targeted error messages would improve usability.
    *   **Impact**: Medium. This would improve user experience and prevent common errors, without changing the core logic.