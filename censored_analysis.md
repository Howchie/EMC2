# Deep Analysis of Censored RT Pathway in EMC2

## 1. Introduction

This document analyzes the implementation of censored and truncated response time (RT) handling in the EMC2 codebase, specifically tracing the pathway from R-level design functions to the low-level C++ likelihood calculations. It highlights the data structures, mathematical logic, and potential points of brittleness to assist in extracting this logic for a Pull Request.

## 2. Data Pathway: R to C++

### 2.1. User Input (`design()` in `R/design.R`)
The user defines censoring and truncation boundaries either as arguments to `design()` or as columns in the data.
*   **Arguments**: `LT` (Lower Truncation), `UT` (Upper Truncation), `LC` (Lower Censoring), `UC` (Upper Censoring).
*   **Data Columns**: If present in `data`, these columns override the arguments.
*   **Defaults**: `LT=0`, `UT=Inf`, `LC=0`, `UC=Inf`.

### 2.2. Data Preparation (`design_model()` in `R/design.R`)
The `design_model` function constructs the Data Augmented Design Matrix (`dadm`).
*   It ensures `LT`, `UT`, `LC`, `UC` exist as columns in `dadm`.
*   If scalar values were provided, they are replicated to match the number of trials.

### 2.3. Data Compression (`compress_dadm()` in `R/design.R`)
This step groups identical trials to optimize likelihood calculation.
*   **Mechanism**: Creates a unique key by pasting `designs` (parameter mappings), `subjects`, `R` (response), `lR` (accumulator role), and `rt`.
*   **Brittleness Risk**: The uniqueness key **does not explicitly include** `LT`, `UT`, `LC`, or `UC` columns (unless they are part of `Fcov` or `Ffun`).
    *   **Implication**: If trials have identical parameters and RTs (e.g., multiple trials censored at `-Inf`) but *different* truncation/censoring bounds, they may be incorrectly grouped. The C++ code will then use the bounds of the first trial in the group for all of them.

### 2.4. Data Simulation (`make_data()` in `R/make_data.R`)
Used for generating synthetic data.
*   **`make_missing()`**: Replaces `rt` with `-Inf` (left censored) or `Inf` (right censored) based on bounds.
*   **Flags**: `LCdirection` (sets `-Inf`), `UCdirection` (sets `Inf`).

## 3. C++ Implementation (`src/particle_ll.cpp`)

The core logic resides in `src/particle_ll.cpp` (implementation) and `src/utils.h` (definitions).

### 3.1. Entry Point
The function `calc_ll` is exported to R. It dispatches to model-specific likelihood functions. For race models (LBA, RDM, LNR), it calls:
*   `c_log_likelihood_race_cens_trunc`

### 3.2. Likelihood Function (`c_log_likelihood_race_cens_trunc`)
This function iterates over unique trial groups. It reads `LT`, `UT`, `LC`, `UC` from the `dadm` columns (using the index of the first row in the group).

#### Batching Strategy
The function splits trials into two processing paths:
1.  **Finite RT Batch**: Trials with finite, positive RTs within truncation bounds (`LT < RT < UT`) and a known winner. These are processed using vectorized operations for speed.
2.  **Censored/Complex Cases**: Trials with `-Inf`, `Inf`, `NA` (missing), or finite RTs outside standard bounds. These are processed individually using numerical integration.

### 3.3. Mathematical Analysis

The likelihood calculation uses the GSL (GNU Scientific Library) for numerical integration.

**Definitions:**
*   $f_k(t)$: PDF of the $k$-th accumulator (winner) at time $t$.
*   $S_j(t)$: Survivor function ($1 - F_j(t)$) of the $j$-th accumulator (loser).
*   $L(t, k) = f_k(t) \prod_{j \ne k} S_j(t)$: Joint likelihood density of winner $k$ at time $t$.

**Integration Logic (`integrate_for_kth_winner_cpp`):**
This function integrates the joint likelihood density over a specified interval $[A, B]$.
$$ I(k, A, B) = \int_{A}^{B} f_k(t) \prod_{j \ne k} S_j(t) \, dt $$

**Cases:**

1.  **Fast Censoring (Left)**: `rt == -Inf` (encoded as `R_NegInf` in C++)
    *   **Domain**: $[LT, LC]$
    *   **Likelihood**: $\sum_{k} I(k, LT, LC)$ (if winner unknown) or $I(k_{obs}, LT, LC)$ (if winner known).

2.  **Slow Censoring (Right)**: `rt == Inf` (encoded as `R_PosInf` in C++)
    *   **Domain**: $[UC, UT]$
    *   **Likelihood**: $\sum_{k} I(k, UC, UT)$ (if winner unknown) or $I(k_{obs}, UC, UT)$ (if winner known).
    *   *Optimization*: For single accumulator models, this simplifies to $S(UC)$ (survivor at lower bound).

3.  **Missing Data**: `rt == NA`
    *   **Domain**: The unobserved region is treated as the union of the left and right censored regions.
    *   **Likelihood**: $\sum_{k} (I(k, LT, LC) + I(k, UC, UT))$.

**Truncation Normalization:**
To condition the probability on the fact that the RT falls within $[LT, UT]$, the raw likelihood is divided by the normalization constant $Z$.
$$ Z = \sum_{k} \int_{LT}^{UT} f_k(t) \prod_{j \ne k} S_j(t) \, dt $$
*   **Implementation**: `get_trunc_normaliser_cpp`.
*   **Caching**: The value of $\log(Z)$ is cached (`inv_Z_cache`) keyed by parameters and bounds to avoid re-computation.
*   **Final Log-Likelihood**: $\log(L_{raw}) - \log(Z)$.

## 4. Brittleness and Recommendations

### 4.1. Brittleness in `compress_dadm`
As noted in section 2.3, `compress_dadm` likely fails to distinguish trials that differ *only* by `LT/UT/LC/UC` if those values are not parameters or covariates.
*   **Risk**: Incorrect likelihoods for mixed-boundary designs.
*   **Fix**: Ensure `LT`, `UT`, `LC`, `UC` are added to the uniqueness key in `compress_dadm` (e.g., via `paste(...)`).

### 4.2. Hardcoded Values
*   The C++ code assumes specific behavior for `Inf` and `-Inf`. Ensure data preparation consistently produces these values for censored trials.

### 4.3. Dependencies
The extraction will require the following components:
*   **Headers**: `src/utils.h` (contains `RacePdfFun`, `ContextForRaceModels`, and prototypes).
*   **Libraries**: `GSL` (GNU Scientific Library) for `gsl_integration_qags`, etc.
*   **Helpers**: `log_sum_exp` (likely in `utility_functions.h`).

## 5. Files for Extraction

For the Pull Request, you will need to extract logic from:

1.  **`src/particle_ll.cpp`**:
    *   `c_log_likelihood_race_cens_trunc`
    *   `integrate_for_kth_winner_cpp`
    *   `gsl_f_race_adapter`
    *   `get_trunc_normaliser_cpp`
    *   `calc_ll` (wrapper)

2.  **`src/utils.h`**:
    *   Struct definitions (`ContextForRaceModels`, `gsl_race_params`).
    *   Typedefs (`RacePdfFun`, `RaceCdfFun`).

3.  **`R/design.R`**:
    *   `compress_dadm` (needs modification for robustness).
    *   `design_model` (ensure `LT`... columns are preserved).

4.  **`R/make_data.R`**:
    *   `make_missing` (simulation logic).
