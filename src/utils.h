#ifndef PARTICLE_LL_UTILS_H
#define PARTICLE_LL_UTILS_H
#include <RcppArmadillo.h>
#include "utility_functions.h" // For any existing utilities we might use
// --- Internal Helper Function Declarations (not directly exported to R via Rcpp::export) ---
// Function pointer types for model's dfun and pfun
// Added LogicalVector is_ok for bound checking and void* context for extra args
typedef Rcpp::NumericVector (*RacePdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
										  Rcpp::LogicalVector winner,
										  Rcpp::LogicalVector is_ok,
                                          void* model_specific_context);
typedef Rcpp::NumericVector (*RaceCdfFun)(Rcpp::NumericVector rt,
                                          Rcpp::NumericMatrix pars,
										  Rcpp::LogicalVector winner,
										  Rcpp::LogicalVector is_ok,
                                          void* model_specific_context);


// Integrand function equivalent to R's f_race_integrand
// p_trial_this_winner_first: (n_accumulators x n_model_params)
// model_dfun, model_pfun: function pointers to the model's density and CDF
double f_race_integrand_cpp(
    double t,
    const Rcpp::NumericMatrix& p_trial_this_winner_first, // Parameters for the winner (1st row) and losers
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    void* model_specific_context // Added context
);

Rcpp::NumericVector f_race_integrand_batch_cpp(
    const Rcpp::NumericVector& rts_batch,
    const Rcpp::List& pars_all_trials_ordered, // List of NumericMatrix, each n_acc x n_params
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    int n_acc,
    void* model_specific_context // Added context
);

double c_log_likelihood_race_cens_trunc(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, covering all dadm rows for that particle
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions, structured for all accumulators
    RacePdfFun model_dfun,                  // Pointer to the model's PDF adapter function
    RaceCdfFun model_pfun,                  // Pointer to the model's CDF adapter function
	const int n_trials,
	const Rcpp::LogicalVector winner,
	const Rcpp::IntegerVector expand,  // Vector for expanding unique LLs to full trial count
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector ok_params,   // Parameter validity for each row in 'pars' matrix
    int n_acc,                              // Number of accumulators in the race (must be > 0 if data exists)
    void* model_context_for_funcs           // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
);

// Truncation correction factor helper
// Returns the factor itself (not log)
double get_trunc_corr_factor_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    double LT, double UT,
    int n_acc,
    double epsilon,
    void* model_specific_context
);

// Helper to order parameters (winner first)
Rcpp::NumericMatrix order_pars_for_winner_cpp(
    const Rcpp::NumericMatrix& p_all_acc, // n_acc x n_params
    int k_idx, // 1-based index of the winner
    int n_acc
);

// Structures
struct ContextForRaceModels {
    double min_lik_for_pdf;
    bool use_posdrift;
};

// Adapter prototypes
Rcpp::NumericVector lba_dfun_adapter(Rcpp::NumericVector rt,
                                     Rcpp::NumericMatrix pars,
                                     Rcpp::LogicalVector is_ok,
                                     Rcpp::LogicalVector winner,
                                     void* context);
Rcpp::NumericVector lba_pfun_adapter(Rcpp::NumericVector rt,
                                     Rcpp::NumericMatrix pars,
                                     Rcpp::LogicalVector is_ok,
                                     Rcpp::LogicalVector winner,
                                     void* context);
Rcpp::NumericVector rdm_dfun_adapter(Rcpp::NumericVector rt,
                                     Rcpp::NumericMatrix pars,
                                     Rcpp::LogicalVector is_ok,
                                     Rcpp::LogicalVector winner,
                                     void* context);
Rcpp::NumericVector rdm_pfun_adapter(Rcpp::NumericVector rt,
                                     Rcpp::NumericMatrix pars,
                                     Rcpp::LogicalVector is_ok,
                                     Rcpp::LogicalVector winner,
                                     void* context);
Rcpp::NumericVector lnr_dfun_adapter(Rcpp::NumericVector rt,
                                     Rcpp::NumericMatrix pars,
                                     Rcpp::LogicalVector is_ok,
                                     Rcpp::LogicalVector winner,
                                     void* context);
Rcpp::NumericVector lnr_pfun_adapter(Rcpp::NumericVector rt,
                                     Rcpp::NumericMatrix pars,
                                     Rcpp::LogicalVector is_ok,
                                     Rcpp::LogicalVector winner,
                                     void* context);
									 
									 // Static adapter for LBA dfun
Rcpp::NumericVector lba_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    // Pass use_posdrift from context to dlba_c
    return dlba_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for LBA pfun
Rcpp::NumericVector lba_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    // Pass use_posdrift from context to plba_c
    return plba_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for RDM dfun
Rcpp::NumericVector rdm_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    return drdm_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for RDM pfun
Rcpp::NumericVector rdm_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    return prdm_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}
// Static adapter for RDM dfun
Rcpp::NumericVector rdmswtn_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    return drdmswtn_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for RDM pfun
Rcpp::NumericVector rdmswtn_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context);
    return prdmswtn_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for LNR dfun
Rcpp::NumericVector lnr_dfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context); // Reusing same context struct
    return dlnr_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// Static adapter for LNR pfun
Rcpp::NumericVector lnr_pfun_adapter(Rcpp::NumericVector rt,
                                            Rcpp::NumericMatrix pars,
                                            Rcpp::LogicalVector is_ok,
                                            Rcpp::LogicalVector winner,
                                            void* context) {
    ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(context); // Reusing same context struct
    return plnr_c(rt, pars, winner, ctx->min_lik_for_pdf, is_ok);
}

// ---- END: Context Structs and Static Adapters ----

// ---- GSL-based helper structures and adapters ----
struct gsl_race_params {
    const Rcpp::NumericMatrix* p_trial;
	Rcpp::LogicalVector winner;
	Rcpp::LogicalVector isok;
    RacePdfFun model_dfun;
    RaceCdfFun model_pfun;
    int n_acc;
    void* model_specific_context;
};

double gsl_f_race_adapter(double t, void *p);

// Helper functions
Rcpp::NumericMatrix order_pars_for_winner_cpp(const Rcpp::NumericMatrix& p_all_acc,
                                             int k_idx,
                                             int n_acc);

Rcpp::NumericVector f_race_integrand_batch_cpp(const Rcpp::NumericVector& rts_batch,
                                              const Rcpp::List& pars_all_trials_ordered,
                                              Rcpp::NumericVector (*model_dfun)(Rcpp::NumericVector,
                                                                               Rcpp::NumericMatrix,
                                                                               Rcpp::LogicalVector,
                                                                               bool,
                                                                               void*),
                                              Rcpp::NumericVector (*model_pfun)(Rcpp::NumericVector,
                                                                               Rcpp::NumericMatrix,
                                                                               Rcpp::LogicalVector,
                                                                               bool,
                                                                               void*),
                                              int n_acc,
                                              void* model_specific_context);

double integrate_for_kth_winner_cpp(int k_winner_idx,
                                    const Rcpp::NumericMatrix& p_all_acc,
									Rcpp::LogicalVector isok,
									double low,
                                    double upp,
                                    Rcpp::NumericVector (*model_dfun)(Rcpp::NumericVector,
                                                                     Rcpp::NumericMatrix,
                                                                     Rcpp::LogicalVector,
                                                                     bool,
                                                                     void*),
                                    Rcpp::NumericVector (*model_pfun)(Rcpp::NumericVector,
                                                                     Rcpp::NumericMatrix,
                                                                     Rcpp::LogicalVector,
                                                                     bool,
                                                                     void*),
                                    int n_acc,
                                    double epsilon,
                                    void* model_specific_context);

double get_trunc_corr_factor_for_kth_winner_cpp(int k_winner_idx,
                                                const Rcpp::NumericMatrix& p_all_acc,
                                                Rcpp::NumericVector (*model_dfun)(Rcpp::NumericVector,
                                                                                 Rcpp::NumericMatrix,
                                                                                 Rcpp::LogicalVector,
                                                                                 bool,
                                                                                 void*),
                                                Rcpp::NumericVector (*model_pfun)(Rcpp::NumericVector,
                                                                                 Rcpp::NumericMatrix,
                                                                                 Rcpp::LogicalVector,
                                                                                 bool,
                                                                                 void*),

												LogicalVector winner_trial, 
												LogicalVector isok_trial, 
                                                double LT, double UT,
                                                int n_acc,
                                                double epsilon,
                                                void* model_specific_context);

// Exported wrappers
Rcpp::NumericVector calc_ll(Rcpp::NumericMatrix p_matrix, Rcpp::DataFrame data,
                            Rcpp::NumericVector constants, Rcpp::List designs,
                            Rcpp::String type_rcpp, Rcpp::List bounds,
                            Rcpp::List transforms, Rcpp::List pretransforms,
                            Rcpp::CharacterVector p_types, double min_ll,
                            Rcpp::List trend);
							
double get_trunc_normaliser_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
	LogicalVector isok_trial, 
    double LT, double UT,
    int n_acc,
    double epsilon = 1e-8,
    void* model_specific_context = nullptr);


#endif