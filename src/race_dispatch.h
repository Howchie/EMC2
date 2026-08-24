#ifndef EMC2_RACE_DISPATCH_H
#define EMC2_RACE_DISPATCH_H

#include <Rcpp.h>
#include <string>

#include "race_contract.h"
#include "col_registry.h"

struct RaceModelAdapter {
  RacePdf1Fun pdf1_ptr = nullptr;
  RaceCdf1Fun cdf1_ptr = nullptr;
  RaceRawFun model_dfun_raw = nullptr;   // fast-path: write log-density to pre-allocated buffer
  RaceRawFun model_pfun_raw = nullptr;   // fast-path: write log-survivor to pre-allocated buffer
  RaceLogSAtTFun logS_at_t_ptr = nullptr; // batch: log-survivor at scalar t for truncation norms
  emc2col::ColSpec col_spec = {nullptr, 0, ""}; // kernel column contract (col_registry.h)
  ContextForRaceModels ctx;
};

RaceModelAdapter resolve_race_model_adapter(const std::string& type_std,
                                            const std::string& caller);

void configure_corr_drift_context(RaceModelAdapter& adapter,
                                  const Rcpp::CharacterVector& keep_names,
                                  const std::string& caller);

void configure_rdmswtn_corr_context(
    RaceModelAdapter& adapter, const Rcpp::CharacterVector& keep_names,
    const std::string& caller);

bool is_stop_signal_type(const std::string& type_std);

#endif
