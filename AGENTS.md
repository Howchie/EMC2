# Repository requirements

- Every model constructor in this package (including race, diffusion,
  discrete-choice, signal-detection, and MRI models) MUST use a compiled C++
  simulator by default. A pure-R generator may exist only as an explicit
  opt-out/reference path (`options(emc2.cpp_rfun = FALSE)`).
- Adding or changing a model requires an exported C++ simulator entry point,
  its Rcpp binding, and a regression test that exercises the default path.
- Keep comments short and local; do not duplicate implementation details in
  prose when the code and tests already make the behavior clear.
