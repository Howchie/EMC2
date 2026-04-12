/* Separate compilation unit for qk15 — cannot share a TU with qk21
 * because both define static arrays named xgk/wg/wgk. */
#include "gsl_bundled/integration/qk15.c"
