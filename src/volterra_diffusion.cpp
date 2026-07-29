// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Volterra reducible-diffusion solvers (leaky accumulation).
//
// This translation unit exists purely to give the Volterra type-2 solver
// headers a home in the build.  The headers define non-inline free functions,
// so they must be included in exactly ONE .cpp -- this one.  Nothing else in
// the package includes them, and none of the exported entry points below are
// referenced by the model dispatch / design machinery yet: they are reachable
// from R (via compileAttributes) for development and validation only.
//
// See utils_reducible_diffusion.h for the shared Abel/Stieltjes machinery,
// model_OU_Volterra.h for the Ornstein-Uhlenbeck first-passage solver, and
// model_BM_Volterra.h for the driftless-Brownian-motion reduction.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "model_OU_Volterra.h"
#include "model_BM_Volterra.h"
