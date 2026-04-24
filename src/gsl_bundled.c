/*
 * gsl_bundled.c — single compilation unit for the minimal GSL subset used by EMC2.
 *
 * Only the integration (QAGS/QAGIU/QAG/QNG/workspace) and error-handler routines are
 * included; all statistical GSL functions (erfc, gaussian_tail, …) are already
 * reimplemented as emc2_* variants in tools.cpp and are NOT compiled here.
 *
 * Source: GNU Scientific Library 2.7.1 (LGPL v2.1+), bundled to avoid a
 * system-level GSL installation requirement.
 *
 * Copyright notices for each file are preserved in the corresponding source
 * files under src/gsl_bundled/.
 */

/* Suppress the `#include <config.h>` that each GSL .c file starts with.
 * Our config.h stub lives in gsl_bundled/ and is found via -Igsl_bundled. */

/* --- Error handler -------------------------------------------------------- */
#include "gsl_bundled/err/strerror.c"
#include "gsl_bundled/err/stream.c"
#include "gsl_bundled/err/error.c"

/* --- Integration workspace ------------------------------------------------ */
#include "gsl_bundled/integration/workspace.c"

/* --- Gauss-Kronrod generic kernel ----------------------------------------- */
/* qk15.c and qk21.c are in separate TUs (gsl_bundled_qk15/21.c) because they
 * both define static arrays named xgk/wg/wgk. */
#include "gsl_bundled/integration/qk.c"

/* --- QAGS + QAGIU (qags.c pulls in its own helpers via relative #includes) - */
#include "gsl_bundled/integration/qags.c"
/* QAG and QNG are compiled in separate TUs (gsl_bundled_qag/qng.c) to avoid
 * helper symbol collisions from nested #includes inside the original sources. */
