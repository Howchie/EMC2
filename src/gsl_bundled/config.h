/* Minimal config.h for bundled GSL — replaces the autoconf-generated one.
 * Only defines what the integration/ and err/ source files actually test.
 */
#ifndef GSL_BUNDLED_CONFIG_H
#define GSL_BUNDLED_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_C99_INLINE 1

/* From GSL configure.ac — used in workspace.c free() guard */
#define RETURN_IF_NULL(x) if (!(x)) { return ; }

/* Silence the PACKAGE/VERSION macros sometimes expected by GSL internals */
#ifndef PACKAGE
#  define PACKAGE "gsl"
#endif
#ifndef VERSION
#  define VERSION "2.7.1"
#endif

#endif /* GSL_BUNDLED_CONFIG_H */
