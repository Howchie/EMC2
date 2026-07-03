lines <- readLines("src/fncs_seven.cpp")
idx1 <- grep("gsl_integration_qag\\(&F, xmin\\[0\\], xmax\\[0\\]", lines)
for (i in idx1) {
  lines[i] <- "    size_t neval;\n    int status = gsl_integration_qng(&F, xmin[0], xmax[0], abstol, reltol, &val, &abserr, &neval);"
}
idx2 <- grep("GSL_INTEG_GAUSS15, ws, &val, &abserr\\);", lines)
for (i in idx2) {
  lines[i] <- ""
}
idx3 <- grep("gsl_integration_workspace\\* ws = ensure_gsl_workspace", lines)
for (i in idx3) {
  lines[i] <- ""
}
writeLines(lines, "src/fncs_seven.cpp")
