with open('src/wald_functions.h', 'r') as f:
    text = f.read()

# First, remove my forward declaration that I inserted at line 430
lines = text.split('\n')
new_lines = []
for line in lines:
    if line.strip() == "inline double dlnorm_std(double x, double meanlog, double sdlog, bool log_p);":
        continue
    new_lines.append(line)

text = '\n'.join(new_lines)

# Re-add the default argument in the definition
text = text.replace("inline double dlnorm_std(double x, double meanlog, double sdlog, bool log_p) {",
                    "inline double dlnorm_std(double x, double meanlog, double sdlog, bool log_p = false) {")

# Add forward declaration correctly:
marker = "inline double lnorm_log_surv_std(double x, double meanlog, double sdlog);"
fwd = "inline double dlnorm_std(double x, double meanlog, double sdlog, bool log_p = false);"
text = text.replace(marker, fwd + "\n" + marker)

with open('src/wald_functions.h', 'w') as f:
    f.write(text)

