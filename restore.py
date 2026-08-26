import re

with open('recovered.txt', 'r') as f:
    lines = f.readlines()

code_to_insert = []
for line in lines:
    m = re.match(r'^\d+:\s(.*)', line)
    if m:
        code_to_insert.append(m.group(1) + '\n')

# Only take the lines up to 658
# In recovered.txt, the last few lines are 657, 658, 659, 660.
# We want from 430 up to 658 inclusive.
final_code = ""
for line in code_to_insert:
    if line.strip() == "// A natural probability difference is rejected once it retains less than":
        break
    final_code += line

with open('src/wald_functions.h', 'r') as f:
    orig = f.read()

marker = "// Shared acceptance constants for the guarded natural race kernels."
if marker in orig:
    new_text = orig.replace(marker, final_code + "\n" + marker)
    with open('src/wald_functions.h', 'w') as f:
        f.write(new_text)
    print("Successfully restored!")
else:
    print("Marker not found!")
