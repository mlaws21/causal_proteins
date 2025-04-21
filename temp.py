import re
from pathlib import Path

# Input file path
input_file = "prion411_effect.log"
output_file = "sorted_causal_effects.html"

# Read raw data
with open(input_file, "r") as f:
    lines = f.readlines()[1:]  # skip timestamp

# Parse lines
data = []
pattern = re.compile(r"treatment_(\w+)\s+\[(\w+)\]: Causal Effect:\s+([-+]?[0-9]*\.?[0-9]+)")

for line in lines:
    match = pattern.search(line)
    if match:
        treatment, label, effect = match.groups()
        data.append((treatment, label, float(effect), line.strip()))

# Sort by causal effect
data.sort(key=lambda x: x[2], reverse=True)

# Generate HTML output
html_lines = [
    "<html><head><title>Sorted Causal Effects</title></head><body style='font-family:monospace;'>",
    "<h2>Sorted by Causal Effect (Descending)</h2><ul>"
]

for treatment, label, effect, full_line in data:
    color = "red" if label.lower() == "pathogenic" else "green"
    html_lines.append(f"<li style='color:{color}'>{full_line}</li>")

html_lines.append("</ul></body></html>")

# Save HTML
Path(output_file).write_text("\n".join(html_lines))
print(f"Saved sorted output to {output_file}")
