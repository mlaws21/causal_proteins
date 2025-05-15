import re
import matplotlib.pyplot as plt
import sys
file_path = sys.argv[1]
# 1) Read and parse the file, collecting point estimates by label
# file_path = 'prion_effect.log'  # <-- update to your actual path/filename
benign_scores = []
pathogenic_scores = []

pattern = re.compile(
    r'treatment_[^\[]+\[(Benign|Pathogenic)\]: '
    r'Causal Effect: (?P<score>-?\d+\.\d+)'
)

with open(file_path, 'r') as f:
    for line in f:
        m = pattern.search(line)
        if not m:
            continue
        label = m.group(1)
        score = float(m.group('score'))
        if label == 'Benign':
            benign_scores.append(score)
        else:
            pathogenic_scores.append(score)

# 2) Create side-by-side vertical boxplots
plt.figure(figsize=(8, 6))
plt.boxplot([benign_scores, pathogenic_scores], labels=['Benign', 'Pathogenic'])
plt.ylabel('Causal Effect', fontsize=12)
plt.title('Distribution of Causal Effect Estimates', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# 3) Save the figure (instead of plt.show())
plt.savefig("bp.png", format="png", dpi=300)
plt.savefig("bp.eps", format="eps", dpi=300)
