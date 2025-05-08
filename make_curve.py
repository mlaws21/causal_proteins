import re
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score

# 1) Read and parse the file
file_path = 'prion506_effect.log'  # <-- replace with your filename
labels = []
scores = []

# pattern captures “Benign” vs “Pathogenic” and the causal‐effect score
pattern = re.compile(
    r'treatment_[^\[]+\[(?P<label>Benign|Pathogenic)\]: '
    r'Causal Effect: -?\d+\.\d+ \('
    r'(?P<low>-?\d+\.\d+),'
)

with open(file_path, 'r') as f:
    for line in f:
        m = pattern.search(line)
        if not m:
            continue
        label = m.group('label')
        score = float(m.group('low'))
        # treat Pathogenic as the “positive” class
        labels.append(1 if label == 'Pathogenic' else 0)
        scores.append(score)

# 2) Compute precision and recall
precision, recall, _ = precision_recall_curve(labels, scores)
avg_prec = average_precision_score(labels, scores)

# 3) Plot
plt.figure(figsize=(8, 6))
plt.plot(recall, precision, lw=2,
         label=f'PR curve (AP = {avg_prec:.2f})')

plt.xlabel('Recall', fontsize=14)
plt.ylabel('Precision', fontsize=14)
plt.title('Precision–Recall Curve', fontsize=16)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='lower left')
plt.tight_layout()

plt.savefig("prc.png", format="png", dpi=300)
plt.savefig("prc.eps", format="eps", dpi=300)


