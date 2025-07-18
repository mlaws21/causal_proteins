
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import re
from sklearn.metrics import precision_recall_curve, average_precision_score
import random
random.seed(42)
import numpy as np

def ordering(project_name, log_fn):
# Load input
    with open(f"outputs/{project_name}/effect.log", "r") as f:
        lines = f.readlines()[1:-2]  # skip timestamp

    # Parse lines
    pattern = re.compile(r"(\w+)\s+\[(\w+)\]: Causal Effect:\s+([-+]?[0-9]*\.?[0-9]+)")
    data = []

    for line in lines:
        match = pattern.search(line)
        if match:
            treatment, label, effect = match.groups()
            data.append((line.strip(), float(effect), label.lower()))

    # Sort by effect
    data.sort(key=lambda x: x[1], reverse=True)

    # Parameters
    font = ImageFont.load_default()  # or use truetype if you want custom font
    line_height = 20
    padding = 10
    img_width = 400
    img_height = padding * 2 + line_height * len(data)

    # Create image
    img = Image.new("RGB", (img_width, img_height), color="white")
    draw = ImageDraw.Draw(img)

    # Draw each line
    for i, (line, _, label) in enumerate(data):
        y = padding + i * line_height
        color = "red" if label == "pathogenic" else "green"
        draw.text((padding, y), line, fill=color, font=font)

    # Save image
    out_name = f"outputs/{project_name}/results/ordering.png"
    img.save(out_name)
    log_fn(f"Saved image as {out_name}")

def generate_boxplot(project_name):
    file_path = f"outputs/{project_name}/effect.log"
    # 1) Read and parse the file, collecting point estimates by label
    # file_path = 'prion_effect.log'  # <-- update to your actual path/filename
    benign_scores = []
    pathogenic_scores = []

    pattern = re.compile(
        r'[^\[]+\[(Benign|Pathogenic)\]: '
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
    plt.savefig(f"outputs/{project_name}/results/bp.png", format="png", dpi=300)
    plt.savefig(f"outputs/{project_name}/results/bp.eps", format="eps", dpi=300)

def generate_pr_curve(project_name):
    # 1) Read and parse the file
    # file_path = 'prionsequence_effect.log'  # <-- replace with your filename
    
    file_path = f"outputs/{project_name}/effect.log"
    
    labels = []
    scores = []

    # pattern captures “Benign” vs “Pathogenic” and the causal‐effect score
    pattern = re.compile(
        r'[^\[]+\[(?P<label>Benign|Pathogenic)\]: '
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

    plt.savefig(f"outputs/{project_name}/results/prc.png", format="png", dpi=300)
    plt.savefig(f"outputs/{project_name}/results/prc.eps", format="eps", dpi=300)


def generate_summary(project_name, threshold=0.0):

    file_path = f"outputs/{project_name}/effect.log"

    labels = []    # 1=Pathogenic, 0=Benign
    lower_bound = []
    upper_bound = []
    point_est = []

    pattern = re.compile(
        r'[^\[]+\[(?P<label>Benign|Pathogenic)\]:\s*'      # label
        r'Causal Effect:\s*'
        r'(?P<point>-?\d+\.\d+)\s*'                                  # point estimate
        r'\(\s*'
        r'(?P<low>-?\d+\.\d+)\s*,\s*'                                # lower CI
        r'(?P<high>-?\d+\.\d+)\s*'                                   # upper CI
        r'\)'
    )

    with open(file_path) as f:
        for line in f:
            m = pattern.search(line)
            if not m:
                continue
            labels.append(1 if m.group('label')=='Pathogenic' else 0)
            lower_bound.append(float(m.group('low')))
            point_est.append(float(m.group('point')))
            upper_bound.append(float(m.group('high')))

    # 3) Dataset summary
    N = len(labels)
    num_pos = sum(labels)
    num_neg = N - num_pos


    # 4) Confusion matrix at `threshold`
    preds = [1 if s>=threshold else 0 for s in lower_bound]

    TP = sum(1 for p,l in zip(preds,labels) if p==1 and l==1)
    FP = sum(1 for p,l in zip(preds,labels) if p==1 and l==0)
    FN = sum(1 for p,l in zip(preds,labels) if p==0 and l==1)
    TN = sum(1 for p,l in zip(preds,labels) if p==0 and l==0)



    # 5) Average rank of Pathogenic items
    #    rank 1 = highest score, rank N = lowest score


    combo = list(zip(point_est, labels))

    scombo = sorted(combo, key=lambda x: x[0], reverse=True)
    num_pos = sum(labels)

    # print(scombo)
    total_rank = 0
    for i,ele in enumerate(scombo, 1):
        if ele[1] == 1:
            total_rank += i

    if num_pos == 0:
        avg_path_rank = 0.0
    else:
        avg_path_rank = total_rank / num_pos


    random_avgs = []
    for _ in range(10000):
        total_rank = 0
        
        random.shuffle(labels)
        for i,ele in enumerate(labels, 1):
            if ele == 1:
                total_rank += i
        if num_pos == 0:
            random_avgs.append(0.0)
        else:
            random_avgs.append(total_rank / num_pos)
        
    random_avgs = np.array(random_avgs)

    ci_lower = np.percentile(random_avgs, 2.5)
    ci_upper = np.percentile(random_avgs, 97.5)
        
        
    with open(f"outputs/{project_name}/results/summary.txt", "w") as f:
        print(f"Total mutations: {N}", file=f)
        print(f"Pathogenic: {num_pos} ({num_pos/N*100:.1f}%)", file=f)
        print(f"Benign:     {num_neg} ({num_neg/N*100:.1f}%)\n", file=f)
        print("Confusion matrix @ threshold =", threshold, file=f)
        print(f"  TP: {TP}", file=f)
        print(f"  FP: {FP}", file=f)
        print(f"  FN: {FN}", file=f)
        print(f"  TN: {TN}\n", file=f)
        print(f"Average rank of Pathogenic mutations: {avg_path_rank:.2f} (out of {N})", file=f)
        print(f"Expected rank of Pathogenic mutations: {np.mean(random_avgs):.2f} ({ci_lower:.2f}, {ci_upper:.2f}) (out of {N})", file=f)
