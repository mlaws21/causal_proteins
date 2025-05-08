import re
import random
random.seed(42)
import numpy as np


# 1) Parameters
file_path = 'prion429none_effect.log'    # your input file
threshold = 0.0                # decide: score >= threshold ⇒ predict Pathogenic

# 2) Parse file
labels = []    # 1=Pathogenic, 0=Benign
lower_bound = []
upper_bound = []
point_est = []

pattern = re.compile(
    r'treatment_[^\[]+\[(?P<label>Benign|Pathogenic)\]:\s*'      # label
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
print(f"Total mutations: {N}")
print(f"Pathogenic: {num_pos} ({num_pos/N*100:.1f}%)")
print(f"Benign:     {num_neg} ({num_neg/N*100:.1f}%)\n")

# 4) Confusion matrix at `threshold`
preds = [1 if s>=threshold else 0 for s in lower_bound]

TP = sum(1 for p,l in zip(preds,labels) if p==1 and l==1)
FP = sum(1 for p,l in zip(preds,labels) if p==1 and l==0)
FN = sum(1 for p,l in zip(preds,labels) if p==0 and l==1)
TN = sum(1 for p,l in zip(preds,labels) if p==0 and l==0)

print("Confusion matrix @ threshold =", threshold)
print(f"  TP: {TP}")
print(f"  FP: {FP}")
print(f"  FN: {FN}")
print(f"  TN: {TN}\n")

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


avg_path_rank = total_rank / num_pos


random_avgs = []
for _ in range(10000):
    total_rank = 0
    
    random.shuffle(labels)
    for i,ele in enumerate(labels, 1):
        if ele == 1:
            total_rank += i
    random_avgs.append(total_rank / num_pos)
    
random_avgs = np.array(random_avgs)

ci_lower = np.percentile(random_avgs, 2.5)
ci_upper = np.percentile(random_avgs, 97.5)
    
print(f"Average rank of Pathogenic mutations: {avg_path_rank:.2f} (out of {N})")   
print(f"Expected rank of Pathogenic mutations: {np.mean(random_avgs):.2f} ({ci_lower:.2f}, {ci_upper:.2f}) (out of {N})")   


        
# sorted_idxs = sorted(range(N), key=lambda i: point_est[i], reverse=True)
# print(sorted_idxs)
# print(point_est)
# ranks = [0]*N
# for rank, idx in enumerate(sorted_idxs, start=1):
#     ranks[idx] = rank

# print(ranks)
# path_ranks = [ranks[i] for i in range(N) if labels[i]==1]
# print(sorted(path_ranks))
# avg_path_rank = sum(path_ranks) / len(path_ranks)
# print(f"Average rank of Pathogenic mutations: {avg_path_rank:.1f} (out of {N})")
