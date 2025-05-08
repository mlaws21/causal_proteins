# def calc_summary(data, cutoff=0):
#     """
#     Calculates classification summary statistics based on a manually chosen cutoff.

#     Parameters:
#     - data: List of tuples in the form (mutation, point_estimate, lower_CI, upper_CI, label)
#             where label is either 'Benign' or 'Pathogenic'
#     - cutoff: The decision boundary. If lower_CI >= cutoff, we predict "Pathogenic".

#     Returns:
#     - A dictionary with TP, FP, TN, FN, FPR, TPR, accuracy, precision
#     """

#     benign_pos = 0  # False Positives
#     benign_neg = 0  # True Negatives
#     path_pos = 0    # True Positives
#     path_neg = 0    # False Negatives

#     for i in data:
#         mutation, point, bottom, top, lab = i

#         if bottom < cutoff:
#             # Predict "Benign"
#             if lab == "Benign":
#                 benign_neg += 1  # True Negative
#             elif lab == "Pathogenic":
#                 path_neg += 1    # False Negative
#         else:
#             # Predict "Pathogenic"
#             if lab == "Benign":
#                 benign_pos += 1  # False Positive
#             elif lab == "Pathogenic":
#                 path_pos += 1    # True Positive

#     TP = path_pos
#     FP = benign_pos
#     TN = benign_neg
#     FN = path_neg

#     # Derived metrics
#     FPR = FP / (FP + TN) if (FP + TN) > 0 else 0.0
#     TPR = TP / (TP + FN) if (TP + FN) > 0 else 0.0
#     Accuracy = (TP + TN) / (TP + FP + TN + FN) if (TP + FP + TN + FN) > 0 else 0.0
#     Precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0

#     # Print results
#     print(f"Cutoff: {cutoff}")
#     print(f"True Positives (Pathogenic predicted Pathogenic): {TP}")
#     print(f"False Positives (Benign predicted Pathogenic): {FP}")
#     print(f"True Negatives (Benign predicted Benign): {TN}")
#     print(f"False Negatives (Pathogenic predicted Benign): {FN}")
#     print(f"False Positive Rate (FPR): {FPR:.3f}")
#     print(f"True Positive Rate (Recall/TPR): {TPR:.3f}")
#     print(f"Accuracy: {Accuracy:.3f}")
#     print(f"Precision: {Precision:.3f}")

#     return {
#         "TP": TP, "FP": FP, "TN": TN, "FN": FN,
#         "FPR": FPR, "TPR": TPR,
#         "accuracy": Accuracy,
#         "precision": Precision
#     }


# import numpy as np
# import matplotlib.pyplot as plt

# def sweep_cutoffs_for_precision(data, thresholds, min_recall=0.0):
#     results = []

#     for cutoff in thresholds:
#         metrics = calc_summary(data, cutoff)
#         metrics["cutoff"] = cutoff
#         results.append(metrics)

#     return results

# def plot_precision_with_min_recall(results, min_recall=0.0):
#     filtered = [r for r in results if r["TPR"] >= min_recall]

#     if not filtered:
#         print(f"No cutoffs meet the minimum recall of {min_recall}")
#         return None

#     cutoffs = [r["cutoff"] for r in filtered]
#     precisions = [r["precision"] for r in filtered]
#     recalls = [r["TPR"] for r in filtered]

#     best_idx = np.argmax(precisions)
#     best_cutoff = cutoffs[best_idx]

#     plt.figure(figsize=(12, 6))
#     plt.plot(cutoffs, precisions, marker='o', label="Precision")
#     plt.axvline(best_cutoff, color='red', linestyle='--', label=f"Best cutoff = {best_cutoff:.3f}")
    
#     # Annotate each point with its precision
#     for x, y in zip(cutoffs, precisions):
#         plt.text(x, y + 0.01, f"{y:.2f}", fontsize=8, ha='center')

#     plt.xlabel("Cutoff")
#     plt.ylabel("Precision")
#     plt.title(f"Precision vs Cutoff (min recall = {min_recall})")
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig("curr_precison_curve.png", format="png", dpi=300)

#     print(f"Best cutoff for precision (recall ≥ {min_recall}): {best_cutoff:.4f}")
#     print(f"Precision: {precisions[best_idx]:.3f} | Recall: {recalls[best_idx]:.3f}")
    
#     return best_cutoff


# def plot_cutoff(data):

#     bottoms = [x[2] for x in data]
#     thresholds = np.linspace(0, max(bottoms), 100)
#     results = sweep_cutoffs_for_precision(data, thresholds, min_recall=0.4)
#     best_cutoff = plot_precision_with_min_recall(results, min_recall=0.4)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def calc_summary(data, cutoff=0):
    benign_pos = 0
    benign_neg = 0
    path_pos = 0
    path_neg = 0

    for i in data:
        mutation, point, bottom, top, lab = i

        if bottom < cutoff:
            if lab == "Benign":
                benign_neg += 1  # TN
            elif lab == "Pathogenic":
                path_neg += 1    # FN
        else:
            if lab == "Benign":
                benign_pos += 1  # FP
            elif lab == "Pathogenic":
                path_pos += 1    # TP

    TP = path_pos
    FP = benign_pos
    TN = benign_neg
    FN = path_neg

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0

    return precision, recall

def generate_precision_recall_curve(data, save_path="precision_recall_curve.png"):
    precisions = []
    recalls = []
    
    thresholds = sorted([x[2] for x in data])
    # thresholds = np.linspace(0, max(bottoms), 100)
    # print(thresholds)
    for t in thresholds:
        p, r = calc_summary(data, t)
        precisions.append(p)
        recalls.append(r)
        
    plt.figure(figsize=(8, 6))
    
    for r, p, t in zip(recalls, precisions, thresholds):
        plt.text(r, p + 0.015, f'{t:.2f}', fontsize=8, ha='center')

    # Plot
    plt.plot(recalls, precisions, marker='o')
    plt.xlabel("Recall (TPR)")
    plt.ylabel("Precision")
    plt.title("Precision–Recall Curve")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(save_path)
    print(f"Saved precision-recall curve to: {save_path}")
