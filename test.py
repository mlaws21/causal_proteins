import re
import numpy as np

def parse_and_evaluate(text, num_simulations=100000):
    # Parse all treatments
    pattern = r"treatment_(\S+) \[(\w+)\]: Causal Effect: ([\-\d\.]+)"
    matches = re.findall(pattern, text)

    # Build structured list
    data = []
    for name, label, effect in matches:
        data.append((name, label, float(effect)))

    # Sort by causal effect (decending)
    data.sort(key=lambda x: x[2], reverse=True)

    # Find positions of pathogenic entries
    pathogenic_positions = [i for i, (_, label, _) in enumerate(data) if label == 'Pathogenic']

    if not pathogenic_positions:
        raise ValueError("No pathogenic entries found.")

    avg_pathogenic_position = np.mean(pathogenic_positions)

    n = len(data)
    k = len(pathogenic_positions)
    expected_random_position = (n - 1) / 2  # Expected mean

    # Simulation to estimate confidence interval
    random_avgs = []
    for _ in range(num_simulations):
        sampled_positions = np.random.choice(n, size=k, replace=False)
        random_avgs.append(np.mean(sampled_positions))
    random_avgs = np.array(random_avgs)

    ci_lower = np.percentile(random_avgs, 2.5)
    ci_upper = np.percentile(random_avgs, 97.5)

    return {
        "num_total_entries": n,
        "num_pathogenic": k,
        "avg_pathogenic_position": avg_pathogenic_position,
        "expected_random_position": expected_random_position,
        "random_ci_95": (ci_lower, ci_upper)
    }


# Example usage
with open("prion506_effect.log", "r") as f:
    text = f.read()

results = parse_and_evaluate(text)
print(results)
