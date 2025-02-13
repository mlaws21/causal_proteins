# vars effecting brest cancer (response):
# - age (old/young) old much more likely to have cancer
# - race (white/non-white) white slightly more likely
# - lifestyle (healthy/unhealthy) unhealty more likely

# vars effecting mutations
# - age (old/young) old more likely to have mutation
# - race  (white/non-white) -- what if we split the stuff into two sub pops and weighted sample
# - lifestyle (healthy/unhealthy) unhealty more likely

import pandas as pd
import random



def generate_data(mutation_path, num_datapoints=1000):
    
    mutation_data = pd.read_csv(mutation_path)
    
    
    rows = []
    for _ in range(num_datapoints):
        age = random.choice([True, False])
        race = random.choice([True, False])
        lifestyle = random.choice([True, False])
        alignment = random.uniform(0, 1)
        cancer = random.uniform(0, 1) < (0.05 + (0.3 if age else 0.1) + (0.1 if race else 0.05) + (0.30 if lifestyle else 0.05))
    
            
        pathogenic_chance = 0.35
        
        if lifestyle:
            pathogenic_chance += 0.1
            
        if age:
            pathogenic_chance += 0.1

        if race:
            pathogenic_chance += 0.1
            
        is_pathogenic = random.uniform(0, 1) < pathogenic_chance
        
        # print(is_pathogenic)

        if is_pathogenic:
            to_sample = mutation_data[mutation_data["Clinical_significance_ENIGMA"] == "Pathogenic"]
        else:
            to_sample = mutation_data[mutation_data["Clinical_significance_ENIGMA"] == "Benign"]
        
        print(len(to_sample))
        sample = to_sample.sample(n=1)
        seq = sample["seq"].iloc[0]
        
        idx = sample["Genomic_Coordinate_hg38"].iloc[0]
            

        row = {
                "ID": idx,
                "Old": age,
                "White": race,
                "Unhealthy": lifestyle,
                "Align Score": alignment,
                "Cancer": cancer,
                "Sequence": seq,
                "is_pathogenic": is_pathogenic}
        
        rows.append(row)
        
        
    
    
    
    df = pd.DataFrame(rows)
    
    df.to_csv('data.csv', index=False)
    
generate_data("brca_mutations.csv")