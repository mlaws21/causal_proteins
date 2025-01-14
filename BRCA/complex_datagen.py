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
    

    path_data = mutation_data[mutation_data["Clinical_significance_ENIGMA"] == "Pathogenic"]
    benign_data = mutation_data[mutation_data["Clinical_significance_ENIGMA"] == "Benign"]
    
    path_55p = int(len(path_data) * 0.55)
    benign_55p = int(len(benign_data) * 0.55)
    
    
    # pop 1 is slightly higher chance of disease
    pop1 = pd.concat([path_data[:path_55p], benign_data[benign_55p:]], axis=0, ignore_index=True)
    pop2 = pd.concat([path_data[path_55p:], benign_data[:benign_55p]], axis=0, ignore_index=True)
    

    
    
    # path_data.to_csv('temp.csv', index=False)

    # pop1.to_csv('pop1.csv', index=False)
    # pop2.to_csv('pop2.csv', index=False)

    
    
    
    
    
    rows = []
    for _ in range(num_datapoints):
        age = random.choice([True, False])
        race = random.choice([True, False])
        lifestyle = random.choice([True, False])
        alignment = random.uniform(0, 1)
        cancer = random.uniform(0, 1) < (0.05 + (0.3 if age else 0.1) + (0.1 if race else 0.05) + (0.30 if lifestyle else 0.05))
    
        if race: #white
            # upweight pop1
            population = random.choices([pop1, pop2], [0.65, 0.35])
        
        else:
            population = random.choices([pop1, pop2], [0.35, 0.65])
            
        pathogenic_chance = 0.5
        
        if lifestyle:
            pathogenic_chance += 0.1
        else:
            pathogenic_chance -= 0.1
            
        if age:
            pathogenic_chance += 0.1
        else:
            pathogenic_chance -= 0.1
            
        is_pathogenic = random.uniform(0, 1) < pathogenic_chance
        
        # print(is_pathogenic)
        
        population = population[0]

        if is_pathogenic:
            to_sample = population[population["Clinical_significance_ENIGMA"] == "Pathogenic"]
        else:
            to_sample = population[population["Clinical_significance_ENIGMA"] == "Benign"]
        
        sample = to_sample.sample(n=1, random_state=2)
        seq = sample["seq"].iloc[0]
        
        
            
            

        row = {"Old": age,
               "White": race,
               "Unhealthy": lifestyle,
               "Align Score": alignment,
               "Cancer": cancer,
               "Sequence": seq}
        
        rows.append(row)
        
        
    
    
    
    df = pd.DataFrame(rows)
    
    df.to_csv('data.csv', index=False)
    
generate_data("brca_mutations.csv")