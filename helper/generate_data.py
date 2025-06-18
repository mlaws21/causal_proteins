import random

# vars effecting disease (response):
# - age (old/young) old much more likely to have cancer
# - race (white/non-white) white slightly more likely
# - lifestyle (healthy/unhealthy) unhealty more likely

# vars effecting mutations
# - age (old/young) old more likely to have mutation
# - race  (white/non-white) -- what if we split the stuff into two sub pops and weighted sample
# - lifestyle (healthy/unhealthy) unhealty more likely
def generate_synth(benign_seqs, pathogenic_seqs, num_datapoints=1000):
    random.seed(42)
    
    benign_seqs = list(benign_seqs.itertuples(index=False, name=None))
    pathogenic_seqs = list(pathogenic_seqs.itertuples(index=False, name=None))

    rows = []
    
    for _ in range(num_datapoints):
        age = random.choice([True, False])
        race = random.choice([True, False])
        lifestyle = random.choice([True, False])
    
            
        pathogenic_chance = 0.35
        
        if lifestyle:
            pathogenic_chance += 0.1
            
        if age:
            pathogenic_chance += 0.1

        if race:
            pathogenic_chance += 0.1
            
        is_pathogenic = random.uniform(0, 1) < pathogenic_chance
        

        if is_pathogenic:
            to_sample = pathogenic_seqs
        else:
            to_sample = benign_seqs
    
        sample = random.choice(to_sample)
            
        idx = sample[0]
        seq = sample[1]
        clinical_sig = sample[2]
        
        row = {
                "ID": idx,
                "Age": age,
                "Race": race,
                "Lifestyle": lifestyle,
                "Sequence_Score": None,
                "Align_Score": None,
                "Disease": None,
                "Sequence": seq,
                "Ground": is_pathogenic}
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    return df
  
def split_and_unique(df):
    
    benign = df[df["clinical_sig"] == False].drop_duplicates(subset="sequence")
    disease = df[df["clinical_sig"] == True].drop_duplicates(subset="sequence")
    
    return benign, disease


def main():
    pass
    
    
    
if __name__ == "__main__":
    main()