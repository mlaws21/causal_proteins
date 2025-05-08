import pandas as pd


good_swaps = [("e", "d"), ("q", "e"), ("y", "h"), ("k", "r"), ("i", "v"), ("l", "i"), ("l", "m"), ("y", "f"), ("y", "w")]

best_swaps = [("i", "v"), ("y", "f")]

path_muts = ['e211d', 'q217r', 'p105t', 'm232r', 'r208h', 'm232t', 'p39l', 'v203i', 'a133v', 'v180i', 'e196k', 't188k', 'a117v', 'v176g', 'v189i', 'e200g', 'h187r', 'm129v', 'q212p', 'g131e', 'g131v', 'p102l', 'e200k', 'e200k_m129v', 'a116v', 'p105l', 's17g', 't183a', 'v210i', 'e211q', 'r208c', 'f198s', 'g114v', 'p105s', 'd202n', 'e196a', 'd178n', 'i215v']

ref = "MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"

ref = ref.lower()

# known
benign_muts = [
    "a2v", "g54s", "g127v", "r136s", "g142s", 
    "n171s", "e219k", "l234f", "d18e", "g142c", "h177q"
]

for swap1, swap2 in best_swaps:
    for i, char in enumerate(ref):
        if i == 0:
            continue
        if char == swap1:
            benign_muts.append(f"{swap1}{i+1}{swap2}")
        if char == swap2:
            benign_muts.append(f"{swap2}{i+1}{swap1}")
            
print(len(set(benign_muts)))

print(len(set(path_muts)))


def mutate_protein(ref, mutations, log_fn=print):
    # mutations should be a list of XNX format mutations
    # you can also give mutation groups which are XNX mutations 
    # separate by "_": XNX_XNX...
    # stop codon is denoted "-"
    
    ref = ref.lower()
    
    for mut in mutations:
        
        outchar = mut[0]
        inchar = mut[-1]

        position = int(mut[1:-1])
        # Char doesnt match
        # assert ref[position-1] == outchar
        if ref[position-1] != outchar:
            log_fn("ERROR")
            log_fn(mut)
            log_fn(mutations)
            log_fn(position)
            log_fn(f"{ref[position-1]},{outchar}")
            exit(1)
        
        # this is a truncation
        if inchar == "-":
            ref = ref[:position-1] 
        else:
            ref = ref[:position-1] + inchar + ref[position:]
    

    return ref.upper()


rows = []
for mutation in benign_muts:
    rows.append({
        "ID": mutation,
        "Sequence": mutate_protein(ref, mutation.split("_")),    
        "clinical_sig": False   
    })


for mutation in path_muts:
    rows.append({
        "ID": mutation,
        "Sequence": mutate_protein(ref, mutation.split("_")),     
        "clinical_sig": True    
    })
# Create the DataFrame
df = pd.DataFrame(rows)

df.to_csv("prion_input.csv", index=False)

