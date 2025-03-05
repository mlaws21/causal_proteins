import sys

def main():
    
    if len(sys.argv) < 2:
        print("Usage: python prion_quick_mut.py [xNx]")
        exit(1)
        
    full_mut = sys.argv[1]
    
    base = "MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"
    ref = base.lower()
    
    for mut in full_mut.split(","):

        outchar = mut[0]
        inchar = mut[-1]
        
        position = int(mut[1:-1])
        
        # print(position)
        
        assert ref[position-1] == outchar
    
        ref = ref[:position-1] + inchar + ref[position:]
    

    print(ref.upper())
    # print(base)
    
    

if __name__ == "__main__":
    main()
