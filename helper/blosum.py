from Bio.Align import substitution_matrices

def blosum_scores(ref, seqs, log_fn):
    out = []
    for i, x in enumerate(seqs):
        if i % (len(seqs) // 10) == 0:
            log_fn(f"{(i / len(seqs)) * 100:.0f}% complete")
            
        out.append(blosum_align(ref, x))
    return out

def blosum_align(s1, s2, gap_pen=6):
    """function that prints out the best global alignment of two sequences

    Args:
        s1 (string): element to align
        s2 (_type_): element to align
    """
    
    blosum90 = substitution_matrices.load("BLOSUM90")
    
    # s1 is rows s2 is columns
    n = len(s1) + 1
    m = len(s2) + 1
    
    s1 = "-" + s1
    s2 = "-" + s2
    
    # first element of tuple is for score second is for back pointer
    dp = [[0] * n for _ in range(m)]
    val = 0
    
    for i in range(n):
        dp[0][i] = val
        val -= gap_pen
    
    val = 0
    for i in range(m):
        dp[i][0] = val
        val -= gap_pen
    
    for i in range(1, m):
        for j in range(1, n):
            
            ind = blosum90[s1[j], s2[i]]
            diag = dp[i-1][j-1] + ind
            up = dp[i-1][j] - gap_pen
            left = dp[i][j-1] - gap_pen
            
            dp[i][j] = max(diag, up, left)
            
    return dp[-1][-1]


def main():
    x = "MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSSQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"
    y = "MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNIAIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"
    print(blosum_align(x, y))

if __name__ == "__main__":
    main()