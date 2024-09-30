from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import random
random.seed(2)
# global constants
COMPLEMENTS = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]

def read_fasta(filename):
    """Reads a .fasta file

    Args:
        filename (str): File to read

    Returns:
        str: String containing all the bases in the file
    """
    f = open(filename, "r")
    out = ""
    for line in f:
        if line[0] != ">":
            out += line.strip()
    
    return out

def reverse_complement(strand):
    """Converts DNA coding strand to the reverse complement

    Args:
        genome (str): DNA coding strand

    Returns:
        str: reverse complement of input strand
    """

    out = ""
    for i in strand[::-1]:
        out += COMPLEMENTS[i]
    
    return out

def count_nucleotides(strand):
    """Counts the frequency of each nucleotide in a strand

    Args:
        strand (str): String of bases

    Returns:
        dict: Frequency Dictionary
    """
    return Counter(strand)

def count_dinucleotides(strand, loop=False):
    """Counts the frequency of each dinucleotide in a strand

    Args:
        strand (str): String of bases
        loop (bool): True if the first and last base count as a dinucleotide

    Returns:
        dict: Frequency Dictionary
    """
    mydict = {}
    for i in range(len(strand) - 1):
        dinucleotide = strand[i:i+2]

        mydict[dinucleotide] = mydict.get(dinucleotide, 0) + 1
        
    if loop:
        dinucleotide = strand[-1] + strand[0]
        mydict[dinucleotide] = mydict.get(dinucleotide, 0) + 1
    
    return mydict

def find_dinucleotide_ratios(strand, loop=False):
    """Calculates all the dinucleotide ratios and prints them, also prints the max.

    Args:
        strand (str): String of bases
        loop (bool): True if the first and last base count as a dinucleotide

    """
    n = len(strand)
    dinuc_dict = count_dinucleotides(strand, loop=loop)
    nuc_dict = count_nucleotides(strand)
    
    ratios = []
    for k, v in dinuc_dict.items():
        n1, n2 = k[0], k[1]
        p_n1n2 = v / n
        p_n1 = nuc_dict[n1] / n
        p_n2 = nuc_dict[n2] / n
        ratio = p_n1n2 / (p_n1 * p_n2)
        ratios.append((k, ratio))
        print(f"{k}: {ratio}")
        
    print(max(ratios, key=lambda x:x[1]))
    
    

def find_orfs_single_frame(strand, orfs):
    """Finds all ORFs for a single strand at a single reading frame

    Args:
        strand (str): String of bases
        orfs (list): Running list of ORFs for a single frame just pass an empty list

    Returns:
        List[(int, int)]: returns a list of start and end codon positions for each ORF
    """
    # this can be used as a flag, when start is < 0 it denotes we are not in an orf
    # if start > 0 then start denotes the start of the current orf.
    start = -1
    n = len(strand)
    for c in range(0, n, 3):
        codon = strand[c:c+3]
        
        
        if codon == START_CODON and start == -1:
                start = c // 3
                
                
        if codon in STOP_CODONS and start != -1:
            if c // 3 - start > 150:
                orfs.append((start, c // 3))
            # resets the counting and begins looking for new start codon
            start = -1

            
    return orfs

def find_all_orfs(strand):
    """Finds the ORFs for all reading frames and for the reverse complement.

    Args:
        strand (str): String of bases

    Returns:
        List[(int, int)]: returns a list of start and end codon positions for each ORF
    """
    orfs = []
    reverse_c = reverse_complement(strand)
    for i in range(0, 3):
        orfs = find_orfs_single_frame(strand[i:], orfs)
        orfs = find_orfs_single_frame(reverse_c[i:], orfs)
        
    return orfs
            
def plot_orf_lens(orfs):
    """Plots the lengths of the ORFs

    Args:
        orfs List[(int, int)]: A list of start and end codon positions for each ORF
    """
    
    # gets the lengths out
    lens = [end - start for start, end in orfs]
    frequencies = Counter(lens)

    # generating graph
    plt.scatter(frequencies.keys(), frequencies.values(), s=10)
        
    plt.xlabel("Number of Codons")
    plt.ylabel("Frequency")
    plt.title("Distribution of ORF lengths")
    plt.xscale
    
    plt.show()
    
    
def plot_orf_vs_random(orfs, rand_orfs):
    """Plots the lengths of genetic open reading frames against the lengths if it was generated totally randomly

    Args:
        orfs List[(int, int)]: A list of start and end codon positions for each ORF from genetic data
        rand_orfs List[(int, int)]: A list of start and end codon positions for each ORF from random data
    """

    lens = [end - start for start, end in orfs]
    frequencies = Counter(lens)
    
    rand_lens = [end - start for start, end in rand_orfs]
    rand_frequencies = Counter(rand_lens)

    plt.scatter(frequencies.keys(), frequencies.values(), s=10, color="blue",label="DNA")
    plt.scatter(rand_frequencies.keys(), rand_frequencies.values(), s=10, color="red", label="Random")
    
        
    plt.xlabel("Number of Codons")
    plt.ylabel("Frequency")
    plt.title("Distribution of ORF lengths in Random and Genetic Data")
    plt.legend()
    
    plt.show()
        
        
def slow_count_nucleotides(strand, isN=False):
    """Counts the frequency of each nucleotide in a strand using .count()

    Args:
        strand (str): String of bases


    Returns:
        dict: Frequency Dictionary
    """
    
    possible = ["A", "G", "C", "T"]
    if isN:
        possible.append("N")
    
    for base in possible:
        print(f"{base}: {strand.count(base)}")
        

def slow_count_dinucleotides(strand, isN=False):
    """Counts the frequency of each dinucleotide in a strand using .count()

    Args:
        strand (str): String of bases
        isN (bool): True if their are Ns in the sequence

    Returns:
        dict: Frequency Dictionary
    """
    bases = ["A", "G", "C", "T"]
    if isN:
        bases.append("N")
    
    
    possible = []
    for i in bases:
        for j in bases:
            possible.append(i + j)
            
    di_strand = [strand[x:x+2] for x in range(len(strand) - 1)]
    for base in possible:
        print(f"{base}: {di_strand.count(base)}")


def find_length_confidence(rand_orfs, quantile=0.999):
    """_summary_

    Args:
        rand_orfs List[(int, int)]: A list of start and end codon positions for each ORF from random data
        quantile (float, optional): Percentage of random occurances to account for. Defaults to 0.999.

    Returns:
        int: Number x such that any length greater than x emperically occurs randomly with only 1 - quantile probablility
    """
    rand_lens = [end - start for start, end in rand_orfs]
    
    rand_lens.sort()
    return rand_lens[int(quantile * len(rand_lens)) + 1]


    
def main():
    """Main method for hw 1, uncomment the lines you need to generate specific output.
    """

    ### Reading .fasta ###
    strand = read_fasta("brca.fasta")
    # print(strand[983:985])
    mutated = strand[:983] + strand[985:]
    print(len(mutated))
    
    orfs = find_orfs_single_frame(strand, []) 
    print([y-x for x, y in orfs])
    
    orfs = find_orfs_single_frame(strand[1:], []) 
    print([y-x for x, y in orfs])
    
    orfs = find_orfs_single_frame(strand[2:], []) 
    print([y-x for x, y in orfs])

    
    # human_strand = read_fasta("humanChr19.fasta")
    
    ### Generating random DNA strand ###
    # rand_strand = "".join([["A", "C", "G", "T"][random.randint(0, 3)] for _ in range(len(strand))])

    ### Generate reverse complement ###
    #reverse_c = reverse_complement(strand)
    
    ### Counting nucleotides and dinucleotides and dinucleotides ratios ###
    # print(count_nucleotides(strand))
    # print(count_dinucleotides(strand))
    # find_dinucleotide_ratios(strand)
    
    ### Sorting dinucleotide counts ###
    # print(sorted([(k, v) for k, v in count_dinucleotides(strand).items()], key=lambda x:x[1]))
    # print(sorted([(k, v) for k, v in count_dinucleotides(human_strand).items()], key=lambda x:x[1]))
    
    ### Finding ORFs ###
    # orfs = find_all_orfs(strand)
    # random_orfs = find_all_orfs(rand_strand)
    
    ### Plotting ORF Lengths and determining confidence ###
    # plot_orf_lens(orfs)
    # plot_orf_vs_random(orfs, random_orfs)
    # print(find_length_confidence(random_orfs))
    
    ### .count() functions ###
    # slow_count_nucleotides(strand)
    # slow_count_dinucleotides(strand)
    

if __name__ == "__main__":
    main()