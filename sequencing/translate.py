CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
COMPLEMENTS = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}


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

def dna_to_protein(dna):
    """converts a single DNA strand in a single reading frame (starting at 0) to an amino acid chain (protein)

    Args:
        dna (string): string rep of DNA strand

    Returns:
        string: string rep of a protein as a chain of amino acids
    """
    
    protein = ""
    for i in range(0, len(dna) - 3, 3):
        protein += CODON_TABLE[dna[i:i+3]]
    return protein
    

def dna_to_all_proteins(dna):
    """converts a DNA strand from a database (for example) into all proteins it could encode from all reading frames and
    the reverse complement. 

    Args:
        dna (string): string rep of DNA strand

    Returns:
        List[string]: list of string reps for all proteins that can be made as a chain of amino acids
    """
    
    proteins = []
    reverse_c = reverse_complement(dna)
    for i in range(0, 3):
        proteins.append(dna_to_protein(dna[i:]))
        proteins.append(dna_to_protein(reverse_c[i:]))

    
    return proteins

def main():
    ## Example test case
    dna = read_fasta("brca.fasta")
    # print(dna[:100])
    print(dna_to_all_proteins(dna[:100]))
    
if __name__ == "__main__":
    main()