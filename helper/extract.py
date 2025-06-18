import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import random



def extract_full_gene(chr_fasta, start, end):
    with open(chr_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # if record.id == "chr17":
            # Extract the sequence (convert to 0-based indexing for Python slicing)
            full = record.seq
            extracted_sequence = record.seq[start-1:end]
            # reverse_complement = extracted_sequence.reverse_complement()
            # print(f"Extracted sequence ({start}-{end}):")
        
    return extracted_sequence

def parse_hgvs(hgvs):
    
    muts = hgvs.split(",")
    
    out = []
    for mut in muts:
        chrom, start_pos, change = mut.split(":")
        
        ref, alt = change.split(">")
        
        spl = start_pos.split(".")
        
        start = spl[-1]

        assert spl[0] != "c" # FIXME

        out.append(
            {
                "ref": ref,
                "alt": alt,
                "start": int(start)
            }
        )
        
    
    return out
    
def extract_exon_locations(filename):
    """Extracts exon locations (joined intervals) from a GenBank file.

    Args:
        filename (string): Path to the GenBank (.gb) file.

    Returns:
        list[list[tuple[int, int]]]: A list of exon ranges for each CDS or exon feature.
                                     Each feature is represented as a list of (start, end) tuples.
    """
    exon_locations = {}

    # Parse the GenBank file
    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                # Look for CDS or exon features
                if feature.type in ["CDS", "exon"]:
                    if "join" in str(feature.location):  # Handles joined intervals
                        # Extract all parts of the join
                        exon_ranges = [
                            (int(part.start+1), int(part.end))
                            for part in feature.location.parts
                        ]
                        # print(feature.qualifiers.get("protein_id", []))
                    else:
                        # Single interval
                        exon_ranges = [(int(feature.location.start), int(feature.location.end))]

                    
                    exon_locations[feature.qualifiers.get("protein_id", [])[0]] = exon_ranges[::-1]
    
    return exon_locations


def reconstruct_cds(seq, exon_ranges):
    """Reconstructs the coding sequence (CDS) from the genomic sequence and exon ranges.

    Args:
        seq (str): The full genomic sequence as a string.
        exon_ranges (list[tuple[int, int]]): List of exon ranges (start, end) as 0-based coordinates.

    Returns:
        str: The reconstructed CDS sequence.
    """
    
    if seq is None:
        return None
    cds = ""
    for start, end in exon_ranges:
        
        # FIXME this is start-1 if we are dealing with BRCA but idk why it isnt with prion
        cds += seq[start:end]  # Extract each exon and concatenate
    return cds


def transcribe_and_translate(cds):
    """Transcribes a DNA CDS to RNA and translates it to a protein sequence.

    Args:
        cds (str): Coding DNA sequence (CDS).

    Returns:
        tuple: Transcribed mRNA sequence and translated protein sequence.
    """
    
    if cds is None:
        return None
    # Convert to Seq object
    dna_seq = Seq(cds)
    
    # print(cds)
    
    dna_seq = dna_seq
    
    # print(dna_seq)

    # Transcribe (DNA → RNA)
    rna_seq = dna_seq.transcribe()
    # print(f"Transcribed mRNA: {rna_seq}")

    # Translate (RNA → Protein)
    protein_seq = rna_seq.translate(to_stop=True)  # Stop translation at stop codon
    # print(f"Translated Protein: {protein_seq}")

    return str(protein_seq)



# TODO look at HGVS format

# want mutation passed in at HGVS format
def mutate_gene(gene, mutations, diagnosis, exons, isoform, expected_len, start):
    seqs = []
    
    
    # print("len(unique(muts))", len(set(mutations)), len(diagnosis))

    for mut, diag in zip(mutations, diagnosis):
        
        if mut == "":  
            ref_protein = transcribe_and_translate(reconstruct_cds(gene, exons[isoform]))
            
            seqs.append(("", ref_protein, diag))
            continue

        hgvs = parse_hgvs(mut)
        # we can remove these one day when i am not lazy
        if len(hgvs) > 1:
            continue
        
        
        row = hgvs[0]
        
        
        if len(row["ref"]) != len(row["alt"]):
            continue
        
        if row["start"] <= start:
            continue
        
        mut_end = row["start"] + len(row["ref"]) - start
        mut_start = row["start"] - start
        
        mutated = gene[:mut_start] + row["alt"] + gene[mut_end:]

        # print(len(mutated))
        
        mut_protein = transcribe_and_translate(reconstruct_cds(mutated, exons[isoform]))
        
        # truncation
        if len(mut_protein) != expected_len:
            continue
        
        seqs.append((mut, mut_protein, diag))

        
    return seqs

def parse_xnx(xnx, base):

    ref = base.lower()
    
    for mut in xnx.split("_"):

        outchar = mut[0]
        inchar = mut[-1]
        
        position = int(mut[1:-1])
        assert ref[position-1] == outchar
        
        if inchar == "-":
            # print("uere")
            ref = ref[:position-1] 
        else:
            ref = ref[:position-1] + inchar + ref[position:]
    

    return ref.upper()

def mutate_gene_xnx(gene, mutations, diagnosis, exons, isoform, expected_len, start):
    seqs = []
    
    
    # print("len(unique(muts))", len(set(mutations)), len(diagnosis))

    for mut, diag in zip(mutations, diagnosis):
        
        ref_protein = transcribe_and_translate(reconstruct_cds(gene, exons[isoform]))
        if mut == "":  
            
            seqs.append(("", ref_protein, diag))
            continue

        mut = mut.replace("*", "-")
        mut_protein = parse_xnx(mut, ref_protein)
        
        # truncation
        if len(mut_protein) != expected_len:
            continue
        
        seqs.append((mut, mut_protein, diag))

        
    return seqs
    
    
def brca_driver():
    
    start = 43044295
    end = 43170327
    isoform = "NP_001394531.1"
    expected_len = 1863 
    
    gene = extract_full_gene("hg38_chr17.fasta", start, end)

    mutation_data = pd.read_csv("brca1_data.csv")
    
    hgvs_muts = mutation_data["Genomic_Coordinate_hg38"]
    diagnosis = mutation_data["Clinical_significance_ENIGMA"] == "Pathogenic"
    
    exons = extract_exon_locations("brca1_hg38_chr17.gb")
    
    gene_muts = mutate_gene(gene, hgvs_muts, diagnosis, exons, isoform, expected_len, start)
    
    
    # print((set(gene_muts)))
    
    
    # [print(x) for x in mutated_proteins]
    # print(len(mutated_proteins), len(set(mutated_proteins)))
    benign = set([(x, y) for x, y, z in gene_muts if not z])
    cancer = set([(x, y) for x, y, z in gene_muts if z])
    
    
    # TODO this is so ugly
    unique_benign = []
    seen = set()
    for i in benign:
        if i[1] not in seen:
            unique_benign.append(i)
            seen.add(i[1])

    unique_cancer = []
    seen = set() 
    for i in cancer:
        if i[1] not in seen:
            unique_cancer.append(i)
            seen.add(i[1])
                
    
    
    # benign = set([y for x, y, z in gene_muts if not z])
    # cancer = set([y for x, y, z in gene_muts if z])
    
    print(len(unique_benign), len(unique_cancer))

    
    df = generate_synth(unique_benign, unique_cancer)
    df.to_csv('brca.csv', index=False)
    
    # print(benign, cancer)
    
    
# first go to NCBI and download hg38 chromosome gene is on
# search the gene in NCBI
# go to the genomic context section and grab location start and end and enter them below
# go to "Genomic regions, transcripts, and products"
# --> click Go to nucleotide: GenBank
# look through .gb file for isoform you like
# get the len of that isoform and enter below too

def prion_driver():
    
    start = 4686456
    end = 4701588
    isoform = "NP_001073590.1"
    expected_len = 253 
    
    gene = extract_full_gene("hg38_chr20.fasta", start, end)

    # print(len(gene))
    
    # mutation_data = pd.read_csv("brca1_data.csv")
    
    # hgvs_muts = mutation_data["Genomic_Coordinate_hg38"]
    # diagnosis = mutation_data["Clinical_significance_ENIGMA"] == "Pathogenic"
    
    exons = extract_exon_locations("prnp_hg38_chr20.gb")
    
    mutation_data = pd.read_csv("prion_muts.csv")
    
    muts = mutation_data["mutation"]
    diag = mutation_data["diagnosis"] == "pathogenic"
    
    
    
    # hgvs_muts = ["", "chr20:4688936:C>T", "chr20:4699570:C>T"]
    # diagnosis = [False, True, False]
    # print(exons)
    
    gene_muts = mutate_gene_xnx(gene, muts, diag, exons, isoform, expected_len, start)
    
    
    # print((set(gene_muts)))
    
    
    # [print(x) for x in mutated_proteins]
    # print(len(mutated_proteins), len(set(mutated_proteins)))
    benign = set([(x, y) for x, y, z in gene_muts if not z])
    cancer = set([(x, y) for x, y, z in gene_muts if z])
    
    
    # TODO this is so ugly
    unique_benign = []
    seen = set()
    for i in benign:
        if i[1] not in seen:
            unique_benign.append(i)
            seen.add(i[1])

    unique_cancer = []
    seen = set() 
    for i in cancer:
        if i[1] not in seen:
            unique_cancer.append(i)
            seen.add(i[1])
    
    df = generate_synth(unique_benign, unique_cancer)
    df.to_csv('prion_uniprot.csv', index=False)