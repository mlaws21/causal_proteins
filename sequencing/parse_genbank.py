def find_protein_seqs_genbank(filename):
    """Reads a .fasta file

    Args:
        filename (str): File to read

    Returns:
        str: String containing all the bases in the file
    """
    f = open(filename, "r")
    proteins = []
    curr = ""
    isoforms = []
    for line in f:
        ls = line.strip()
        # print(ls)
        
        if "note" in ls:
            isoforms.append(ls)
        
        if curr != "":
            curr += ls
            
        if "\"" in ls:
            if curr != "":
                proteins.append(curr.split("\"")[-2])
            curr = ""
        if "translation" in ls:
            curr += ls
    
    f.close()
    return proteins, isoforms

def find_cds(filename):
    f = open(filename, "r")
    cds = []
    alive = (False, False)
    curr = ""
    for line in f:
        ls = line.strip()

        if "CDS             " in ls:
            alive = (True, True)
            
        if alive[0] and "translation" in ls:
            alive = False, True
        
        if not alive[0] and "\"" in ls:
            curr += ls
            cds.append(curr)
            alive = (False, False)
            curr = ""
        if alive[0] or alive[1]:
            curr += ls
            
    return cds
def main():
    
    # proteins, isoforms = parse_genbank("BRCA_datasets/sequence.gb")
    # print(proteins)
    
    # cds = find_cds("BRCA_datasets/sequence.gb")
    
    # cds = find_cds("prion_datasets/sequence.gb")
    
    proteins, _ = find_protein_seqs_genbank("prion_datasets/sequence.gb")
    
    # print(cds[0])
    
    print([len(x) for x in proteins])
    
    #transcript variant 3
    print(proteins[0])
    
    print(proteins[0][0:104] + "L" + proteins[0][105:])
    
    
    # for i in proteins:
    #     if len(i) == 1863:
    #         print(i)
    #         break

if __name__ == "__main__":
    main()