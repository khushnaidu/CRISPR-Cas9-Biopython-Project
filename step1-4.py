from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62
genome_list = []

# Opens the file that has the test cases in it.
# Reads in each genome into its own list, per line
with open('gene.fna') as f:
    line = f.readline()
    all_genomes = []
    while line:
        if("<" not in line):
            all_genomes.append(line)
            line = f.readline().rstrip("\n")
        else:
            line = f.readline().rstrip("\n")

print(len(all_genomes))
# Now gpo through each list genome 
# and run algorthim to pick out sequnces of 20 nucleotides preceding a PAM.
all_genomes_RGNA = []
for a_genome in all_genomes:
    # print(len(a_genome))
    all_GRNA_sequnces_current = []
    current_genome = ''.join(a_genome)
    # print(len(current_genome))
    # print(type(current_genome))
    find_index = 0
    while find_index is not -1:
        find_index = current_genome.find("NGG")
        # print(find_index)
        if (find_index > 20):
            all_GRNA_sequnces_current.append(current_genome[find_index-20:find_index])
            current_genomoe = current_genome[find_index+3:]
            # print(len(current_genome))
    #     else:
    #         # print("need to fix")
    # print("moving on")

print(all_genomes_RGNA)
                        
