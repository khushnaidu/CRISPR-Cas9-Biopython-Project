from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62
genome_list = []

# Opens the file that has the test cases in it.
# Reads in each genome into its own list, per line
with open('test.fna') as f:
    line = f.readline()
    all_genome = ""
    while line:
        if line.find(">") != -1:
            genome_list.append(all_genome)
            all_genome = ""
            line = f.readline()
        else:
            line = f.readline().rstrip("\n")
            all_genome += line
    genome_list.append(all_genome)
genome_list = [x for x in genome_list if x]

# Now gpo through each list genome 
# and run algorthim to pick out sequnces of 20 nucleotides preceding a PAM.
all_genomes_RGNA = []
for a_genome in genome_list:
    all_GRNA_sequnces_current = []
    find_index = 0
    while find_index is not -1:
        find_index = a_genome.find("NGG")
        print(find_index)
        if (find_index > 20):
            all_GRNA_sequnces_current.append(a_genome[find_index-20:find_index])
            a_genome = a_genome[find_index+3:]
            print(a_genome)
            print(find_index)
    all_genomes_RGNA.append(all_GRNA_sequnces_current)

print(all_genomes_RGNA)
                        
