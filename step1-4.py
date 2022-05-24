from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62
genome_list = []

# Opens the file that has the test cases in it.
# Reads in each genome into its own list, per line
# BROKEN
with open('gene.fna') as f:
    line = f.readline().rstrip("\n")
    all_genome = []
    counter = 0
    while line:
        if(line != ""):
            ## This is the line we want to record as the genome
            if line.find(">") == -1:
                all_genome.append(line)
                line = f.readline().rstrip("\n")
            else:
                counter += 1
                g = ''.join(all_genome)
                genome_list.append(g)
                all_genome = []
                line = f.readline().rstrip("\n")
    g = ''.join(all_genome)
    genome_list.append(g)

# Get rid of extra entries that are ''
genome_list = [x for x in genome_list if x]
# Now gpo through each list genome 
# and run algorthim to pick out sequnces of 20 nucleotides preceding a PAM.
all_genomes_RGNA = []
for a_genome in genome_list:
    all_GRNA_sequnces_current = []
    find_index = 0
    while find_index != -1:
        find_index = a_genome.find("NGG")
        if (find_index > 20):
            all_GRNA_sequnces_current.append(a_genome[find_index-20:find_index])
            a_genome = a_genome[find_index+3:]
    all_genomes_RGNA.append(all_GRNA_sequnces_current)

# Convert to complement
for seq in range(len(all_genomes_RGNA)):
    for gRNA in range(len(all_genomes_RGNA[seq])):
        new_string = ""
        for n in all_genomes_RGNA[seq][gRNA]:
            if(n == "A"):
                new_string += "T"
            if(n == "T"):
                new_string += "A"
            if(n == "G"):
                new_string += "C"
            if(n == "C"):
                new_string += "G"
        all_genomes_RGNA[seq][gRNA] = new_string
    new_string = ""

print(all_genomes_RGNA)
                        
