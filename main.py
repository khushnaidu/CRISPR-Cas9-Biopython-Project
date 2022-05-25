from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

matrix = matlist.blosum62
genome_list = []

# Opens the file that has the test cases in it.
# Reads in each genome into its own list, per line
# BROKEN
with open('/Users/biswajeetnaidu/PycharmProjects/CRISPR-Cas9 Project/venv/lib/test.fna.txt') as f:
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
all_genomes_GRNA = []
for a_genome in genome_list:
    all_GRNA_sequnces_current = []
    find_index = 0
    while find_index != -1:
        find_index = a_genome.find("CGG")
        if (find_index > 20):
            all_GRNA_sequnces_current.append(a_genome[find_index-20:find_index])
            a_genome = a_genome[find_index+3:]
    all_genomes_GRNA.append(all_GRNA_sequnces_current)

# Convert to complement
for seq in range(len(all_genomes_GRNA)):
    for gRNA in range(len(all_genomes_GRNA[seq])):
        new_string = ""
        for n in all_genomes_GRNA[seq][gRNA]:
            if(n == "A"):
                new_string += "T"
            if(n == "T"):
                new_string += "A"
            if(n == "G"):
                new_string += "C"
            if(n == "C"):
                new_string += "G"
        all_genomes_GRNA[seq][gRNA] = new_string
    new_string = ""

print(all_genomes_GRNA)

gRNA_score = [] #array to store score
for seq in range(len(all_genomes_GRNA)):
    for gRNA in range(len(all_genomes_GRNA[seq])):
        TCscore = 0
        GCscore = 0
        CumulativeScore = 0
        index = 0
        for n in all_genomes_GRNA[seq][gRNA]:
            index = index + 1
            if (index >= 15):
                if (n == 'T' or n == 'C'):
                    TCscore+= 1  #Assigning score based on number of C's and T's adjacent to PAM
            if (n == 'G' or n == 'C'):
                GCscore+= 2  # Assigning a score for each G and C found (weighting)
        CumulativeScore = GCscore + TCscore # Higher the score, the worse the quality of gRNA
        gRNA_score.append(CumulativeScore)

print(gRNA_score)

f = open("output.txt", "w")
f.write("The metrics list we used to make this score were: \n")
f.write("\n")
f.write("-The counting of the Ts and Cs adjacent to the PAM site(Because the nucleotides adjacent to the PAM site contain significantly lower counts of Cs and Ts\n")
f.write("-A measure of higher GC content in the whole sequence(Functional gRNAs show lower counts of GC content\n")
f.write("\n")
f.write("Low Goodness Score = Good Quality gRNA\n")
f.write("\n")
f.write("---------------------------------------------------------------------\n")
f.write("\n")
f.write("      SEQUENCE      \t \tGOODNESS SCORE \n")
f.write("\n")
for seq in range(len(all_genomes_GRNA)):
    for gRNA in range(len(all_genomes_GRNA[seq])):
        f.write(all_genomes_GRNA[seq][gRNA])
        f.write("\t \t \t  ")
        f.write(str(gRNA_score[gRNA]))
        f.write("\n")



