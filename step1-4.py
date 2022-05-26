from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import json

matrix = matlist.blosum62

with open('config.json', 'r') as json_file:
	json_load = json.load(json_file)

# print(json_load)

# Opens the file that has the test cases in it.
# Reads in each genome into its own list, per line
def read_in_data():
    genome_list = []
    with open(json_load['file_name']) as f:
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
        genome_list = [x for x in genome_list if x]
        return genome_list

def check_for_target_sites(all_genomes_gRNA):
    target_seqeunce_sites_count = 0
    for seq in range(len(all_genomes_gRNA)):
        for gRNA in range(len(all_genomes_gRNA[seq])):
            if all_genomes_gRNA[seq][gRNA].find(json_load["target_sequence"]) != -1:
                target_seqeunce_sites_count += 1
    return target_seqeunce_sites_count


# Now go through each list genome 
# and run algorthim to pick out sequnces of 20 nucleotides preceding a PAM.
def find_viable_gRNA(genome_list):
    all_genomes_gRNA = []
    for a_genome in genome_list:
        all_GRNA_sequnces_current = []
        counter1 = 0
        counter2 = 1
        for i in a_genome:
            if(counter2 < len(a_genome) and a_genome[counter1] == "G" and a_genome[counter2] == "G" and counter1 >= 20):
                all_GRNA_sequnces_current.append(a_genome[counter1-20:counter2])
            if(counter2 == len(a_genome)):
                break
            counter1+=1
            counter2+=1
        all_genomes_gRNA.append(all_GRNA_sequnces_current)    
    return all_genomes_gRNA

def convert_to_RNA(all_genomes_gRNA):
    for seq in range(len(all_genomes_gRNA)):
        for gRNA in range(len(all_genomes_gRNA[seq])):
            new_string = ""
            for n in all_genomes_gRNA[seq][gRNA]:
                if(n == "A"):
                    new_string += ""
                if(n == "T"):
                    new_string += "A"
                if(n == "G"):
                    new_string += "C"
                if(n == "C"):
                    new_string += "G"
            all_genomes_gRNA[seq][gRNA] = new_string
        new_string = ""
    return all_genomes_gRNA

# Convert to complement RNA
def convert_to_complement_dna(all_genomes_gRNA):
    for seq in range(len(all_genomes_gRNA)):
        for gRNA in range(len(all_genomes_gRNA[seq])):
            new_string = ""
            for n in all_genomes_gRNA[seq][gRNA]:
                if(n == "A"):
                    new_string += "T"
                if(n == "T"):
                    new_string += "A"
                if(n == "G"):
                    new_string += "C"
                if(n == "C"):
                    new_string += "G"
            all_genomes_gRNA[seq][gRNA] = new_string
        new_string = ""
    return all_genomes_gRNA

def scoring_gRNA(all_genomes_gRNA):
    gRNA_score = [] #array to store score
    for seq in range(len(all_genomes_gRNA)):
        for gRNA in range(len(all_genomes_gRNA[seq])):
            TCscore = 0
            GCscore = 0
            CumulativeScore = 0
            index = 0
            for n in all_genomes_gRNA[seq][gRNA]:
                index = index + 1
                if (index >= 15):
                    if (n == 'T' or n == 'C'):
                        TCscore+= 1  #Assigning score based on number of C's and T's adjacent to PAM
                if (n == 'G' or n == 'C'):
                    GCscore+= 2  # Assigning a score for each G and C found (weighting)
            CumulativeScore = GCscore - TCscore # Higher the score, the better the quality of gRNA
            gRNA_score.append(CumulativeScore)
    return gRNA_score
                        
# for seq in range(len(all_genomes_gRNA)):
#     for gRNA in range(len(all_genomes_gRNA[seq])):
#         if(all_genomes_gRNA[seq][gRNA].find("AAAA") != -1):
#             print("Found AAAA")
#             all_genomes_gRNA[seq][gRNA] = all_genomes_gRNA[seq][gRNA][find_index:find_index+4]
#         if(all_genomes_gRNA[seq][gRNA].find("CCCC") != -1):
#             print("Found CCCC")
#             all_genomes_gRNA[seq][gRNA] = all_genomes_gRNA[seq][gRNA][find_index:find_index+4]
#         if(all_genomes_gRNA[seq][gRNA].find("GGGG") != -1):
#             print("Found GGGG")
#             all_genomes_gRNA[seq][gRNA] = all_genomes_gRNA[seq][gRNA][find_index:find_index+4]
#         if(all_genomes_gRNA[seq][gRNA].find("UUUU") != -1):
#             print("Found UUUU")
#             all_genomes_gRNA[seq][gRNA] = all_genomes_gRNA[seq][gRNA][find_index:find_index+4]

def write_to_output_file(all_genomes_gRNA,gRNA_score):
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
    for seq in range(len(all_genomes_gRNA)):
        for gRNA in range(len(all_genomes_gRNA[seq])):
            f.write(all_genomes_gRNA[seq][gRNA])
            f.write("\t \t \t  ")
            f.write(str(gRNA_score[gRNA]))
            f.write("\n")

if __name__ == "__main__":
    text = input("Would you like to do a general anaylsis for validaity of gRNA in the genome file or look for a specific target sequence? 1 or 2.")
    genomes_split_by_seq = read_in_data()
    potential_DNA_per_seq = find_viable_gRNA(genomes_split_by_seq)
    if(text == "1"):
        ## This list contains a nested [[]]. Each inner list is all gRNA sequences per genome. 
        final_genome_list_gRNA = convert_to_RNA(potential_DNA_per_seq)
        scores_gRNA = scoring_gRNA(final_genome_list_gRNA)
        write_to_output_file(final_genome_list_gRNA,scores_gRNA)
    if(text == "2"):
        converted =convert_to_complement_dna(potential_DNA_per_seq)
        x = check_for_target_sites(converted)
        f = open("target_sequence_results.txt", "w")
        f.write("Total number of valid target sequences were found regarding " + json_load["diesease_name"] + ": " + str(x))
        

