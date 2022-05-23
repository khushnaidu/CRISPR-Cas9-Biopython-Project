from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62

with open('gene.fna') as f:
    firstSeq = f.readlines()
firstSeq.pop(0)
for i in range(len(firstSeq)):
    firstSeq[i] = firstSeq[i].rstrip("\n")

print(firstSeq)