# Python version: Python 3.8.5
# Author: Alba Burillo Navarro
# OS: Windows 10
# Program name: IDENTIFYING NEISSERIA SPECIES
# Description: The aim of this script is to be able to differentiate between the species of the genus Neisseria,
# starting from unknown sequences of the rplF (50S ribosomal protein L6) [3] and rpmH (50S ribosomal protein L34)
# [4] genes, obtaining the molecular weight of the proteins.

# BIOLOGICAL PROBLEM: Population studies have used genomic techniques and confirmed that the main bacterial genera
# isolated in the oral cavity are: Streptococcus, Actinomyces, Veillonella, Fusobacterium, Porphyromonas, Prevotella,
# Treponema, Neisseria, Haemophilus, Eubacteria, Lactobacterium, Capnophycum, Pseudomonas and Propionibacterium. [1]
# Among them we can find different pathogenic species, within the genus Neisseria we have two of great importance,
# Neisseria meningitidis and Neisseria gonorrhoeae.

# In the MALDI-TOF Mass Spectrometry technique, which generates unique protein mass spectrum for each species based
# on their ribosomal proteins, there are problems for the differentiation of these two species from the rest of the
# non-pathogenic species of the Neisseria spp [2]. Due to the great similarity of their spectrum, and as a consequence
# to the great similarity of the species to each other.


# We need to install different modules:
# pip install biopython
from Bio import SeqIO
# pip install scipy
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

# We have a FASTA file (L34.fasta) with the DNA sequence of the rpmH gene,
# which encodes the 50S ribosomal protein L34 in different microorganisms.
# (Including different species of Neisseria)

# Let's take a look at the original DNA file:
for seq_record in SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34.fasta",
                              "fasta"):
    # We can see the ID of the sequence and the microorganism to which it belongs
    print(seq_record.description)
    # We can see the DNA sequence
    print(repr(seq_record.seq))
    # We can see the length of the sequence
    print(len(seq_record))

# We transform the DNA sequence into RNA
file_in = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34.fasta"
file_out = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34_RNA.fasta"
with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        # We use the .transcribe() argument to obtain the RNA sequence
        seq_record.seq = seq_record.seq.transcribe()
        # write new fasta file
        r = SeqIO.write(seq_record, f_out, 'fasta')

# Let's take a look at the L34_RNA file
for seq_record in SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34_RNA.fasta",
                              "fasta"):
    # We can see the ID of the sequence and the microorganism to which it belongs
    print(seq_record.description)
    # We can see the RNA sequence
    print(repr(seq_record.seq))
    # We can see the length of the sequence
    print(len(seq_record))

# We obtain the protein sequence
file_in = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34.fasta"
file_out = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34_PROT.fasta"
with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        seq_record.seq = seq_record.seq.translate(to_stop=True)
        # write new fasta file
        r = SeqIO.write(seq_record, f_out, 'fasta')

# Let's take a look at the L34_PROT file
for seq_record in SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34_PROT.fasta",
                              "fasta"):
    print(seq_record.description)
    # We can see the PROTEIN sequence
    print(repr(seq_record.seq))
    print(len(seq_record))


# We create a function to obtain the molecular weigth of a protein in Daltons (Da).
# The molecular weigth  is calculated by the addition of average masses of amino acids
# in the protein and the average isotopic mass of one water molecule.
def mol_weigth(proteins):
    # We create an empty list called listrw to store the roundweigths of all the proteins in the fasta file
    listrw = []
    for i in proteins:
        # The weight of each amino acid
        weights = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
                   'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
                   'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.1,
                   'S': 87.03, 'T': 102.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06}
        weight = sum(weights[p] for p in i)
        # We add the molecular weight of one water molecule (H20) = 18.01524 Da)
        weight = (weight + 18.01524)
        # We round to two decimal places
        roundweight = round(weight, 2)
        listrw.append(roundweight)
        # Print the sum of the molecular weights
        # print("The molecular weight of this protein is:", roundweight, "Da")
    return listrw


# We are going to obtain the molecular weigth of the L34 proteins in the L34_PROT.fasta
# We create a list with the L34 proteins
L34proteins = [seq_record.seq for
               seq_record in SeqIO.parse
               (r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L34_PROT.fasta", "fasta")]


# The molecular weight of the 50S ribosomal protein L34 in the different species of the Neisseria genus
# is the same (5051.96 Da / ~ 5052 Da), which allows us to differentiate them from other genera of bacteria.
# We create the function neisseria_genus() to identify the L34 protein sequences corresponding to Neisseria genus.


def neisseria_genus(listrw):
    for i in listrw:
        if i == 5051.96:
            print("The molecular weight of this protein is:", i, "Da")
            print("Belongs to the genus Neisseria")
        else:
            print("The molecular weight of this protein is:", i, "Da")
            print("Does NOT belong to the genus Neisseria")


# We apply the function mol_weigth() to the list of proteins and safe the results in a object called L34_MW
L34_MW = mol_weigth(L34proteins)
# We apply the function neisseria_genus() to the object
neisseria_genus(L34_MW)

# We have a FASTA file (L6.fasta) with the DNA sequence of the rplF gene, which encodes the
# 50S ribosomal protein L6 in different species of Neisseria.

# Let's take a look at the original DNA file:
for seq_record in SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6.fasta",
                              "fasta"):
    # We can see the ID of the sequence and the microorganism to which it belongs
    print(seq_record.description)
    # We can see the DNA sequence
    print(repr(seq_record.seq))
    # We can see the length of the sequence
    print(len(seq_record))

# We transform the DNA sequence into RNA
file_in = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6.fasta"
file_out = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6_RNA.fasta"
with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        # We use the .transcribe() argument to obtain the RNA sequence
        seq_record.seq = seq_record.seq.transcribe()
        # write new fasta file
        r = SeqIO.write(seq_record, f_out, 'fasta')

# Let's take a look at the L6_RNA file
for seq_record in SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6_RNA.fasta",
                              "fasta"):
    # We can see the ID of the sequence and the microorganism to which it belongs
    print(seq_record.description)
    # We can see the RNA sequence
    print(repr(seq_record.seq))
    # We can see the length of the sequence
    print(len(seq_record))

# We obtain the protein sequence
file_in = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6.fasta"
file_out = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6_PROT.fasta"
with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        seq_record.seq = seq_record.seq.translate(to_stop=True)
        # write new fasta file
        r = SeqIO.write(seq_record, f_out, 'fasta')

# Let's take a look at the L6_PROT file
for seq_record in SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6_PROT.fasta",
                              "fasta"):
    print(seq_record.description)
    # We can see the PROTEIN sequence
    print(repr(seq_record.seq))
    print(len(seq_record))

# We are going to obtain the molecular weigth of the L6 proteins
L6proteins = [seq_record.seq for seq_record in
              SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\L6_PROT.fasta",
                          "fasta")]


# The molecular weight of the 50S ribosomal protein L6 is different in the species of the genus Neisseria,
# which allows us to differentiate between them.
# We create the function ident_neisseria() to identify the species of Neisseria.


def neisseria_species(listrw):
    for i in listrw:
        if i == 18879.89:
            print("The molecular weight of this protein is:", i, "Da")
            print("Neisseria flavescens")
        elif i == 18897.87:
            print("The molecular weight of this protein is:", i, "Da")
            print("Neisseria meningitidis")
        elif i == 18916.86:
            print("The molecular weight of this protein is:", i, "Da")
            print("Neisseria gonorrhoeae")
        elif i == 18907.91:
            print("The molecular weight of this protein is:", i, "Da")
            print("Neisseria mucosa")
        elif i == 18894.93:
            print("The molecular weight of this protein is:", i, "Da")
            print("Neisseria lactamica")
        elif i == 18936.95:
            print("The molecular weight of this protein is:", i, "Da")
            print("Neisseria sicca")
        elif i == 18893.9:
            print("The molecular weight of this protein is:", i, "Da")
            print("Neisseria elongata")
        else:
            print("The molecular weight of this protein is:", i, "Da")
            print("Unidentified")


# We apply the function mol_weigth() to the list of proteins and safe the results in a object called L34_MW
L6_MW = mol_weigth(L6proteins)
# We apply the function neisseria_genus() to the object
neisseria_species(L6_MW)

# We are going to perform a hierarchical cluster to see how the proteins are distributed
# The more similar two proteins are, the more similar is their sequence and therefore their molecular weight.
# In the L34 proteins dendrogram we can see that the 0,1,2 and 3 molecular weights are cluster together,
# they are the proteins that correspond to the Neisseria genus.
X = [[i] for i in L34_MW]
linked = linkage(X, 'single')
plt.figure(figsize=(10, 7))
dendrogram(linked,
           orientation="left",
           distance_sort='descending',
           show_leaf_counts=True)
plt.show()

# In the L6 protein dendrogram we can see the similarities between
# the Neisseria species.

X = [[i] for i in L6_MW]
linked = linkage(X, 'single')
plt.figure(figsize=(10, 7))
dendrogram(linked,
           orientation="left",
           distance_sort='descending',
           show_leaf_counts=True)
plt.show()

# We can apply this script to files with unknown sequences, DNA,
# RNA or proteins of both genes, and also search for other genes
# to differentiate more species based on this code.


# We transform the Unknown L6 DNA sequence (unknown_L6.fasta) into RNA
# In this file we have 4 different unknown sequences and we are going to see if if
# any of them belongs to the genus Neisseria and what specie it is.

file_in = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\unknown_L6.fasta"
file_out = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\unknown_L6_RNA.fasta"
with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        # We use the .transcribe() argument to obtain the RNA sequence
        seq_record.seq = seq_record.seq.transcribe()
        # write new fasta file
        r = SeqIO.write(seq_record, f_out, 'fasta')


# We obtain the protein sequence
file_in = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\unknown_L6.fasta"
file_out = r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\unknown_L6_PROT.fasta"
with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        seq_record.seq = seq_record.seq.translate(to_stop=True)
        # write new fasta file
        r = SeqIO.write(seq_record, f_out, 'fasta')


# We are going to obtain the molecular weigth of the L6 proteins
unknown_L6 = [seq_record.seq for seq_record in
              SeqIO.parse(r"C:\Users\Albaburillo\PycharmProjects\pythonProject2\neisseria\unknown_L6_PROT.fasta",
                          "fasta")]

unknown_L6_MW = mol_weigth(unknown_L6)
# We apply the function neisseria_genus() to the object
neisseria_species(unknown_L6_MW)


X = [[i] for i in unknown_L6_MW]
linked = linkage(X, 'single')
plt.figure(figsize=(10, 7))
dendrogram(linked,
           orientation="left",
           distance_sort='descending',
           show_leaf_counts=True)
plt.show()

# BIBLIOGRAPHY

# 1.Serrano-Coll HA, Sánchez-Jiménez M, Cardona-Castro N. Conocimiento de la microbiota de la cavidad oral a través
# de la metagenómica. Rev. CES Odont 2015; 28(2): 112-118

# 2.Morel, F., Jacquier, H., Desroches, M., Fihman, V., Kumanski, S., Cambau, E., Decousser, J. y Berçot, B. 2018,
# "Use of Andromas and Bruker MbibliographyALDI-TOF/MS MS in the identification of Neisseria", European Journal of
# Clinical Microbiology & Infectious Diseases, vol. 37, no. 12, pp. 2273-2277.

# 3.Bennett JS, Watkins ER, Jolley KA, Harrison OB, Maiden MCJ. Identifying Neisseria species by use of the 50S
# ribosomal protein L6 (rplF) gene. J Clin Microbiol. 2014;52(5):1375–81.

# 4.Malakhova, M.M., Maier, T., Kubanova, A.A., Govorun, V.M., Svistunova, T.S., Gazarian, A.O., Borovskaya, A.D.,
# Kostrzewa, M., Ilina, E.N., Vereshchagin, V.A. y Kruglov, A.N. 2009, "Direct Bacterial Profiling by Matrix-Assisted
# Laser Desorption−Ionization Time-ofFlight Mass Spectrometry for Identification of Pathogenic Neisseria",
# Journal of Molecular Diagnostics, The, vol. 11, no. 1, pp. 75-86.
