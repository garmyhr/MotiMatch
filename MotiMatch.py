import sys
import numpy as np
import pandas as pd
from math import log2

# Global lists representing values in the first column
# Used for printing
NUCLEOTIDES = ['A', 'C', 'G', 'T']
NUCLEOTIDES_SUM = ['A', 'C', 'G', 'T', 'Sum']

# Main
def main():

    M = sys.argv[1]          # Name of the motif-file
    G = sys.argv[2]          # Name of the genome-file
    N = int(sys.argv[3])     # Number of hits wanted
    P = float(sys.argv[4])   # Total pseudocounts

    # a)
    # Loading sequences into memory
    genome_seq = get_single_sequence(G)
    motif_sequences = get_multiple_sequences(M)

    # Used for printing
    X = len(motif_sequences[0])
    COL = [i for i in range(1, X+1)]

    # Calculations
    background_freq = count_rel_freq(genome_seq)
    PFM = calculate_PFM(motif_sequences); consensus = find_consensus(PFM)
    PPM, PWM, IC = calculate_matrices(PFM, P, background_freq)

    # b)
    print("Genome length: %d" % len(genome_seq))

    # c)
    print("Number of sequences: %d" % len(motif_sequences))
    print("Motif lenght: %d" % X)

    # f)
    print("Motif consensus: %s" % consensus)

    # d)
    print("\nBackground frequencies:"); print(background_freq)

    # e)
    print("\nPositional Frequency Matrix"); print(pd.DataFrame(data=PFM, index=NUCLEOTIDES, columns=COL))

    # h)
    print("\nPositional Probability Matrix"); print(pd.DataFrame(data=PPM, index=NUCLEOTIDES, columns=COL))

    # i)
    print("\nPositional Weight Matrix"); print(pd.DataFrame(data=matrix_with_sum(PWM), index=NUCLEOTIDES_SUM, columns=COL))

    # j)
    IC = matrix_with_sum(IC)
    print("\nInformation Content"); print(pd.DataFrame(data=IC, index=NUCLEOTIDES_SUM, columns=COL))

    # k)
    total_information = sum(IC[4,:])
    print("\nTotal Information Content: %f" % total_information)

    # l)
    # m)
    top_motifs = scan_genome(genome_seq, PWM, X, N)

    # n)
    print("\nTop %d sequences found:\n" % N)
    print_results(top_motifs, N)

    # Used for testing the input motifs
    # check = []
    # for motif in motif_sequences:
    #     score = score_sequence(PWM, motif, X)
    #     add = Motif(score, motif, 0)
    #     check.append(add)
    # check.sort(key=lambda x: x.score, reverse=True)
    # check[0].print_motif()

    sys.exit()

# Calculates the Positional Frequency Matrix of a set of sequences
def calculate_PFM(motif_sequences):

    # X = Length of the motif
    # Y = Number of motifs found in motif_sequences
    X = len(motif_sequences[0])
    Y = len(motif_sequences)

    PFM = np.zeros((4, X))

    # Building the sequences with letter at position i from all motifs
    positional_sequences = []
    for i in range(0, X):
        letters_in_position = []
        for motif in motif_sequences:
            letters_in_position.append(motif[i])
        positional_sequences.append(letters_in_position)

    # Counting instances of each nucleotide in the positional sequences
    # and adding them to PFM
    for i in range(X):
        A = positional_sequences[i].count('A')
        C = positional_sequences[i].count('C')
        G = positional_sequences[i].count('G')
        T = positional_sequences[i].count('T')
        PFM[0, i] = A
        PFM[1, i] = C
        PFM[2, i] = G
        PFM[3, i] = T

    return PFM

# Finds the consensus-sequence based on the Positional Frequency Matrix
def find_consensus(PFM):

    m, n = np.shape(PFM)

    # Find consensus-string using the PFM
    consensus = ""
    for i in range(0, n):
        counts = PFM[0:, i]
        largest = max(counts)

        for j in range(len(counts)):
            if counts[j] == largest:
                consensus += NUCLEOTIDES[j]

    return consensus

# Parameters:
# PFM - Postitional Frequency Matrix
# P - value of pseudocounts
# background_freq - values for the background frequencies
#
# Returns:
# PPM - Positional Probability Matrix
# PWM - Positional Weight Matrix
# IC - Information Content
def calculate_matrices(PFM, P, background_freq):

    # Creating multiple matrices to return
    m, n = np.shape(PFM)
    PPM = np.zeros((m,n))
    PWM = np.zeros((m,n))
    IC = np.zeros((m,n))

    # Add pseudocounts
    # Compute relative frequencies of each nucleotide in each position, creating the Positional Probability Matrix
    for i in range(0, n):
        values = PFM[0:, i]
        total = sum(values)
        for j in range(0, m):
            p_i = background_freq[NUCLEOTIDES[j]]
            PPM[j, i] = (PFM[j, i] + (p_i * P)) / (total + P)


    # Compute log-odds score for each nucleotide, creating the Positional Weight Matrix
    for i in range(0, n):
        for j in range(0, m):
            p_i = background_freq[NUCLEOTIDES[j]]
            PWM[j, i] = log2(PPM[j, i] / p_i)

    # Compute informational content of each position in the motif
    for i in range(0, n):
        for j in range(0, m):
            IC[j, i] = PPM[j, i] * PWM[j, i]

    return PPM, PWM, IC

# Parameters:
# genome - the entire sequence of the genome to be scanned
# PWM - Positional Weight Matrix
# X - Length of the motif
# N - int representing how many matches to return
#
# Returns:
# top_sequences - a list containing the N highest scoring motifs
def scan_genome(genome, PWM, X, N):

    L = len(genome)
    top_sequences = []

    # For every position in genome
    for i in range(L-X):

        sequence = genome[i : X+i]
        score = score_sequence(PWM, sequence, X)
        motif = Motif(score, sequence, i)

        if len(top_sequences) < N:
            top_sequences.append(motif)
        else:
            for j in range(N):
                if (top_sequences[j].score < motif.score):
                    top_sequences[j] = motif; break;

    top_sequences.sort(key=lambda x: x.score, reverse=True)
    return top_sequences

# Calculates the score of a given sequence, based on a PWM
def score_sequence(PWM, sequence, X):

    values = []
    for i in range(X):
        if sequence[i] == 'A':
            values.append(PWM[0, i])
        elif sequence[i] == 'C':
            values.append(PWM[1, i])
        elif sequence[i] == 'G':
            values.append(PWM[2, i])
        elif sequence[i] == 'T':
            values.append(PWM[3, i])

    return sum(values)

# Counts the frequency of all instances in a sequence
def count_rel_freq(sequence):
    freq = {}
    for item in sequence:
        if(item in freq):
            freq[item] += 1.0
        else:
            freq[item] = 1.0

    for key in freq:
        freq[key] = freq[key]/len(sequence)

    return freq

# Loads a single sequence from FASTA-file into memory as an array
def get_single_sequence(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    sequence = ""
    for i in range(1, len(lines)):
        lines[i] = lines[i].strip('\n')
        sequence += lines[i]

    return list(sequence)

# Loads multiple sequences from FASTA-file into memory as a 2D-array
def get_multiple_sequences(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    sequences = []

    for i in range(1, len(lines), 2):
        cur = lines[i].strip('\n')
        sequences.append(list(cur))

    return sequences

# Takes a matrix and adds a bottom row with the sum of each column
def matrix_with_sum(matrix):
    m, n = np.shape(matrix)
    RetMat = np.zeros((m+1,n))

    for i in range(0, n):
        values = matrix[0:4, i]
        total = sum(values)
        for j in range(0, m+1):

            if j == 4:
                RetMat[j, i] = total
            else:
                RetMat[j, i] = matrix[j, i]

    return RetMat

# Prints all the result to terminal
def print_results(motifs, N):

    arr = np.zeros((N, 3), dtype='object')
    for i in range(N):
        motif = motifs[i]
        arr[i, 0] = motif.position
        arr[i, 1] = motif.score
        arr[i, 2] = ''.join(motif.sequence)

    print(pd.DataFrame(data=arr, index=range(1,N+1), columns=['Position', 'Score', 'Sequence']))

# A Motif class, used to store additional information
class Motif:
    def __init__(self, score, sequence, position):
        self.score = score
        self.sequence = sequence
        self.position = position

if __name__ == "__main__":
    main()
