#*************************************************************************
#
#   Program:   calc_percentage_hydrophobic_for_codon
#   File:      calc_percentage_hydrophobics_SEQ.py
#  
#   Version:   
#   Date:      06.06.2017
#   Function:  Determine the percentage of amino acids that arise as a consequence of single point mutations of a single codon that is hydrophobic
#  
#   Copyright:  (c) UCL, Faiza Javaid, 2017
#   Author:     Faiza Javaid
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   Email:      faiza.javaid.12@ucl.ac.uk
#              
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If
#   someone else breaks this code, I don't want to be blamed for code
#   that does not work!
#
#   The code may not be sold commercially or included as part of a
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description: 
#   This script would be used to calculate , from a given DNA sequence, what aa arise as a result of single point mutations of a codon and what percentage of those aa are hydrophobic 
#
#*************************************************************************
#



from calc_percentage_hydrophobics_COD import find_mutant_codons
#The above function will return new codons formed as a result of mutation of a single codon
from calc_percentage_hydrophobics_COD import calc_percentage_hydrophobic_for_codon
#The above function uses the previous function to dermine what percentage of aa formed from mutating the original codon are hydrophobic

codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}



amino_acid_hydrophobicity = {
    'F': 2.8, 'M': 1.9, 'D':-3.5, 'L':3.8, 'A': 1.8,
    'R':-4.5, 'N':-3.5, 'I': 4.5, 'C':2.5, 'E': 3.5,
    'Q':-3.5, 'G':-0.4, 'H':-3.2, 'K':-3.9,'V': 4.2, 
    'W':-0.9, 'Y':-1.3, 'P':-1.6, 'T':-0.7,'S':-0.8
}    



DNA_sequence = "ATGGCGACT"

def sequence_to_codons(DNA_sequence): #Converts the DNA sequence to codons and returns as a list
      """ Function: Returns 3-leter codons from a given DNA sequence
      Input: DNA sequence- a string containing the bases A, T, G, C
      Output: a list of strings containing 3-letter codons
      >>> sequence_to_codons('ATGGCGACT')
      ['ATG', 'GCG', 'ACT']
      """

      codons =list()
      for i in range(0, len(DNA_sequence), 3):
            codon= DNA_sequence[i:i+3]
            codons.append(codon)
      return codons








if __name__ =="__main__":
      import doctest
      doctest.testmod()





codons = sequence_to_codons(DNA_sequence)
for codon in codons:
    perc_hphob = calc_percentage_hydrophobic_for_codon(codon)
    print(codon, ':', perc_hphob)
