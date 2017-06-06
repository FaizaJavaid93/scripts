#*************************************************************************
#
#   Program:   calc_percentage_hydrophobic_for_codon
#   File:      calc_percentage_hydrophobics_COD.py
#
#   Version:
#   Date:      06.06.2017
#   Function:
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



#Populate a dictionary with all the amino acid codons and respective amino acids

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

#Populate a dictionary with the amino acid and their hydrophobicity value 

amino_acid_hydrophobicity = {
    'F': 2.8, 'M': 1.9, 'D':-3.5, 'L':3.8, 'A': 1.8,
    'R':-4.5, 'N':-3.5, 'I': 4.5, 'C':2.5, 'E': 3.5,
    'Q':-3.5, 'G':-0.4, 'H':-3.2, 'K':-3.9,'V': 4.2, 
    'W':-0.9, 'Y':-1.3, 'P':-1.6, 'T':-0.7,'S':-0.8
}    




def find_mutant_codons(codon):
      """ Function: Returns the mutant codons possible from a codon in as a consequence of single base changes
      
      Input: codon - a string containing a 3-letter codon
      Output: a list of strings containing the (mutant) 3-letter codons

      >>> find_mutant_codons('ATG')
      ['ATG', 'AAG', 'ATA', 'TTG', 'ATT', 'CTG', 'ACG', 'ATC', 'GTG', 'AGG']
      """

      new_codons = list()
      bases =list(codon)
      new_codons.append(codon)
      
      for this_base in ['A','T','C','G']:
            if this_base!= bases[0]:
                  new_codons.append(this_base + bases[1] + bases[2])  
            if this_base!= bases[1]:
                  new_codons.append(bases[0] + this_base + bases[2])
            if this_base!= bases[2]:
                  new_codons.append(bases[0] + bases[1] + this_base)

      return new_codons  


#check = find_mutant_codons(codon)
#print check      





def calc_percentage_hydrophobic_for_codon(codon):
      """ Function: Returns a value for the percentage of amino acids that are hydrophobic as a consequence of single base mutations in a codon 

      Input: codon - a string containing a 3-lette codon 
      Output: float - a percentage for the number of hydrophobic amino acids

      >>> calc_percentage_hydrophobic_for_codon ('ACG')
      0.20000000000000001
      """

      new_codons= find_mutant_codons(codon)

      num_hydrophobic = 0
      for codon in new_codons:
            aa= codontable[codon]
            if aa == '_':
                  continue
            else:
                  hydrophobicity= amino_acid_hydrophobicity[aa]
            if hydrophobicity >0:
                  num_hydrophobic +=1

      return(num_hydrophobic/10.0)
           



if __name__ =="__main__":
      import doctest
      doctest.testmod()

codon= 'ACG'
run_test = calc_percentage_hydrophobic_for_codon(codon)
print run_test
