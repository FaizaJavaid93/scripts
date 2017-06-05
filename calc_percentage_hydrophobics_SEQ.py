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



DNA_sequence = "ATGGCGACT"

def sequence_to_codons(DNA_sequence): #Converts the DNA sequence to codons and returns as a list
      """ Returns the number of hydrophobic amino acids from new codons
      >>> sequence_to_codons ('ATGGCGACT')
      1.5
      """

      codons =list()
      for i in range(0, len(DNA_sequence), 3):
            codon= DNA_sequence[i:i+3]
            codons.append(codon)
      return codons



def print_codons(DNA_sequence): #Returns the codons on a separate line 
      codons = sequence_to_codons(DNA_sequence)
      for codon in codons:
            return codon 



def find_mutant_codons(codon):
      new_codons = list()
      codons = sequence_to_codons(codon)

      for codon in codons:

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






def calc_percentage_hydrophobic_for_codon(codon):
      new_codons= find_mutant_codons(codon)

      num_hydrophobic = 0
      for codon in new_codons:
            print codon
            aa= codontable[codon]
            print aa
            if aa == '_':
                  continue
            else:
                  hydrophobicity= amino_acid_hydrophobicity[aa]
            print hydrophobicity
            
            if hydrophobicity >0:
                  num_hydrophobic +=1

      return(num_hydrophobic/10.0)



if __name__ =="__main__":
      import doctest
      doctest.testmod()



run_test = calc_percentage_hydrophobic_for_codon(DNA_sequence)
print run_test


