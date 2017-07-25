#*************************************************************************
#
#   Program:   calc_percentage_hydrophobic_for_codon
#   File:      calc_percentage_hydrophobics_COD.py
#
#   Version:
#   Date:      25.07.2017
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



#Populate a dictionary with all the amino acid codons and respective amino acids. This will help to assign the codons and determine what aa they code for 


import re

codontable = {
    'ata':'I', 'atc':'I', 'att':'I', 'atg':'M',
    'aca':'T', 'acc':'T', 'acg':'T', 'act':'T',
    'aac':'N', 'aat':'N', 'aaa':'K', 'aag':'K',
    'agc':'S', 'agt':'S', 'aga':'R', 'agg':'R',
    'cta':'L', 'ctc':'L', 'ctg':'L', 'ctt':'L',
    'cca':'P', 'ccc':'P', 'ccg':'P', 'cct':'P',
    'cac':'H', 'cat':'H', 'caa':'Q', 'cag':'Q',
    'cga':'R', 'cgc':'R', 'cgg':'R', 'cgt':'R',
    'gta':'V', 'gtc':'V', 'gtg':'V', 'gtt':'V',
    'gca':'A', 'gcc':'A', 'gcg':'A', 'gct':'A',
    'gac':'D', 'gat':'D', 'gaa':'E', 'gag':'E',
    'gga':'G', 'ggc':'G', 'ggg':'G', 'ggt':'G',
    'tca':'S', 'tcc':'S', 'tcg':'S', 'tct':'S',
    'ttc':'F', 'ttt':'F', 'tta':'L', 'ttg':'L',
    'tac':'Y', 'tat':'Y', 'taa':'_', 'tag':'_',
    'tgc':'C', 'tgt':'C', 'tga':'_', 'tgg':'W'
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

      >>> find_mutant_codons('atg')
      ['atg', 'aag', 'ata', 'ttg', 'att', 'ctg', 'acg', 'atc', 'gtg', 'agg']
      """

      new_codons = list()
#Convert each codon (a string of 3 letters) into a list so you can iterate through and replace each position in the codon with a new base (with access via the index)
      bases =list(codon)
      new_codons.append(codon)

#The new base will either be at position 1, 2 or 3 provided that the new base (this_base) is not the same as the base it is trying to replace

      for this_base in ['a','t','c','g']:
            if this_base!= bases[0]:
                  new_codons.append(this_base + bases[1] + bases[2])  
            if this_base!= bases[1]:
                  new_codons.append(bases[0] + this_base + bases[2])
            if this_base!= bases[2]:
                  new_codons.append(bases[0] + bases[1] + this_base)
#Return a list called new_codons populated with the new codons

      return new_codons  







def calc_percentage_hydrophobic_for_codon(codon):
      """ Function: Returns a value for the percentage of amino acids that are hydrophobic as a consequence of single base mutations in a codon 

      Input: codon - a string containing a 3-lette codon 
      Output: float - a percentage for the number of hydrophobic amino acids

      >>> calc_percentage_hydrophobic_for_codon ('gcg')
      0.20000000000000001
      """
#Get mutant codons 
#Determine the amino acid and hydrophobicity correspoding to the codon
      new_codons= find_mutant_codons(codon)
      if '-'not in codon and 'n' not in codon and 'r' not in codon:
          
            aa= codontable[codon]
            hphob= amino_acid_hydrophobicity[aa]
     
#Create a counter

            num_hydrophobic = 0
            for codon in new_codons:
                  if '-' not in codon and 'n' not in codon and 'r' not in codon: 
                        aa= codontable[codon]
                        if aa == '_':
                              continue
                        else:
                              hydrophobicity= amino_acid_hydrophobicity[aa]
#Determine if new aa is MORE  hydrophobic than current - values over the hphob value of original codon i.e. positive values indicate that aa is hydrophobic 

                              if hydrophobicity >hphob:
                                    num_hydrophobic +=1
#Get percentage - use a float y

            return    (num_hydrophobic/10.0)
           






if __name__ =="__main__":
      import doctest
      doctest.testmod()




f= open("mutations.txt")

total_codon= 0
exp_hydrophobic = 0
not_hydrophobic= 0
obs_hydrophobic = 0
for line in f: 


#Find out if the new codon which formed by mutation is more or less hydrophobic than the original codon in the gemline sequence. Make a count of the number of times you get a more vs a  less hydrophobic 

#From the file mutations.txt, pull out the lines of interest, and from that extract the information of interest. i.e. germline codon that has been mutated and new resultant amino acid 

      if re.match (r'\s+\d+', line):
            codon=  line.split(":")[1].split(" ")[1]
            if '-'not in codon and 'n' not in codon and 'r' not in codon:


#Print the hydrophobicity of the original germline amino acid 
# Keep a running count of the number of codons analysed 
                  total_codon +=1
                  aa= codontable[codon]
                  old_hphob= amino_acid_hydrophobicity[aa]
                  print old_hphob
                  perc_hydrophobic = calc_percentage_hydrophobic_for_codon(codon)
                  print perc_hydrophobic
                  exp_hydrophobic = perc_hydrophobic + exp_hydrophobic

                  if re.match (r'\s+\d+', line):
                        new_aa=  line.split(":")[1].split(" ")[4]
                        new_aa= new_aa.strip('\n')


#print the hydrophobicity of the new amino acid
                        if '?' not in new_aa:
                              new_hphob = amino_acid_hydrophobicity[new_aa]
                              print new_hphob
                              if new_hphob > old_hphob: 
                                    obs_hydrophobic +=1
                                    
                              else: 
                                    not_hydrophobic +=1
                                    


print ("The total number of codons is : %s" %(total_codon))
print("The number of non-hydrophobic amino acids is: %s" %(not_hydrophobic))
print ("The number of observed hydrophobic amino acids is: %s"%(obs_hydrophobic))
print ("The number of expected hydrophobic amino acids is: %s"%(exp_hydrophobic))

            
