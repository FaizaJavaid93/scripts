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

amino_acid_composition = {
    'F': 2.8, 'M': 1.9, 'D':-3.5, 'L':3.8, 'A': 1.8,
    'R':-4.5, 'N':-3.5, 'I': 4.5, 'C':2.5, 'E': 3.5,
    'Q':-3.5, 'G':-0.4, 'H':-3.2, 'K':-3.9,'V': 4.2, 
    'W':-0.9, 'Y':-1.3, 'P':-1.6, 'T':-0.7,'S':-0.8
}    


#Define function - we want to calculate the percentage of hydrophobic residues 

#def calc_percentage_hydrophobic_for_codon(codon):



#Convert the codons in the DNA to amino acid sequence
#A function to convert the sequence into codons
# break down each codon into a single base and store result in a variable 

def convert_sequence_to_codons(DNA_string):
    pro='';
    base_list = ['A','T','G','C']
    for i in range(0, len(DNA_string), 3):
        codon = DNA_string[i:i+3] # as we have not yet defined what a codon is         
        print codon
        aminoacid = codontable[codon]
        print aminoacid
        bases = list(codon)  #list to break ATG down into A, T, G
        print (bases) 
        
        
        for item in base_list:
          #iterate through the 4 possible bases
            print item
         
            for base in bases: #iterate through the codon 
                #print base
                if item ==base:
                   continue #if item is the same as the base, then ignore it, loop does not apply to it
                else:
                    for base in bases: 
                          #SOLVE THIS/// TRY STORING EACH AS A VARIABLE 
                        print bases
                        break
                
#                   for base in bases:
 #                      bases[0]= item
       #             print base
            break # breaks the for items loop - returns only one item from the base_list
        break # breaks the for loop - returns only one codon 
    
#        bases[0] = 'T'
 #       print bases
  #      new_codon = ''.join(bases)
   #     print new_codon
   #     mutant_1= codontable[new_codon]
   #     print mutant_1
          
     #   first_base = bases[0:1]  #Store each base as a variable
     #   print first_base
      #  second_base = bases[1:2]
       # print second_base
      #  third_base = bases[2:3]
     #   print third_base     
# pro+=aminoacid
            

         

    print pro #P rints the protein sequence as a string of aminoacids 
    return codon
    
sequence = "ATGGCGACTGTCGAACCGGAAACCACCCCTACTCCTAATCCCCCGACTACAGAAGAGGAGAAAACGGAATCTAATCAGGAGGTTGCTAACCCAGAACACTATATTAAACATCCCCTACAGAACAGATGGGCACTCTGGTTTTTTAAAAATGATAAAAGCAAAACTTGGCAAGCAAACCTGCGGCTGATCTCCAAGTTTGATACTGTTGAAGACTTTTGGGCTCTGTACAACCATATCCAGTTGTCTAGTAATTTAATGCCTGGCTGTGACTACTCACTTTTTAAGGATGGTATTGAGCCTATGTGGGAAGATGAGAAAAACAAACGGGGAGGACGATGGCTAATTACATTGAACAAACAGCAGAGACGAAGTGACCTCGATCGCTTTTGGCTAGAGACACTTCTGTGCCTTATTGGAGAATCTTTTGATGACTACAGTGATGATGTATGTGGCGCTGTTGTTAATGTTAGAGCTAAAGGTGATAAGATAGCAATATGGACTACTGAATGTGAAAACAGAGAAGCTGTTACACATATAGGGAGGGTATACAAGGAAAGGTTAGGACTTCCTCCAAAGATAGTGATTGGTTATCAGTCCCACGCAGACACAGCTACTAAGAGCGGCTCCACCACTAAAAATAGGTTTGTTGTTTAA"

base_list = ['A','T','G','C']

cod =convert_sequence_to_codons(sequence)





# function to calculate the percentage of hydrophobic residues

# function to mutate each base in the codon and return the resultant amino acid 
        
