#*************************************************************************
#
#   Program:   get the antibody data from a file 
#   File:      get_human_heavyy.py
#
#   Version:
#   Date:      10.07.2017
#   Function:  Read in a file, step through the file one line at a time to extract data (including accession number, organism name, nucleotide sequence and DNA sequence) for gamma heavy chain antibodies. The data is then stored in a file called antibody.out.  
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


f = open("emblig-20160419-28088.xml")

for line in f:
    if "<chain>" in line:
        chain_type   = ''
        organism     = ''
        aa_sequence  = ''
        nuc_sequence = ''
        accession    = ''
    elif "<organism>" in line:
        organism = line.split(">")[1].split("<")[0]
    elif "<chain_type>" in line:
        chain_type = line.split(">")[1].split("<")[0]
    elif "<aa_sequence>" in line:
        aa_sequence = line.split(">")[1].split("<")[0]
    elif "<nuc_sequence" in line:
        nuc_sequence = line.split(">")[1].split("<")[0]
    elif "<accession>" in line:
        accession = line.split(">")[1].split("<")[0]
    elif "</chain>" in line:
        if organism == "homo sapiens" and chain_type == "heavy_gamma":
            print(">" + accession + "|protein")
            print(aa_sequence)
            print(">" + accession + "|dna")
            print(nuc_sequence)
            print()
