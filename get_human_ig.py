#*************************************************************************
#
#   Program:   get the antibody data from a file 
#   File:      get_human_ig.py
#
#   Version:
#   Date:      07.07.2017
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


import re
import json
import sys

with open("emblig-20160419-28088.xml") as f:
	lines = f.readlines()

d = {}

for i, line in enumerate(lines):
	if "<chain_type>heavy_gamma</chain_type>" in line:
		d[i] = {}
		
		j=i
		while True:
			if "<organism>" in lines[j]:
				organism = lines[j]
				organism = organism.split(">")[1].split("<")[0]
				
				d[i]["organism"] = organism
				break
			j -= 1	
	
		j=i
		while True:
			if "<accession>" in lines[j]:
				accession = lines[j]
				accession = accession.split(">")[1].split("<")[0]
				
				d[i]["accession"] = accession
				break
			j -= 1

		j=i
		while True:
			if "<aa_sequence>" in lines[j]:
				aa_sequence = lines[j]
				aa_sequence = aa_sequence.split(">")[1].split("<")[0]
				
				d[i]["aa_sequence"] = aa_sequence
				break
			j += 1

		j=i
		while True:
			if "<nuc_sequence" in lines[j]:
				nuc_sequence = lines[j]
				nuc_sequence = nuc_sequence.split(">")[1].split("<")[0]
				
				d[i]["nuc_sequence"] = nuc_sequence
				break
			j += 1

print d



with open("antibody.out", "w") as o:
	json.dump(d, o)

sys.exit()


