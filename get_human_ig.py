import re
import json
import sys

with open("antibody.dat") as f:
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


