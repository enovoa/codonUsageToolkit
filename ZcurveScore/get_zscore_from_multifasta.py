#!/usr/bin/python

from ComputeZCurveVector import *
import re
import sys

# Functions
def get_zscore_from_multivectors_file(infile, outfile,my_params):
	for line in infile.readlines():
		#Generate myseq and name_seq
		zvector=line.strip().replace("\n","")
		print zvector
		zscore=compute_z_curve_score_from_zvectors(zvector,my_params)
		outfile.write(str(zscore))
		outfile.write("\n")


def get_zscore_from_multifasta(infile, outfile,my_params):
	counter=0
	for line in infile.readlines():
		#Generate myseq and name_seq
		line=line.strip().replace("\n","")
		if re.match(">", line):
			if counter==0: #so it's the first sequence
				counter=1
				name_seq=line
				myseq=""
			else:
				#So now we are with a new sequence
				myseq=myseq.upper()   
				a = re.compile("[ATCG]*")
				m = a.match(myseq)
				if len(m.group())==len(myseq):  # so no other chars than 'ATCG'
					zscore=compute_z_curve_score(myseq,my_params)
					outfile.write(str(zscore))
					outfile.write("\n")
				
				#Continue with the new sequence
				name_seq=line
				myseq=""    
		else:
			myseq=myseq+line

	# So last sequence
	myseq=myseq.upper()   
	a = re.compile("[ATCG]*")
	m = a.match(myseq)
	if len(m.group())==len(myseq):  # so no other chars than 'ATCG'
		zscore=compute_z_curve_score(myseq,my_params)
		outfile.write(str(zscore))
		outfile.write("\n")

			

def get_zvectors_from_multifasta(infile, outfile):
	counter=0
	for line in infile.readlines():
		#Generate myseq and name_seq
		line=line.strip().replace("\n","")
		if re.match(">", line):
			if counter==0: #so it's the first sequence
				counter=1
				name_seq=line
				myseq=""
			else:
				#So now we are with a new sequence
				myseq=myseq.upper()   
				a = re.compile("[ATCG]*")
				m = a.match(myseq)
				if len(m.group())==len(myseq):  # so no other chars than 'ATCG'
					zvector=compute_z_curve_vector(myseq)
					for i in zvector:
						outfile.write(str(i))
						outfile.write("\t")
					outfile.write("\n")
				
				#Continue with the new sequence
				name_seq=line
				myseq=""    
		else:
			myseq=myseq+line
	# So last sequence
	myseq=myseq.upper()   
	a = re.compile("[ATCG]*")
 	m = a.match(myseq)
	if len(m.group())==len(myseq):  # so no other chars than 'ATCG'
		zvector=compute_z_curve_vector(myseq)
		for i in zvector:
			outfile.write(str(i))
			outfile.write("\t")
		outfile.write("\n")

