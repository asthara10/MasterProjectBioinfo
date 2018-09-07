def ParseFasta(file):
	"""Parses a fasta and retrieves all sequences.
	In this case sequences are returned in a list as there are no sequence names."""
	import sys
	infile = open(file)
	seqs = []
	seq = ""
	
	for line in infile:
		line = line.strip()
		if line.startswith(">"):
			if len(seq) >= 1:
				seqs.append(seq)
				seq = ""
		else:
			seq = seq + line
		top = "---" + seq[:5] + seq[-6:] + "---"
		sys.stdout.write(top)
	seqs.append(seq)
		
	return seqs
	
def addAdapters(fragments, adaptREshort, adaptRElong, adaptCas, ends5, ends3):
	"""Adds the corresponding adapters to the ends of each fragment.
	ends5 = SbfF, SbfR, ApaF, ApaR
	ends3 = SbfF, SbfR, ApaF, ApaR"""
	import sys
	constructs = []
	
	for seq in fragments:
		new = ""
		if (seq[:2] == ends5[0]) and (seq[-6:] != ends3[0]) and (seq[-6:] != ends3[1]) and (seq[-1:] != ends3[2]) and (seq[-1:] != ends3[3]):
			new = adaptRElong + seq + adaptCas9 + complement(reverse(seq))
		if (seq[:2] == ends5[1]) and (seq[-6:] != ends3[0]) and (seq[-6:] != ends3[1]) and (seq[-1:] != ends3[2]) and (seq[-1:] != ends3[3]):
			new = reverse(adaptRElong) + seq + reverse(adaptCas9) + complement(reverse(seq))
		if (seq[:5] == ends5[2]) and (seq[-6:] != ends3[0]) and (seq[-6:] != ends3[1]) and (seq[-1:] != ends3[2]) and (seq[-1:] != ends3[3]):
			new = adaptREshort + seq + adaptCas9 + complement(reverse(seq))
		if (seq[:5] == ends5[3]) and (seq[-6:] != ends3[0]) and (seq[-6:] != ends3[1]) and (seq[-1:] != ends3[2]) and (seq[-1:] != ends3[3]):
			new = reverse(adaptREshort) + seq + reverse(adaptCas9) + complement(reverse(seq))
		if (seq[-6:] == ends3[0]) and (seq[:2] != ends5[0]) and (seq[:2] != ends5[1]) and (seq[:5] != ends5[2]) and (seq[:5] != ends5[3]):
			new = complement(reverse(seq)) + adaptCas9 + seq + adaptREshort
		if (seq[-6:] == ends3[1]) and (seq[:2] != ends5[0]) and (seq[:2] != ends5[1]) and (seq[:5] != ends5[2]) and (seq[:5] != ends5[3]):
			new = complement(reverse(seq)) + reverse(adaptCas9) + seq + severse(adaptREshort)
		if (seq[-1:] == ends3[2]) and (seq[:2] != ends5[0]) and (seq[:2] != ends5[1]) and (seq[:5] != ends5[2]) and (seq[:5] != ends5[3]):
			new = complement(reverse(seq)) + adaptCas9 + seq + adaptRElong
		if (seq[-1:] == ends3[3]) and (seq[:2] != ends5[0]) and (seq[:2] != ends5[1]) and (seq[:5] != ends5[2]) and (seq[:5] != ends5[3]):
			new = complement(reverse(seq)) + reverse(adaptCas9) + seq + reverse(adaptRElong)
		if len(new) >= 1:
			sys.stdout.write("yes")
			constructs.append(new)
		
	return constructs

def Concatenate(fragments, n=10):
	"""Concatenates n times each fragment simulating the SMRT sequencing.
	Default coverage of 10."""
	import random
	rolling = []
	
	for seq in fragments:
		sign = random.randing(0, 1)
		if sign:
			rolling.append(seq*(n-random.randint(0, 2)))
		else:
			rolling.append(seq*(n+random.randint(0, 2)))
		
	return rolling
	
def addMutation(fragments, s=6, d=6, i=3):
	"""Adds mutation to the sequences.
	For long reads we assume 15% overall base mutation, being probability of substition and deletion 6% and probability of insertion 3%"""
	import random
	nucls = ["A", "C", "T", "G"]
	simulation = []
	
	for seq in fragments:
		newseq = ""
		for nt in seq:
			which = random.randint(1, 15)
			if which <= 6:
				muts = random.randint(0, 100)
				if muts <= s:
					newseq = newseq + nucls[random.randing(0, 4)]
			elif which <= 12:
				mutd = random.randint(0, 100)
				if mutd <= d:
					continue
			else:
				muti = random.randint(0, 100)
				if muti <= i:
					newseq = newseq + nt + nucls[random.randing(0, 4)]
		simulation.append(newseq)
		
	return simulation

def reverse(seq):
	"""Returns reverse of a sequence"""
	return seq[::-1]
	
def complement(seq):
	"""Returns complement of a sequence"""
	dict = {"A":"T", "C": "G", "T":"A", "G":"C"}
	seq_list = list(seq)
	seq_list = [dict[base] for base in seq_list]
	return ''.join(seq_list)
	
def saveFasta(data, name):
	outfile = open(name, "w")
	for seq in data:
		outfile.write(">\n")
		outfile.write(seq + "\n")
	
	
if __name__ == "__main__":
	import argparse
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", dest="filename", help="input fasta file", action="store")
	parser.add_argument("-as", dest="shortREadapter", help="restriction enzyme adapter sequence without sticky ends", action="store")
	parser.add_argument("-al", dest="longREadapter", help="restriction enzyme adapter sequence with sticky ends", action="store")
	parser.add_argument("-ac", dest="cas9adapter", help="Cas9 adapter sequence", action="store")
	parser.add_argument("-l", dest="fiveprime", help="possible starts in 5' end of a digested sequence", nargs='+')
	parser.add_argument("-r", dest="threeprime", help="possible ends in 3' end of a digested sequence", nargs='+')
	parser.add_argument("-o", dest="output", help="output fasta file name", action="store")
	args = parser.parse_args()
	
	saveFasta(addMutation(Concatenate(addAdapters(ParseFasta(args.filename), args.shortREadapter, args.longREadapter, args.cas9adapter, args.fiveprime, args.threeprime))), args.output)
