def findLTR(infiles, LTRfile):
	"""Returns a list with the names of all reads that contain the 5' LTR sequence at the beginning.
	infiles is a list of fasta files for each sample
	LTRfile is a fasta file containing the LTR seuqence"""
	names = []
	LTR = open(LTRfile)
	LTR.readline()
	LTRseq = LTR.readline().upper()
	for infile in infiles:
		fasta = open(infile)
		for line in fasta:
			if line.startswith(">"):
				name = line
				line = fasta.readline()
				if line[:10] in LTRseq:
					names.append(name[1:])
	return names

def writeNames(namelist, out):
	"""wirtes a text file with the names of the reads provided in the namelist argument."""
	with open(out, 'w') as outf:
		for name in namelist:
			nameR2 = name.split()
			toChange = nameR2[1].split(":")
			toChange = "2:" + toChange[1] + ":" + toChange[2] + ":" + toChange[3]
			nameR2 = nameR2[0] + " " + toChange
			outf.write("%s" %(name))
			outf.write("%s\n" %(nameR2))

def getSample(infiles):
	"""Returns a dictionary with sample names as keys and a list of reads names as values.
	infiles is the list of text files that contain the names of the reads in that sample."""
	samples = {}
	for infile in infiles:
		samples[infile] = []
		txt = open(infile)
		for line in txt:
			line = line.strip()
			samples[infile].append(line)
	return samples
			
		
def getPosition(infile, output, seqnames, samples):
	"""Outputs a file containing the sample name, the read name, the chromosome and the start and end positions.
	infile is the .sai file from the mapping
	output is the name of the output file
	seqnames are the names of the sequences containing the LTR from funcion findLTR()
	samples are the names of the samples and reads from function getSamples()"""
	seqs = seqnames
	sam = open(infile)
	outfile = open(output, "w")
	for line in sam:
		if not line.startswith("@"):
			line = line.strip()
			line = line.split("\t")
			for sample in samples:
				if (line[0] in samples[sample]) and (line[2] != "*"):
					outfile.write("%s\t%s\t%s\t%s\t%s\n" % (sample, line[0], line[2], line[3], line[7]))



if __name__ == "__main__":
	import argparse
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", dest="fastas", help="fasta file names", nargs='+')
	parser.add_argument("-t", dest="txts", help="text file names", nargs='+')
	parser.add_argument("-l", dest="ltr", help="LTR sequence fasta file names", action="store")
	parser.add_argument("-s", dest="sai", help=".sai file name containing mapping", action="store")
	parser.add_argument("-o", dest="out", help="output file name", action="store")
	args = parser.parse_args()
	
	#getPosition(args.sai, args.out, findLTR(args.fastas, args.ltr), getSample(args.txts))
	writeNames(findLTR(args.fastas, args.ltr), args.out)


