def GetData(df, cut):
	"""Retrieves the start and end data of each cut in the data frame.
	The cut variable corresponds at the number of the nucleotide that will be left before the cut."""
	import pandas
	data = [df.loc[:,'seqnames'], df.loc[:,'start'], df.loc[:,'end'], df.loc[:,'strand']]
	i = 0
	cuts = []
	newcuts = []
	first = True
	
	while i < len(data[0]):
		if first:
			chr = data[0][i]
			if data[3][i] == "+":
				newcuts.append(data[1][i] + cut-1)
			else:
				newcuts.append(data[2][i] - cut-1)
			first = False
		else:
			newchr = data[0][i]
			if newchr == chr:
				if data[3][i] == "+":
					newcuts.append(data[1][i] + cut-1)
				else:
					newcuts.append(data[2][i] - cut-1)
			else:
				cuts.append(newcuts)
				newcuts = []
				chr = newchr
				if data[3][i] == "+":
					newcuts.append(data[1][i] + cut-1)
				else:
					newcuts.append(data[2][i] - cut-1)
		cuts.append(newcuts)

	return cuts
	
def ConvineCuts(LoL):
	"""Convines the cut places of multiple enzymes.
	Argument LoL is the List of Lists of the multiple enzymes cuts.
	Sorts from smaller to larger."""
	final = []
	i = 0
	
	while i < len(LoL[0]):
		new = []
		for L in LoL:
			new.append(L[i])
		new.append(0)
		new.sort()
		final.append(new)
		i+=1

	return final
