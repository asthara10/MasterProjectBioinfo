def GetData(df, cut):
	"""Retrieves the start and end data of each cut in the data frame.
	The cut variable corresponds at the number of the nucleotide that will be left before the cut."""
	import pandas
	i = 0
	cuts = []
	newcuts = []
	first = True
	
	while i < len(df.loc[:,'seqnames']):
		if first:
			chr = df.loc[i,'seqnames']
			if df.loc[i,'strand'] == "+":
				newcuts.append(df.loc[i,'start'] + cut-1)
			else:
				newcuts.append(df.loc[i,'end'] - cut-1)
			first = False
		else:
			newchr = df.loc[i,'seqnames']
			if newchr == chr:
				if df.loc[i,'strand'] == "+":
					newcuts.append(df.loc[i,'start'] + cut-1)
				else:
					newcuts.append(df.loc[i,'end'] - cut-1)
			else:
				cuts.append(newcuts)
				newcuts = []
				chr = newchr
				if df.loc[i,'strand'] == "+":
					newcuts.append(df.loc[i,'start'] + cut-1)
				else:
					newcuts.append(df.loc[i,'end'] - cut-1)
		cuts.append(newcuts)
		i+=1

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
