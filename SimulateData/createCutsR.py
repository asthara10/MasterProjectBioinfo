def GetData(df, cut):
	"""Retrieves the start and end data of each cut in the data frame.
	The cut variable corresponds at the number of the nucleotide that will be left before the cut."""
	import pandas
	cuts = []
	newcuts = []
	first = True
	
	for i in range(len(df.loc[:,'seqnames'])):
		if first:
			chr = df.loc[i,'seqnames']
			if (df.loc[i,'strand'] == "+") and (df.loc[i,'start'] != 0):
				newcuts.append(df.loc[i,'start'] + cut-1)
			elif (df.loc[i,'strand'] == "-") and (df.loc[i,'start'] != 0):
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
				if (df.loc[i,'strand'] == "+") and (df.loc[i,'start'] != 0):
					newcuts.append(df.loc[i,'start'] + cut-1)
				elif (df.loc[i,'strand'] == "-") and (df.loc[i,'start'] != 0):
					newcuts.append(df.loc[i,'end'] - cut-1)
	cuts.append(newcuts)

	return cuts
	
def ConvineCuts(LoL):
	"""Convines the cut places of multiple enzymes.
	Argument LoL is the List of Lists of the multiple enzymes cuts.
	Sorts from smaller to larger."""
	final = []
	longest = max(map(lambda x: len(x), LoL))
	
	for i in range(longest):
		new = []
		for L in LoL:
			try:
				new = new + L[i]
			except:
				continue
		new.append(0)
		new = set(new)
		new = list(new)
		new.sort()
		final.append(new)

	return final
