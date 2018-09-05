### Loading packages

library(Biostrings)
library(stringi)
library(BSgenome.Sscrofa.UCSC.susScr3)
library(reticulate)

### Retrieving genome

genome <- BSgenome.Sscrofa.UCSC.susScr3
print("genome retrieved")

### Saving enzymes

Cas9FA <- DNAString("GGTGACCCTCCTCCAGTACGAGG")
Cas9FC <- DNAString("GGTGACCCTCCTCCAGTACGCGG")
Cas9FT <- DNAString("GGTGACCCTCCTCCAGTACGTGG")
Cas9FG <- DNAString("GGTGACCCTCCTCCAGTACGGGG")
Cas9RA <- DNAString("GGTCATCCACGTACTGGAGGAGG")
Cas9RC <- DNAString("GGTCATCCACGTACTGGAGGCGG")
Cas9RT <- DNAString("GGTCATCCACGTACTGGAGGTGG")
Cas9RG <- DNAString("GGTCATCCACGTACTGGAGGGGG")

### Obtain cuts

CutsCas9FA <- vmatchPattern(Cas9FA, genome)
CutsCas9FC <- vmatchPattern(Cas9FC, genome)
CutsCas9FT <- vmatchPattern(Cas9FT, genome)
CutsCas9FG <- vmatchPattern(Cas9FG, genome)
CutsCas9RA <- vmatchPattern(Cas9RA, genome)
CutsCas9RC <- vmatchPattern(Cas9RC, genome)
CutsCas9RT <- vmatchPattern(Cas9RT, genome)
CutsCas9RG <- vmatchPattern(Cas9RG, genome)

### Save as data frame

RangesFA = as(CutsCas9FA, "data.frame")
RangesFC = as(CutsCas9FC, "data.frame")
RangesFT = as(CutsCas9FT, "data.frame")
RangesFG = as(CutsCas9FG, "data.frame")
RangesRA = as(CutsCas9RA, "data.frame")
RangesRC = as(CutsCas9RC, "data.frame")
RangesRT = as(CutsCas9RT, "data.frame")
RangesRG = as(CutsCas9RG, "data.frame")
print("patterns found")
for (chr in names(genome)) {
  if (! chr %in% RangesFA$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesFA <- rbind(RangesFA, add)
  }
}
for (chr in names(genome)) {
  if (! chr %in% RangesFC$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesFC <- rbind(RangesFC, add)
  }
}
for (chr in names(genome)) {
  if (! chr %in% RangesFT$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesFT <- rbind(RangesFT, add)
  }
}
for (chr in names(genome)) {
  if (! chr %in% RangesFG$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesFG <- rbind(RangesFG, add)
  }
}
for (chr in names(genome)) {
  if (! chr %in% RangesRA$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesRA <- rbind(RangesRA, add)
  }
}
for (chr in names(genome)) {
  if (! chr %in% RangesRC$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesRC <- rbind(RangesRC, add)
  }
}
for (chr in names(genome)) {
  if (! chr %in% RangesRT$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesRT <- rbind(RangesRT, add)
  }
}
for (chr in names(genome)) {
  if (! chr %in% RangesRG$seqnames) {
    add <- data.frame("seqnames"=chr, "start"=0, "end"=0, "width"=0, "strand"="+")
    RangesRG <- rbind(RangesRG, add)
  }
}

RangesFA = RangesFA[order(RangesFA$seqnames),]
RangesFC = RangesFC[order(RangesFC$seqnames),]
RangesFT = RangesFT[order(RangesFT$seqnames),]
RangesFG = RangesFG[order(RangesFG$seqnames),]
RangesRA = RangesRA[order(RangesRA$seqnames),]
RangesRC = RangesRC[order(RangesRC$seqnames),]
RangesRT = RangesRT[order(RangesRT$seqnames),]
RangesRG = RangesRG[order(RangesRG$seqnames),]

### Run python script to obtain cuts

cutCas9 = 17
source_python('./createCutsR.py')
print("python imported")
listFA <- GetData(RangesFA, cutCas9)
listFC <- GetData(RangesFC, cutCas9)
listFT <- GetData(RangesFT, cutCas9)
listFG <- GetData(RangesFG, cutCas9)
print("GetDataF completed")
listRA <- GetData(RangesRA, cutCas9)
listRC <- GetData(RangesRC, cutCas9)
listRT <- GetData(RangesRT, cutCas9)
listRG <- GetData(RangesRG, cutCas9)
print("GetDataR completed")
allStarts = ConvineCuts(list(listFA, listFC, listFT, listFG, listRA, listRC, listRT, listRG))
print("all python functions finished")

newGenome = DNAStringSet()
i <- 1
n <- 1
print("starting while loop")
while (i <= length(allStarts)){
  j <- 1
  for (nt in allStarts[[i]]){
    if (j < length(allStarts[[i]])) {
      newGenome <- append(newGenome, DNAStringSet(genome[[i]][nt+1:(allStarts[[i]][[j+1]]-nt)]))
	  n <- n+1
	} else if (j == length(allStarts[[i]])) {
      newGenome <- append(newGenome, DNAStringSet(genome[[i]][nt+1:(length(genome[[i]])-nt)]))
	  n <- n+1
    }
    j <- j+1
  }
  i <- i+1
}
print("end of while loop")

length(newGenome[[1]])

### Digest with restriction enzymes

SbfI <- DNAString("CCTGCAGG")
ApaLI <- DNAString("GTGCAC")
CutsSbf <- vmatchPattern(SbfI, newGenome)
CutsApa <- vmatchPattern(ApaLI, newGenome)
RangesSbf = as(CutsSbf, "data.frame")
RangesApa = as(CutsApa, "data.frame")

RangesSbf$strans=rep(c('+'), each=length(RangesSbf[[1]]))
RangesApa$strans=rep(c('+'), each=length(RangesApa[[1]]))
colnames(RangesSbf) <- c("seqnames", "group", "start", "end", "width", "strand")
colnames(RangesApa) <- c("seqnames", "group", "start", "end", "width", "strand")
RangesSbf = RangesSbf[order(RangesSbf$seqnames),]
RangesApa = RangesApa[order(RangesApa$seqnames),]

cutSbf <- 6
cutApa <- 1
listSbf <- GetData(RangesSbf, cutCas9)
listApa <- GetData(RangesApa, cutCas9)
allStarts2 = ConvineCuts(list(listSbf, listApa))

### Get final target fragments

target = DNAStringSet()
i <- 1
while (i <= length(allStarts)){
	last <- allStarts[[i]][[length(allStarts[[i]])]]
	if (length(allStarts[[i]]) > 1) {
		target <- append(target, DNAStringSet(genome[[i]][1:(allStarts[[i]][[2]])]))
		target <- append(target, DNAStringSet(genome[[i]][last+1:(length(genome[[i]])-last)]))
	} else {
		target <- append(target, DNAStringSet(genome[[i]][1:(length(genome[[i]]))]))
	}
  i <- i+1
}

length(target[[1]])
writeDNAStringSet(target, "./DigestedFragments.fasta")

