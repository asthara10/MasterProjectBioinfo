library(Biostrings)
library(stringi)
library(reticulate)
use_python("/usr/bin/python3")

### Retrieving genome

genome <- DNAStringSet(c(c1="ACTGACTGGACTANNNNATCAGTGGCATATCATGACTTGAGGCATGACATCAGTACGATACGGACTATGC", c2="ACTACATAGCAGGTACATGACGGATACGATACGGACTTACCATACGTACACCATGACATT", c3="TTTATTTATTTATTTATTTATTTATTT"), use.names = TRUE)

### Saving enzymes

Cas9F <- DNAString("GG")
Cas9R <- DNAString("TAC")

### Obtain cuts

CutsCas9F <- vmatchPattern(Cas9F, genome)
CutsCas9R <- vmatchPattern(Cas9R, genome)

### Save as data frame

RangesF = as(CutsCas9F, "data.frame")
RangesR = as(CutsCas9R, "data.frame")
RangesF$strand=c('+', '+', '+', '+', '+', '+', '+')
RangesR$strand=c('+', '+', '+', '+', '+', '+', '+', '+', '+')
colnames(RangesF) <- c("seqnames", "group", "start", "end", "something", "strand")
colnames(RangesR) <- c("seqnames", "group", "start", "end", "something", "strand")

### Run python script to obtain cuts

source_python('./createCutsR.py')

cutCas9 = 1
FF <- GetData(RangesF, cutCas9)
RR <- GetData(RangesR, cutCas9)
allStarts = ConvineCuts(list(FF, RR))

newGenome = DNAStringSet()
i <- 1
while (i <= length(genome)){
  j <- 1
  for (nt in allStarts[[i]]){
    if (j < length(allStarts[[i]])) {
      newGenome <- append(newGenome, DNAStringSet(genome[[i]][nt+1:(allStarts[[i]][[j+1]]-nt)]))
    } else if (j == length(allStarts[[i]])) {
      newGenome <- append(newGenome, DNAStringSet(genome[[i]][nt+1:(length(genome[[i]])-nt)]))
    }
    j <- j+1
  }
  i <- i+1
}

### Digest with restriction enzymes

SbfI <- DNAString("TA")
ApaLI <- DNAString("GA")
CutsSbf <- vmatchPattern(SbfI, newGenome)
CutsApa <- vmatchPattern(ApaLI, newGenome)
RangesSbf = as(CutsSbf, "data.frame")
RangesApa = as(CutsApa, "data.frame")
RangesSbf$strand=c('+', '+', '+', '+')
RangesApa$strand=c('+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+')
colnames(RangesSbf) <- c("seqnames", "group", "start", "end", "something", "strand")
colnames(RangesApa) <- c("seqnames", "group", "start", "end", "something", "strand")
listSbf <- GetData(RangesSbf, cutCas9)
listApà <- GetData(RangesApa, cutCas9)
allStarts2 = ConvineCuts(list(listSbf, listApa))

### Get final target fragments

target = DNAStringSet()
i <- 1
while (i <= length(allStarts)){
	first <- allStarts[[i]][[1]]
	last <- allStarts[[i]][[length(allStarts[[i]])]]
	target <- append(target, DNAStringSet(genome[[i]][first+1:(allStarts[[i]][[2]]-first)]))
    target <- append(target, DNAStringSet(genome[[i]][last+1:(length(genome[[i]])-last)]))
  i <- i+1
}
