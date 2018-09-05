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
print("starting while loop")
while (i <= length(allStarts)){
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
print("end of while loop")

newGenome
