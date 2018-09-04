### Loading packages

library(Biostrings)
library(stringi)
library(BSgenome.Sscrofa.UCSC.susScr3)
library(reticulate)

### Retrieving genome

genome <- BSgenome.Sscrofa.UCSC.susScr3

### Saving enzymes

Cas9F <- DNAString("GGTGACCCTCCTCCAGTACGNGG")
Cas9R <- DNAString("GGTCATCCACGTACTGGAGGNGG")

### Obtain cuts

CutsCas9F <- vmatchPattern(Cas9F, genome, fixed=FALSE)
CutsCas9R <- vmatchPattern(Cas9R, genome, fixed=FALSE)

### Save as data frame

RangesF = as(CutsCas9F, "data.frame")
RangesR = as(CutsCas9R, "data.frame")

### Run python script to obtain cuts

cutCas9 = 17
source_python('./createCutsR.py')
allStarts = ConvineCuts(list(GetData(RangesF, cutCas9), GetData(RangesR, cutCas9)))

newGenome = DNAStringSet()
i <- 1
while (i <= length(genome)){
  j <- 1
  for (nt in allStarts[[i]]){
    append(newGenome, DNAString(genome[[i]][nt:allStarts[[i]][[j+1]]-1]))
    j <- j+1
  }
  i <- i+1
}

newGenome
