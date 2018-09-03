### Loading packages

library(Biostrings)
library(stringi)
library(BSgenome.Sscrofa.UCSC.susScr3)

### Retrieving genome

genome <- BSgenome.Sscrofa.UCSC.susScr3

### Saving enzymes

Cas9F <- DNAString("GGTGACCCTCCTCCAGTACGNGG")
Cas9R <- DNAString("GGTCATCCACGTACTGGAGGNGG")

### Obtain cuts

CutsCas9F <- vmatchPattern(Cas9F, genome)
CutsCas9R <- vmatchPattern(Cas9R, genome)

### Save as data frame

RangesF = as(CutsCas9F, "data.frame")
RangesR = as(CutsCas9R, "data.frame")

### Run python script to obtain cuts

command = "python"

path = '"./FindStartBases.py"'

cutCas9 = 17
argsF = c(RangesF, cutCas9)
argsR = c(RangesR, cutCas9)
allArgsF = c(path, argsF)
allArgsR = c(path, argsR)

startsF = system2(command, args=allArgsF, stdout=TRUE)
startsR = system2(command, args=allArgsR, stdout=TRUE)

path = '"./CombineEnzymes.py"'

startlist = c(startsF, startsR)
allArgs = c(path, startlist)

allStarts = system2(command, args=allArgs, stdout=TRUE)

newGenome = DNAStringSet()
i <- 1
for (sequence in genome){
  j <- 1
  for (nt in allStarts[i]){
    append(newGenome, DNAString(genome[[i]][nt:allStarts[i][j+1]-1]))
    j <- j+1
  }
  i <- i+1
}

newGenome
