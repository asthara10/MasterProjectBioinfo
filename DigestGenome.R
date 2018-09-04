### Loading packages

library(Biostrings)
library(stringi)
library(BSgenome.Sscrofa.UCSC.susScr3)
library(reticulate)

### Retrieving genome

genome <- BSgenome.Sscrofa.UCSC.susScr3
print("genome retrieved")

### Saving enzymes

Cas9F <- DNAString("GGTGACCCTCCTCCAGTACGNGG")
Cas9R <- DNAString("GGTCATCCACGTACTGGAGGNGG")

### Obtain cuts

CutsCas9F <- vmatchPattern(Cas9F, genome, fixed=FALSE)
CutsCas9R <- vmatchPattern(Cas9R, genome, fixed=FALSE)

### Save as data frame

RangesF = as(CutsCas9F, "data.frame")
RangesR = as(CutsCas9R, "data.frame")
print("patterns found")

### Run python script to obtain cuts

cutCas9 = 17
source_python('./createCutsR.py')
print("python imported")
listF <- GetData(RangesF, cutCas9)
listR <- GetData(RangesR, cutCas9)
print("GetData completed")
allStarts = ConvineCuts(list(listF, listR))
print("all python functions finished")

newGenome = DNAStringSet()
i <- 1
print("starting while loop")
while (i <= length(genome)){
  j <- 1
  for (nt in allStarts[[i]]){
    if (j < length(allStarts[[i]])) {
      newGenome <- append(newGenome, DNAStringSet(genome[[i]][nt:(allStarts[[i]][[j+1]]-1)]))
    } else if (j == length(allStarts[[i]])) {
      newGenome <- append(newGenome, DNAStringSet(genome[[i]][nt:length(genome[[i]])]))
    }
    j <- j+1
  }
  i <- i+1
}
print("end of while loop")

newGenome
