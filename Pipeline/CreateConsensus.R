### Loading packages

library(Biostrings)
library(rlist)
library(seqinr)
library(ape)
library(strataG)
library(msa)

### Allowing arguments

# args = commandArgs(TRUE)
# seqResults = args[1]
# adapterSeqs = args[2]

### Saving vmatchPattern2 function

source("./vmatchPattern2.R")

### Saving fasta files to DNAStringSets

reads <- readDNAStringSet("../DataSimulation/SimulatedReads.fasta")
#seqResults
adapters <- readDNAStringSet("./Adapters.fa")
#adapterSeqs

### Finding adapter positions

adapt1 <- vmatchPattern2(adapters[[1]], reads, max.mismatch=5, with.indels=TRUE)
adapt2 <- vmatchPattern2(adapters[[2]], reads, max.mismatch=5, with.indels=TRUE)

#print("adapt1")
#for (i in c(1:length(adapt1))) {
#  print(length(adapt1[[i]]))
#}
#print("adapt2")
#for (i in c(1:length(adapt2))) {
#  print(length(adapt2[[i]]))
#}

### Defining operation functions

add <- function(a){
  return(a+1)
}
decr <- function(a){
  return(a-1)
}

### Obtaining subreads without adapters

starts <- list()
ends <- list()

for (i in c(1:length(adapt1))){
  stList <- c(start(adapt1[[i]]), start(adapt2[[i]]))
  stList <- sort(stList)
  endList <- c(end(adapt1[[i]]), end(adapt2[[i]]))
  endList <- sort(endList)
  
  if (length(stList) == 0) {
    st <- as.integer(c(0))
    end <- as.integer(c(0))
  } else {
    if ((tail(endList, n=1)+1) > length(reads[[i]])) {
      st <- as.integer(c(1, sapply(endList[1:length(endList)-1], add), tail(endList, n=1)))
    } else {
      st <- as.integer(c(1, sapply(endList, add)))
    }
    if ((head(stList, n=1)-1) < 1) {
      end <- as.integer(c(1, sapply(stList[2:length(endList)], decr), length(reads[[i]])))
    } else {
      end <- as.integer(c(sapply(stList, decr), length(reads[[i]])))
    }
  }

  starts <- list.append(starts, st)
  ends <- list.append(ends, end)
}

#starts <- sort(unlist(starts))
#ends <- sort(unlist(ends))
names <- as.character(c(1:length(adapt1)))
names(starts) <- names
names(ends) <- names

ranges <- IRangesList(start=starts, end=ends)

for (i in c(1:length(reads))) {
  newSet <- extractAt(reads[[i]], ranges[[i]])
  as.DNAbin(as.alignment(newSet))
  mafft(newSet, opts = "--adjustdirection")
}


#allFragments <- reads[ranges]
#print(allFragments)



# consensus
# https://rdrr.io/bioc/msa/man/msaConsensusSequence-methods.html
# mafft
# https://www.rdocumentation.org/packages/ips/versions/0.0-7/topics/mafft 
# https://www.rdocumentation.org/packages/strataG/versions/2.0.2/topics/mafft 
