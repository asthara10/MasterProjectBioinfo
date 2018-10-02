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

#adapt1 <- vmatchPattern2(adapters[[1]], reads, max.mismatch=5, with.indels=TRUE)
adapt2 <- vmatchPattern2(adapters[[2]], reads, max.mismatch=5, with.indels=TRUE) ## This one is the Cas9 adapter, the only one that is present in the CCS

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

for (i in c(1:length(adapt2))){
  stList <- start(adapt2[[i]])
  stList <- sort(stList)
  endList <- end(adapt2[[i]])
  endList <- sort(endList)
  
  if (length(stList) == 0) {
    st <- as.integer(c(1))
    end <- as.integer(c(length(reads[[i]])))
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
names <- as.character(c(1:length(adapt2)))
names(starts) <- names
names(ends) <- names

ranges <- IRangesList(start=starts, end=ends)

sequences <- list()

for (i in c(1:length(reads))) {
  if (length(ranges[[i]]) < 2) {
    sequences <- list.append(sequences, as.character(reads[[i]]))
  } else {
    newSet <- extractAt(reads[[i]], ranges[[i]])
    alignment = as.DNAbin(newSet)
    bin = mafft(alignment, opts = "--adjustdirection")
    align = ape::as.alignment(bin)
    alignForm = seqinr::as.alignment(nb=(length(align$seq)), nam=(align$nam), seq=(align$seq))
    cons = c2s(seqinr::consensus(alignForm, method = "IUPAC", type = "DNA"))
    sequences <- list.append(sequences, cons)
  }
}

write.fasta(sequences, c(1:length(sequences)), "consensusReads.fa")
