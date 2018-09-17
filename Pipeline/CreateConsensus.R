### Loading packages

library(Biostrings)

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

adapt1 <- vmatchPattern2(adapters[[1]], reads, max.mismatch=10, with.indels=TRUE)
adapt2 <- vmatchPattern2(adapters[[2]], reads, max.mismatch=10, with.indels=TRUE)

### Obtaining subreads without adapters


