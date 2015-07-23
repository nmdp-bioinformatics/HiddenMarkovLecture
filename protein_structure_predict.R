## Install and load pacakges
if (!require("HMM", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("HMM", dependencies = TRUE)
  library("HMM", verbose=F, warn.conflicts =F)
}

## load amino acid sequences
#           <> : start
# <end> or end : end
# 20 amino acids

train_seq <- read.table("data/protein_second_train.txt", header=T, sep="\t")
test_seq <- read.table("data/protein_second_test.txt", header=T, sep="\t")

# seperate different sequence chains
seq_list <- function(seq){
  seqs <- list(data.frame(AminoAcid=character(), Structure=character()))
  seq_count <- 0
  numAA <- dim(seq)[1]
  for(id in 1:numAA){
    if(!(seq$AminoAcid[id] %in% c("<>", "end","<end>"))){
      count <- count + 1
      seqs[[seq_count]] <- rbind(seqs[[seq_count]], seq[id,])
    }else { 
      if(id == 1){
        seq_count <- seq_count + 1
      }else if(!(seq$AminoAcid[id-1] %in% c("<>", "end","<end>")) & id < numAA){
        seq_count <- seq_count + 1
        seqs[[seq_count]] <- data.frame(AminoAcid=character(), Structure=character())
      }
      count <- 0
    }
  }
  return(seqs)
}

train_seqs <- seq_list(train_seq)
test_seqs <- seq_list(test_seq)

amino_acids <- c("A","C","D","E","F","G","H","I","K","L",
                 "M","N","P","Q","R","S","T","V","W","Y")
structure <- c("_","e","h")

############### supervised learning and predictions #################
# 3 states:
# state 1- alpha-helix; (h)
# state 2- beta-sheet;  (e)
# state 3- coil;        (_)   
#

protein_seq <- data.frame(AminoAcid=character(), Structure=character())
seq_num = length(train_seqs)

# trainsitional matrix
transMat <- matrix(0, 3,3, dimnames=list(structure, structure))

# prior prob
priorMat <- matrix(0,1,3, dimnames=list("prior",structure))

# emission prob
emissMat <- matrix(0, length(amino_acids), 3, dimnames=list(amino_acids, structure))

for(id in 1:seq_num){
  ss <- as.matrix(train_seqs[[id]]["Structure"])

  for(jd in 2:length(ss)){
    transMat[ss[jd-1], ss[jd]] <- transMat[ss[jd-1], ss[jd]] + 1
  }
  tt <- as.matrix(table(train_seqs[[id]]))
  priorMat <- priorMat + colSums(tt[amino_acids, structure])
  emissMat <- emissMat + tt[amino_acids, structure]
}

## Known parameters - empirical counts 
transMat <- t(apply(transMat, 1, function(x) x/sum(x)))
emissMat <- apply(emissMat, 2, function(x) x/sum(x))
priorMat <- priorMat/sum(priorMat)

### Initialize HMM model

hmm <- initHMM(States=structure, Symbols=amino_acids, startProbs=priorMat, 
               transProbs=transMat, emissionProbs=emissMat)

print(hmm)

p <- posterior(hmm, as.matrix(test_seqs[[1]]["AminoAcid"]))

v <- viterbi(hmm, as.matrix(test_seqs[[1]]["AminoAcid"]))

logForwardProbabilities = forward(hmm,as.matrix(test_seqs[[1]]["AminoAcid"]))
print(exp(logForwardProbabilities))

#######################  unsupervised #####################
## Baum-Welch algorithm (E-M algorithm)
# 
#
hmm2 <- initHMM(States=structure, Symbols=amino_acids)
bw = baumWelch(hmm2, as.matrix(test_seqs[[1]]["AminoAcid"]), 10)

v2 <- viterbi(bw$hmm, as.matrix(test_seqs[[2]]["AminoAcid"]))
