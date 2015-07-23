# Install and load libraries

if (!require("rphast", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("rphast", dependencies = TRUE)
  library("rphast", verbose=F, warn.conflicts =F)
}

if (!require("ggplot2", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("ggplot2", dependencies = TRUE)
  library("ggplot2", verbose=F, warn.conflicts =F)
}

source("http://bioconductor.org/biocLite.R")
biocLite("seqLogo")
biocLite("ggtree")

############## Genome Notations ############
#       hg18 : Human genome                #
#    panTro2 : Chimp Genome                #
#    ponAbe2 : Pongo abelii Genome         #
#    rheMac2 : Rhesus Genome               # 
#        rn4 : Rat Genome                  #
#        mm9 : Mouse Genome                # 
#    cavPor3 : Guinea pig Genome           # 
#    canFam2 : Dog Genome                  #
#    equCab2 : Horse Genome                #
#    dasNov2 : Armadillo Genome            #
############################################

### Compare from 10 species 
seqnames <- c("hg18", "panTro2", "ponAbe2", "rheMac2", "equCab2", 
              "canFam2", "dasNov2", "mm9", "rn4", "cavPor3")

# read the phylogenetic tree model
neutralMod <- read.tm("data/placentalMammals.mod") 
neutralMod$tree <- prune.tree(neutralMod$tree, seqs=seqnames, all.but=TRUE)

## view the phylogenetic tree
require(ggtree)
tree <- read.tree("data/tree.nwk")
ggtree(tree) + theme_tree2()+ geom_text(aes(label=label), size = 5, color = "purple", hjust=-0.3)

# putative binding sites for the neuron-restrictive silencer factor (NRSF)
mafFiles <- list.files("K:/xover/Haplotype/HMM/data/NRSF", pattern="*.maf", full.names=TRUE)
nrsfNames <- sub(".maf$", "", basename(mafFiles)) # remove dir and .maf
# Read a feature file (GFF, BED, or GenePred)
nrsfSites <- read.feat("K:/xover/Haplotype/HMM/data/NRSF/NRSF.gff")

# aligned sequences
msaList <- list()
for (i in 1:length(mafFiles)) {
  smallMsa <- read.msa(mafFiles[i], seqnames=seqnames)
  smallMsa <- strip.gaps.msa(smallMsa)
  smallMsa$offset <- 0
  feat <- nrsfSites[which(nrsfSites$feature == nrsfNames[i]),]
  if (feat$strand == "-")
    smallMsa <- reverse.complement.msa(smallMsa)
  msaList[[nrsfNames[i]]] <- smallMsa
}
aggMsa <- concat.msa(msaList)
motifLen <- unique(sapply(msaList, ncol))

if (length(motifLen) != 1L) warning("all motifs should have same length!\n")

feats <- feat(seqname="hg18", src=as.character(sapply(nrsfNames, rep, motifLen)),
              feature=rep(sprintf("site.%i", 1:motifLen), length(mafFiles)),
              start=1:ncol(aggMsa), end=1:ncol(aggMsa))

mods <- phyloFit(aggMsa, init.mod=neutralMod, no.opt="ratematrix",
                 features=feats,
                 scale.only=TRUE, ninf.sites=10)

nullLike <- numeric()
lr <- numeric()
for (i in 1:motifLen) {
  nullLike[i] <- likelihood.msa(concat.msa(lapply(msaList, `[`, cols=i)), neutralMod)
  lr[i] <- mods[[i]]$likelihood - nullLike[i]
}
barplot(lr, names.arg=1:motifLen, ylab="likelihood ratio")

haveSeqLogo <- require("seqLogo")
if (haveSeqLogo) {
  m <- read.table("data/NRSF/NRSF.mtf")
  pwm <- makePWM(t(m))
  seqLogo(pwm, xfontsize=10)
} else {
  plot(c(0),c(0), type="n", xlab="", ylab="", xaxt="n", yaxt="n")
  text(0, 0, "Need to install seqLogo to see this plot")
}


## --- phylo-HMM building ---
mods[["neutral"]] <- neutralMod

# get transitional probabilities
get.trans.mat <- function(lambda, state.names, motifLen) {
  transMat <- matrix(0, nrow=length(state.names), ncol=length(state.names),
                     dimnames=list(state.names, state.names))
  transMat["neutral", "site.1"] <- lambda
  transMat["neutral", "neutral"] <- 1 - lambda
  for (i in 1:(motifLen-1))
    transMat[sprintf("site.%i", i), sprintf("site.%i", i+1)] <- 1
  transMat[sprintf("site.%i", motifLen), "neutral"] <- 1.0
  transMat
}

lambda <- 0.0001
transMat <- get.trans.mat(lambda, names(mods), motifLen)

# phylo_hmm model
nrsfHmm <- hmm(transMat)

### simulation data
simLength <- 100000
simData <- simulate.msa(mods, simLength, hmm=nrsfHmm, get.features=TRUE)

# score alignment using a general phylo-HMM
hmmScores <- score.hmm(msa=simData$msa, mod=mods, hmm=nrsfHmm, viterbi=TRUE,
                       states=sprintf("site.%i", 1:motifLen))

predicted <- hmmScores$in.states
correct <- simData$feats[substr(simData$feats$feature, 1, 4)=="site",]
numSitePredicted <- coverage.feat(predicted)
numSiteCorrect <- coverage.feat(correct)
numSiteOverlap <- coverage.feat(predicted, correct)
cat(numSitePredicted, numSiteCorrect, numSiteOverlap, "\n")

wholeRegion <- feat("hg18", src=".", feature="all",start=1, end=simLength)
sensitivity <- coverage.feat(correct, predicted)/coverage.feat(correct)
specificity <- coverage.feat(correct, predicted, wholeRegion, not=c(TRUE, TRUE, FALSE))/
  coverage.feat(correct, wholeRegion, not=c(TRUE, FALSE))
ppv <- coverage.feat(correct, predicted) /
  coverage.feat(predicted)
cat(specificity, sensitivity, ppv, "\n")

tracks <- list(as.track.feat(correct, "actual binding sites"),
               as.track.feat(predicted, "predicted binding sites", col="red"),
               as.track.wig(coord=hmmScores$post.prob.wig$coord,
                            score=hmmScores$post.prob.wig$post.prob,
                            name="hmm Posterior probilities",
                            smooth=FALSE, col="red", ylim=c(0,1)))
plot.track(tracks)

