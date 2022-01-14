
library(dada2)
library(ggplot2)
library(phyloseq)
library(Biostrings)

home_dir <- "/users/tg/jwillis/SLL"



### SEE TUTORIAL HERE:   https://benjjneb.github.io/dada2/tutorial.html

# ****************************************************** #

pathMock <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/mock_communities"
length(list.files(pathMock))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs.mocks <- sort(list.files(pathMock, pattern="_R1_001.fastq", full.names = TRUE))
fnRs.mocks <- sort(list.files(pathMock, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: datax_xxx_SAMPLENAME_XXX.fastq
sample.names.mocks <- sapply(strsplit(basename(fnFs.mocks), "_"), `[`, 3)


# ****************************************************** #
#### EXAMINE QUALITY PROFILES OF FORWARD AND REVERSE READS #### 
# ****************************************************** #
plotQualityProfile(fnFs.mocks[1:4]) + 
  scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5))
# From this it seems should truncate forward reads at least at position 250 (trim last 50 nucs), and remove first 25

plotQualityProfile(fnRs.mocks[1:4]) +
  scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5))
# reverse have worse quality at ends, truncate at position 200, and remove first 10



# ****************************************************** #
#### PERFORM FILTERING AND TRIMMING #### 
# ****************************************************** #
filt_path.mocks <- file.path(pathMock, "filtered_mocks") # Place filtered files in filtered_reads/ subdirectory
filtFs.mocks <- file.path(filt_path.mocks, paste0(sample.names.mocks, "_F_filt.fastq.gz"))
filtRs.mocks <- file.path(filt_path.mocks, paste0(sample.names.mocks, "_R_filt.fastq.gz"))


# use maxEE=2 cause # of expected errors is better than average quality scores; http://www.drive5.com/usearch/manual/expected_errors.html
# out.mocks <- filterAndTrim(fnFs.mocks, filtFs.mocks, fnRs.mocks, filtRs.mocks, 
#                      truncLen=c(250,200), maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
#                      trimLeft=c(25,10), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
# 
# mean(out.mocks[,2] / out.mocks[,1])


out.mocks <- filterAndTrim(fnFs.mocks, filtFs.mocks, fnRs.mocks, filtRs.mocks, 
                           truncLen=c(275,225), maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                           trimLeft=c(25,10), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

mean(out.mocks[,2] / out.mocks[,1]) # 0.6093571






# ****************************************************** #
#### MAKE LIST OF LISTS FOR EACH MERGE METHOD ####
# ****************************************************** #

# outputs for these different merging methods were produced using the script
#    /users/tg/jwillis/SLL/code/joiners_mocks.sh
# except for "dada2" which will use the mergePairs() function from dada2 package and therefore is not merged until later
merge_methods <- c("dada2","dada2_5","dada2_10","dada2_20","PEAR","pandaseq_sb","pandaseq_ea","pandaseq_flash",
                   "pandaseq_pear","pandaseq_rdp_mle","pandaseq_stitch","pandaseq_uparse","fqj")
all_mergers.mocks <- list()

for (met in merge_methods) all_mergers.mocks[[met]] <- list()


# ****************************************************** #
#### LEARN THE ERROR RATES #### 
# ****************************************************** #


for (met in merge_methods) {
  print(met)
  if (met == "dada2") {
    all_mergers.mocks[[met]][["errF"]] <- learnErrors(filtFs.mocks, multithread=TRUE, nreads = 2e+6)
    all_mergers.mocks[[met]][["errR"]] <- learnErrors(filtRs.mocks, multithread=TRUE, nreads = 2e+6)
  
    } else {
      pathMock.joined <- sprintf("%s/joined/%s", filt_path.mocks, met)
      mocks.joined <- sort(list.files(pathMock.joined, pattern=".fastq.gz", full.names = TRUE))
      
      all_mergers.mocks[[met]][["err"]] <- learnErrors(mocks.joined, multithread=TRUE, nreads = 2e+6)
  }
}

# show error rates for each possible transition (eg. A->C, A->G, …) 
plotErrors(all_mergers.mocks[["dada2"]][["errF"]], nominalQ=TRUE) # looks all good
plotErrors(all_mergers.mocks[["dada2"]][["errR"]], nominalQ=TRUE) # looks all good
plotErrors(all_mergers.mocks[["PEAR"]][["err"]], nominalQ=TRUE) # looks all good
plotErrors(all_mergers.mocks[["pandaseq_sb"]][["err"]], nominalQ=TRUE) # looks mostly fine
plotErrors(all_mergers.mocks[["pandaseq_ea"]][["err"]], nominalQ=TRUE) # looks not excellent
plotErrors(all_mergers.mocks[["pandaseq_flash"]][["err"]], nominalQ=TRUE) # looks bad
plotErrors(all_mergers.mocks[["pandaseq_pear"]][["err"]], nominalQ=TRUE) # looks alright
plotErrors(all_mergers.mocks[["pandaseq_rdp_mle"]][["err"]], nominalQ=TRUE) # looks all good
plotErrors(all_mergers.mocks[["pandaseq_stitch"]][["err"]], nominalQ=TRUE) # looks mostly fine
plotErrors(all_mergers.mocks[["pandaseq_uparse"]][["err"]], nominalQ=TRUE) # looks alright
plotErrors(all_mergers.mocks[["fqj"]][["err"]], nominalQ=TRUE) # looks all good




# from this post about not converging after 10 rounds or error estimation https://github.com/benjjneb/dada2/issues/77
dada2:::checkConvergence(all_mergers.mocks[["dada2"]][["errF"]])
dada2:::checkConvergence(all_mergers.mocks[["dada2"]][["errR"]])
dada2:::checkConvergence(all_mergers.mocks[["PEAR"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["pandaseq_sb"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["pandaseq_ea"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["pandaseq_flash"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["pandaseq_pear"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["pandaseq_rdp_mle"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["pandaseq_stitch"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["pandaseq_uparse"]][["err"]])
dada2:::checkConvergence(all_mergers.mocks[["fqj"]][["err"]])







# ****************************************************** #
#### DEREPLICATION #### 
# ****************************************************** #
# Combines all identical sequences into unique sequences with abundances
# Reduces computation time by eliminating redundant comparisons
# *** UNIQUE TO DADA2: retains summary of quality info for each unique sequence
#     This is an average of positional qualities for duplicated reads to inform the error model

for (met in merge_methods) {
  print(met)
  if (met == "dada2") {
    all_mergers.mocks[[met]][["derepFs"]] <- derepFastq(filtFs.mocks, verbose=TRUE)
    all_mergers.mocks[[met]][["derepRs"]] <- derepFastq(filtRs.mocks, verbose=TRUE)
    # Name the derep-class objects by the sample names
    names(all_mergers.mocks[[met]][["derepFs"]]) <- sample.names.mocks
    names(all_mergers.mocks[[met]][["derepRs"]]) <- sample.names.mocks
    
  } else {
    pathMock.joined <- sprintf("%s/joined/%s", filt_path.mocks, met)
    mocks.joined <- sort(list.files(pathMock.joined, pattern=".fastq.gz", full.names = TRUE))
    
    all_mergers.mocks[[met]][["derep"]] <- derepFastq(mocks.joined, verbose=TRUE)
    # Name the derep-class objects by the sample names
    names(all_mergers.mocks[[met]][["derep"]]) <- sample.names.mocks
  }
}






# ****************************************************** #
#### SAMPLE INFERENCE #### 
# ****************************************************** #
# Apply core sequence-variant inference algorithm to the dereplicated data
for (met in merge_methods) {
  print(met)
  if (met == "dada2") {
    all_mergers.mocks[[met]][["dadaFs"]] <- dada(all_mergers.mocks[[met]][["derepFs"]], err=all_mergers.mocks[[met]][["errF"]], multithread=TRUE)
    all_mergers.mocks[[met]][["dadaRs"]] <- dada(all_mergers.mocks[[met]][["derepRs"]], err=all_mergers.mocks[[met]][["errR"]], multithread=TRUE)
    
  } else {
    all_mergers.mocks[[met]][["dada"]] <- dada(all_mergers.mocks[[met]][["derep"]], err=all_mergers.mocks[[met]][["err"]], multithread=TRUE)
  }
}





# ****************************************************** #
#### MERGE PAIRED READS #### 
# ****************************************************** #
# Further reduce spurious sequence variants by merging overlapping reads
# F and R reads must be in matching order at time of dereplication

# for the dada2 merging method, will need to merge now, and will give it the same name as for the other merging methods,
# which do not need to perform this step since they began as merged reads
all_mergers.mocks[["dada2"]][["dada"]] <- mergePairs(all_mergers.mocks[["dada2"]][["dadaFs"]], all_mergers.mocks[["dada2"]][["derepFs"]], 
                                                     all_mergers.mocks[["dada2"]][["dadaRs"]], all_mergers.mocks[["dada2"]][["derepRs"]], 
                                                     verbose=TRUE)

# add extra names for different maxMismatch values in the mergePairs() function
merge_methods.extras <- c(merge_methods, "dada2_5","dada2_10","dada2_20")

all_mergers.mocks[["dada2_5"]][["dada"]] <- mergePairs(all_mergers.mocks[["dada2"]][["dadaFs"]], all_mergers.mocks[["dada2"]][["derepFs"]], 
                                                       all_mergers.mocks[["dada2"]][["dadaRs"]], all_mergers.mocks[["dada2"]][["derepRs"]], 
                                                       verbose=TRUE, maxMismatch = 5)

all_mergers.mocks[["dada2_10"]][["dada"]] <- mergePairs(all_mergers.mocks[["dada2"]][["dadaFs"]], all_mergers.mocks[["dada2"]][["derepFs"]], 
                                                        all_mergers.mocks[["dada2"]][["dadaRs"]], all_mergers.mocks[["dada2"]][["derepRs"]], 
                                                        verbose=TRUE, maxMismatch = 10)

all_mergers.mocks[["dada2_20"]][["dada"]] <- mergePairs(all_mergers.mocks[["dada2"]][["dadaFs"]], all_mergers.mocks[["dada2"]][["derepFs"]], 
                                                        all_mergers.mocks[["dada2"]][["dadaRs"]], all_mergers.mocks[["dada2"]][["derepRs"]], 
                                                        verbose=TRUE, maxMismatch = 20)






# ****************************************************** #
#### CONSTRUCT SEQUENCE TABLE #### 
# ****************************************************** #
# Gives a table of sequences from samples -- higher res OTU table than traditional methods
for (met in merge_methods.extras) {
  print(met)
  # Rows are sample names, cols are sequence variants
  all_mergers.mocks[[met]][["seqtab"]] <- makeSequenceTable(all_mergers.mocks[[met]][["dada"]])
  # Inspect distribution of sequence lengths
  print(table(nchar(getSequences(all_mergers.mocks[[met]][["seqtab"]]))))
}

# [1] "dada2"
# 250 261 346 353 366 386 404 405 406 407 408 409 410 411 412 417 418 421 422 423 424 425 426 429 430 431 432 433 441 
#   2   1   1   2   1   1   6  53  10   3   3   3   5   4   2   7   1   1   2   5   6  28   1  19 196  51   2   1   1 
# 
# [1] "PEAR"
# 149 154 185 187 221 223 230 264 301 313 318 334 376 404 405 406 407 408 409 410 411 412 413 416 417 418 421 422 423 424 425 426 429 430 431 432 433 451 
#   1   1   1   1   1   1   1   1   1   1   1   1   1  12  71  10   2   5   3   5   7   2   1   1   7   2   1   5   6   6  25   2  32 215  85   2   2   1 
# 
# [1] "pandaseq_sb"
# 264 301 313 318 334 376 404 405 406 407 408 409 410 411 412 413 416 417 418 421 422 423 424 425 426 429 430 431 432 433 451 463 
#   1   1   1   1   1   1  12  70  10   2   5   3   5   7   2   1   1   7   2   1   5   6   6  25   3  32 222  80   2   2   1   1 
# 
# [1] "pandaseq_ea"
# 264  273  301  404  405  406  408  409  410  411  412  413  416  417  418  421  422  423  424  425  426  429  430  431  432 
#   1    1    1    9   71   99    5    3   10    5    1    1    1   36    1    1    4   19    3   91    2   22 1190  107    2 
# 
# [1] "pandaseq_flash"
# 461 462 463 
#   6   3  26 
# 
# [1] "pandaseq_pear"
# 264 301 313 318 334 376 404 405 406 407 408 409 410 411 412 413 416 417 418 421 422 423 424 425 426 429 430 431 432 433 462 463 
#   1   1   1   1   1   1  12 199  71   2   5   3   7   7   2   1   1  21   2   1   5  14   6  36   3  33 542 105   2   2   1  10 
# 
# [1] "pandaseq_rdp_mle"
# 264 301 313 318 334 376 404 405 406 407 408 409 410 411 412 413 416 417 418 421 422 423 424 425 426 429 430 431 432 433 461 462 463 
#   1   1   1   1   1   1  12  68   9   2   5   3   5   7   2   1   1   7   2   1   5   6   6  24   2  30 211  80   2   2   1   2   7 
# 
# [1] "pandaseq_stitch"
# 264 301 404 405 406 408 409 410 411 412 413 416 417 418 421 422 423 424 425 426 429 430 431 432 462 463 
#   1   1  11  46  10   5   3   3   5   1   1   1   6   1   1   4   5   3   7   3  23 163  55   2   6   6 
# 
# [1] "pandaseq_uparse"
# 264 301 313 318 334 404 405 406 408 409 410 411 412 413 416 417 418 421 422 423 424 425 426 429 430 431 432 
#   1   1   1   1   1  10  44   9   5   3   3   6   2   1   1   5   1   1   4   5   4   7   2  22 148  50   2 
# 
# [1] "fqj"
# 264 301 313 318 334 376 404 405 406 407 408 409 410 411 412 413 416 417 418 421 422 423 424 425 426 429 430 431 432 433 
#   1   1   1   1   1   1  11  67  10   2   5   3   5   6   2   1   1   8   2   1   5   7   6  24   1  30 211  75   2   2 
#
# [1] "dada2_5"
# 250 261 280 346 353 366 386 404 405 406 407 408 409 410 411 412 417 418 421 422 423 424 425 426 429 430 431 432 433 441 
#   2   1   1   1   2   1   1   7  66  10   3   4   3   5   4   2   9   1   1   3   7   8  37   4  24 352 133   8   1   1 
# 
# [1] "dada2_10"
# 250 261 280 346 353 366 386 404 405 406 407 408 409 410 411 412 417 418 421 422 423 424 425 426 427 429 430 431 432 433 441 
#   2   1   1   1   2   1   1   7  72  12   4   4   4   7   4   4  13   2   1   3  13  12  55   7   3  53 543 300  24   3   1 
# 
# [1] "dada2_20"
# 250 261 280 346 353 366 386 404 405 406 407 408 409 410 411 412 413 414 416 417 418 419 420 421 422 423 424 425 426 427 429 430 431 432 433 441 
#   2   1   1   1   2   1   1   9 140  52  14   7  10  12  18  16   3   1   2  21   8   4   7   6   5  27  20  93  26   9  54 576 324  31   5   1



# ****************************************************** #
#### REMOVE CHIMERAS #### 
# ****************************************************** #
# dada() function removes substitution and indel errors, but not chimeras.
# But the accuracy of seqs after denoising makes chimera ID easier than when dealing with OTUs.
# Removes seqs that can be reconstructed as bimeras (2 parent chimera) from more abundant seqs

for (met in merge_methods.extras) {
  print(met)
  # Rows are sample names, cols are sequence variants
  all_mergers.mocks[[met]][["seqtab.nochim"]] <- removeBimeraDenovo(all_mergers.mocks[[met]][["seqtab"]], 
                                                                    method="consensus", multithread=TRUE, verbose=TRUE)
  dim(all_mergers.mocks[[met]][["seqtab"]])
  dim(all_mergers.mocks[[met]][["seqtab.nochim"]])
  sum(all_mergers.mocks[[met]][["seqtab.nochim"]])/sum(all_mergers.mocks[[met]][["seqtab"]])
}
# [1] "dada2"
# Identified 287 bimeras out of 418 input sequences.
# [1] "PEAR"
# Identified 326 bimeras out of 523 input sequences.
# [1] "pandaseq_sb"
# Identified 325 bimeras out of 519 input sequences.
# [1] "pandaseq_ea"
# Identified 222 bimeras out of 1686 input sequences.
# [1] "pandaseq_flash"
# Identified 1 bimeras out of 35 input sequences.
# [1] "pandaseq_pear"
# Identified 390 bimeras out of 1099 input sequences.
# [1] "pandaseq_rdp_mle"
# Identified 314 bimeras out of 509 input sequences.
# [1] "pandaseq_stitch"
# Identified 198 bimeras out of 373 input sequences.
# [1] "pandaseq_uparse"
# Identified 171 bimeras out of 340 input sequences.
# [1] "fqj"
# Identified 311 bimeras out of 493 input sequences.
# [1] "dada2_5"
# Identified 509 bimeras out of 702 input sequences.
# [1] "dada2_10"
# Identified 857 bimeras out of 1160 input sequences.
# [1] "dada2_20"
# Identified 1051 bimeras out of 1510 input sequences.








# ****************************************************** #
#### TRACK READS THROUGH THE PIPELINE #### 
# ****************************************************** #
# check number of reads that remained through each step in pipeline
getN <- function(x) sum(getUniques(x))

for (met in merge_methods.extras) {
  print(met)
  if (met %in% c("dada2","dada2_5","dada2_10","dada2_20")) {
    all_mergers.mocks[[met]][["track"]] <- cbind(out.mocks, sapply(all_mergers.mocks[["dada2"]][["dadaFs"]], getN),
                                                 sapply(all_mergers.mocks[[met]][["dada"]], getN),
                                                 rowSums(all_mergers.mocks[[met]][["seqtab"]]),
                                                 rowSums(all_mergers.mocks[[met]][["seqtab.nochim"]]),
                                                 rowSums(all_mergers.mocks[[met]][["seqtab.nochim"]])/out.mocks[,1])
    colnames(all_mergers.mocks[[met]][["track"]]) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim", "%_of_start")
    rownames(all_mergers.mocks[[met]][["track"]]) <- sample.names.mocks
    
  } else {
    all_mergers.mocks[[met]][["track"]] <- cbind(out.mocks, sapply(all_mergers.mocks[[met]][["dada"]], getN),
                                                 rowSums(all_mergers.mocks[[met]][["seqtab"]]),
                                                 rowSums(all_mergers.mocks[[met]][["seqtab.nochim"]]),
                                                 rowSums(all_mergers.mocks[[met]][["seqtab.nochim"]])/out.mocks[,1])
    colnames(all_mergers.mocks[[met]][["track"]]) <- c("input", "filtered", "denoised", "tabled", "nonchim", "%_of_start")
    rownames(all_mergers.mocks[[met]][["track"]]) <- sample.names.mocks
  }
}





# ****************************************************** #
#### ASSIGN TAXONOMY #### 
# ****************************************************** #
# Package provides native implementation of naive Bayesian classifier from the 
#   Ribosomal Database Project (RDP). 
# assignTaxonomy() takes a set of seqs, a training set of taxonomically classified seqs, outputs
#   taxonomic assignments with at least "minBoot" bootstrap confidence

# The dada folks maintain formatted training fastas for RDP training set, GreenGenes clustered
#   at 97% identity and the Silva ref db: https://benjjneb.github.io/dada2/training.html

for (met in merge_methods.extras) {
  print(met)
  all_mergers.mocks[[met]][["tax.128"]] <- assignTaxonomy(all_mergers.mocks[[met]][["seqtab.nochim"]], 
                                                      sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)
  all_mergers.mocks[[met]][["tax.132"]] <- assignTaxonomy(all_mergers.mocks[[met]][["seqtab.nochim"]],
                                                          sprintf("%s/silva_nr_v132_train_set.fa.gz",path.tut), multithread=TRUE)
}


# ***************** #

# can also do species level assignments based on exact matching between ASVs and sequenced reference strains
#   https://benjjneb.github.io/dada2/assign.html#species-assignment

for (met in merge_methods.extras) {
  # Add species
  print(met)
  # split into 2 tables because some may be too large to load into memory all at once
  halfRows.128 <- round(nrow(all_mergers.mocks[[met]][["tax.128"]]) / 2)
  ts.128.p1 <- addSpecies(all_mergers.mocks[[met]][["tax.128"]][ 1:halfRows.128, ], 
                          sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
  ts.128.p2 <- addSpecies(all_mergers.mocks[[met]][["tax.128"]][ (halfRows.128+1):nrow(all_mergers.mocks[[met]][["tax.128"]]), ], 
                          sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
  
  all_mergers.mocks[[met]][["tax_spec.128"]] <- rbind(ts.128.p1, ts.128.p2)
  
  # split into 2 tables because some may be too large to load into memory all at once
  halfRows.132 <- round(nrow(all_mergers.mocks[[met]][["tax.132"]]) / 2)
  ts.132.p1 <- addSpecies(all_mergers.mocks[[met]][["tax.132"]][ 1:halfRows.132, ], 
                          sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))
  ts.132.p2 <- addSpecies(all_mergers.mocks[[met]][["tax.132"]][ (halfRows.132+1):nrow(all_mergers.mocks[[met]][["tax.132"]]), ], 
                          sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))
  
  all_mergers.mocks[[met]][["tax_spec.132"]] <- rbind(ts.132.p1, ts.132.p2)
}


# ***************** #

# inspect taxonomic assignments
for (met in merge_methods.extras) {
  print(met)
  all_mergers.mocks[[met]][["tax_spec.128.print"]] <- all_mergers.mocks[[met]][["tax_spec.128"]]
  all_mergers.mocks[[met]][["tax_spec.132.print"]] <- all_mergers.mocks[[met]][["tax_spec.132"]]
  rownames(all_mergers.mocks[[met]][["tax_spec.128.print"]]) <- NULL
  rownames(all_mergers.mocks[[met]][["tax_spec.132.print"]]) <- NULL
  
  print("tax.128")
  print(summary(is.na(all_mergers.mocks[[met]][["tax_spec.128.print"]][,"Genus"])))
  print(summary(is.na(all_mergers.mocks[[met]][["tax_spec.128.print"]][,"Species"])))
  print("tax.132")
  print(summary(is.na(all_mergers.mocks[[met]][["tax_spec.132.print"]][,"Genus"])))
  print(summary(is.na(all_mergers.mocks[[met]][["tax_spec.132.print"]][,"Species"])))
  print("")
}

# ***************** #










# ****************************************************** #
#### Make OTU tables #### 
# ****************************************************** #

for (met in merge_methods.extras) {
  
  for (silva in c("tax.128","tax.132")) {
    print(c(met, silva))
    # change NAs to unclassified
    tax.fixNames <- all_mergers.mocks[[met]][[silva]]
    tax.fixNames[ is.na(tax.fixNames) ] <- "unclassified"
    
    tax.fixNames <- t(apply(tax.fixNames, 1, function(x) {
      all.levels <- as.character(x)
      if ("unclassified" %in% all.levels) {
        # at each taxlevel, if unclassified, include the names of all levels above it, so we can see how much of its taxonomy we know,
        # and can differentiate unclassified taxa when possible
        new.row <- sapply(2:ncol(tax.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                        paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                        all.levels[tl]))
        new.row <- c(all.levels[1], new.row)
      } else {
        all.levels
      }
    }))
    colnames(tax.fixNames) <- colnames(all_mergers.mocks[[met]][[silva]])
    
    # now that the tax table has updated names for unclassified taxa, get otu tables
    otu.name <- sprintf("%s.otu.genus", silva) # e.g. tax.132.otu.genus
    otu.name.rel <- sprintf("%s.otu.genus.rel", silva)
    
    tmp_otu <- all_mergers.mocks[[met]][["seqtab.nochim"]]
    
    colnames(tmp_otu) <- unname(tax.fixNames[ colnames(all_mergers.mocks[[met]][["seqtab.nochim"]]), "Genus"])
    non.bact <- tax.fixNames[ tax.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus" ]
    tmp_otu <- as.matrix( tmp_otu[ , ! colnames(tmp_otu) %in% non.bact ] )
    # see this post for how to merge all the columns with the same column names:
    #   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
    tmp_otu <- tmp_otu %*% sapply(unique(colnames(tmp_otu)),"==", colnames(tmp_otu))
    tmp_otu.rel <- t(apply(tmp_otu, 1, function(x) round(100 * x/sum(x), 2)))
    tmp_otu.rel[is.nan(tmp_otu.rel)] <- 0
    
    all_mergers.mocks[[met]][[otu.name]] <- tmp_otu
    all_mergers.mocks[[met]][[otu.name.rel]] <- tmp_otu.rel
    print(dim(tmp_otu))
    print(dim(tmp_otu.rel))
    print('')
  }
  
}







# ****************************************************** #
#### Get scores for mock community representation #### 
# ****************************************************** #

# get values for both even and staggered mock communities:
mock.ref <- getSequences(file.path(path.tut, "HMP_MOCK.v35.fasta"))
mock.genera <- c("Acinetobacter","Actinomyces","Bacillus","Bacteroides","Clostridium_sensu_stricto_1","Deinococcus","Enterococcus",
                 "Escherichia/Shigella","Helicobacter","Lactobacillus","Listeria","Neisseria","Propionibacterium","Pseudomonas",
                 "Rhodobacter","Staphylococcus","Staphylococcus","Streptococcus","Streptococcus","Streptococcus")
# changed Escherichia => Escherichia/Shigella and Clostridium => Clostridium_sensu_stricto_1
# because that is how the names appear in the silva database and thus how they are assigned here
mock.genera.copies <- c(10000,1000,100000,1000,100000,1000,1000,1000000,10000,10000,10000,10000,10000,100000,1000000,100000,1000000,
                        100000,1000000,1000)
mock.genera.copies <- sapply(mock.genera.copies, function(x) round(100 * x/sum(mock.genera.copies), 2))
mock.genera.even <- rep(10000, 20)
mock.genera.even <- sapply(mock.genera.even, function(x) round(100 * x/sum(mock.genera.even), 2))
mock.genera.copies <- t(cbind(mock.genera, mock.genera.even, mock.genera.copies))
colnames(mock.genera.copies) <- mock.genera.copies[1,]
mock.genera.copies <- mock.genera.copies[2:3,]
mock.genera.copies <- apply(mock.genera.copies, 2, function(x) as.numeric(x))
rownames(mock.genera.copies) = c("Even","Staggered")
mock.genera.copies <- t(mock.genera.copies %*% sapply(unique(colnames(mock.genera.copies)), "==", colnames(mock.genera.copies)))


# get vector of all genera that appear in all tables
all_gen <- unique(unname(unlist(sapply(merge_methods.extras, function(x) unique(colnames(all_mergers.mocks[[x]][["tax.128.otu.genus.rel"]]))))))
all_gen <- c(all_gen,unique(unname(unlist(sapply(merge_methods.extras, function(x) unique(colnames(all_mergers.mocks[[x]][["tax.132.otu.genus.rel"]])))))))
all_gen <- sort(unique(all_gen))



# ********************************* #
# function to get scores for each type of mock community
get_mock_comp_score <- function(m, meth, o.name) {
  
  # use sum of counts in given mock community as a weight for score
  mock.counts <- sum(all_mergers.mocks[[meth]][[o.name]][ m, ])
  # use relative abundance values to compare to known percentages in mock communities
  mock.rel <- all_mergers.mocks[[meth]][[sprintf("%s.rel",o.name)]] [ m, ]
  # add known percentage for those genera that do not appear in a given sample
  not_present <- unique(mock.genera[ ! mock.genera %in% names(mock.rel) ])
  
  if (startsWith(m, "HM-782")) {
    # for the even mock community
    sco <- sum(sapply(names(mock.rel), function(x) ifelse(x %in% mock.genera,
                                                   # if one of the genera that should be present, get percent difference
                                                   abs(mock.rel[x] - mock.genera.copies[x,"Even"]) / mock.genera.copies[x,"Even"],
                                                   # else it is a false positive, just add its abundance
                                                   mock.rel[x])))
    sco <- sco + sum(mock.genera.copies[not_present,"Even"])
    
  } else if (startsWith(m, "HM-783")) {
    # for the staggered mock community
    sco <- sum(sapply(names(mock.rel), function(x) ifelse(x %in% mock.genera,
                                                   # if one of the genera that should be present, get percent difference
                                                   abs(mock.rel[x] - mock.genera.copies[x,"Staggered"]) / mock.genera.copies[x,"Staggered"],
                                                   # else it is a false positive, just add its abundance
                                                   mock.rel[x])))
    sco <- sco + sum(mock.genera.copies[not_present,"Staggered"])
    
  } else {
    # for the control communities
    sco <- sum(mock.rel)
  }
  
  return(sco * mock.counts)
}
# ********************************* #



# ********************************* #
# produce scores
for (met in merge_methods.extras) {
  for (silva in c("tax.128", "tax.132")) {
    
    otu.name <- sprintf("%s.otu.genus", silva)
    score.name <- sprintf("%s.score", silva)
    
    total_counts <- sum( all_mergers.mocks[[met]][[otu.name]] )
    
    all_mergers.mocks[[met]][[score.name]] <- sum(sapply(sample.names.mocks, get_mock_comp_score, met, otu.name )) / total_counts
    
    print(c( met, silva, all_mergers.mocks[[met]][[score.name]] ))
    
  }
}
# ********************************* #



t128_scores <- sapply(merge_methods.extras, function(x) all_mergers.mocks[[x]][["tax.128.score"]])
t132_scores <- sapply(merge_methods.extras, function(x) all_mergers.mocks[[x]][["tax.132.score"]])
scores <- as.data.frame(cbind(merge_methods.extras, t128_scores, t132_scores))
scores <- scores[ ! rownames(scores) == "pandaseq_flash", ] # exclude flash cause it has a terrible outlier score
scores.m <- melt(scores, id.vars = "merge_methods.extras")
colnames(scores.m) <- c("Method", "Silva", "Score")
scores.m$Score <- as.numeric(scores.m$Score)

ggplot(scores.m, aes(x=Method, y=Score, fill=Silva)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# So the dada2 built-in merging method is apparently by far the best way
# and there is very little difference between the different maxMismatch values for the mergePairs() 
# function, but still the default of 0 has the lowest scores in both v128 and v132 of silva


# save.image(file = "~/Downloads/dada2_benchmark_mocks.RData")






