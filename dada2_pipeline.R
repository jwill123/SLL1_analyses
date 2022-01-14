

library(dada2)
library(ggplot2)
library(phyloseq)

home_dir <- "/users/tg/jwillis/SLL"

### SEE TUTORIAL HERE:   https://benjjneb.github.io/dada2/tutorial.html


# ****************************************************** #
path1 <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/s1-s700"
path2 <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/s701-s1319"
pathMock <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/mock_communities"
length(list.files(path1))
length(list.files(path2))
length(list.files(pathMock))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path1, pattern="_R1_001.fastq", full.names = TRUE))
fnFs <- sort( c( fnFs, sort(list.files(path2, pattern="_R1_001.fastq", full.names = TRUE)) ) )
fnFs <- sort( c( fnFs, sort(list.files(pathMock, pattern="_R1_001.fastq", full.names = TRUE)) ) )
fnRs <- sort(list.files(path1, pattern="_R2_001.fastq", full.names = TRUE))
fnRs <- sort( c( fnRs, sort(list.files(path2, pattern="_R2_001.fastq", full.names = TRUE)) ) )
fnRs <- sort( c( fnRs, sort(list.files(pathMock, pattern="_R2_001.fastq", full.names = TRUE)) ) )
# Extract sample names, assuming filenames have format: datax_xxx_SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3)


# ****************************************************** #
#### EXAMINE QUALITY PROFILES OF FORWARD AND REVERSE READS #### 
# ****************************************************** #
plotQualityProfile(fnFs[31:34]) + 
  scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5))
# From this it seems should truncate forward reads at least at position 250 (trim last 50 nucs)

plotQualityProfile(fnRs[31:34]) +
  scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5))
# reverse have worse quality at ends, truncate at position 200

# **** If using this workflow on your own data ****:
# Your reads must still overlap after truncation in order to merge them later! The tutorial 
#    is using 2x250 V4 sequence data, so the forward and reverse reads almost completely 
#    overlap and our trimming can be completely guided by the quality scores. If you are 
#    using a less-overlapping primer set, like V1-V2 or V3-V4, your truncLen must be large 
#    enough to maintain 20 + biological.length.variation nucleotides of overlap between them.
# Non-overlapping primer sets are supported as well with mergePairs(..., justConcatenate=TRUE) 
#    when performing merging.





# ****************************************************** #
#### PERFORM FILTERING AND TRIMMING #### 
# ****************************************************** #
# filt_path <- file.path("/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission", "filtered_reads") # Place filtered files in filtered_reads/ subdirectory
# filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
# 
# 
# use maxEE=2 cause # of expected errors is better than average quality scores; http://www.drive5.com/usearch/manual/expected_errors.html
# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
#                      truncLen=c(250,200), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                      compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
# 
# mean(out[,2] / out[,1])


# This removed almost half the reads per sample on average
# So here I output filtered reads into a different folder, using a more relaxed maxEE filter for the reverse reads
filt_path <- file.path("/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission", "filtered_reads") # Place filtered files in filtered_reads/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(275,225), maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                     trimLeft=c(25,10), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

mean(out[,2] / out[,1]) # 0.705 -- not bad

# **** If using this workflow on your own data ****:
# The standard filtering parameters are starting points, not set in stone. For example, 
#    if too few reads are passing the filter, considering relaxing maxEE, perhaps especially
#    on the reverse reads (eg. maxEE=c(2,5)). If you want to speed up downstream computation,
#    consider tightening maxEE. For paired-end reads consider the length of your amplicon when 
#    choosing truncLen as your reads must overlap after truncation in order to merge them later. 
# **** If using this workflow on your own data ****:
# For common ITS amplicon strategies, it is undesirable to truncate reads to a fixed length due 
#    to the large amount of length variation at that locus. That is OK, just leave out truncLen.
#    Make sure you removed the forward and reverse primers from both the forward and reverse reads though!








# ****************************************************** #
#### LEARN THE ERROR RATES #### 
# ****************************************************** #
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# when this is run with the filtFs, stops after sample 149 because the 'nbases' parameter is 1e8 reads
#   by default, stops after this many are read into memory. Obviously there is a lot of data remaining,
#   so if the values look weird after running this, maybe should increase the 'nbases' value

# errF.2 <- learnErrors(filtFs.2, multithread=TRUE, nreads = 2e+6)
# errR.2 <- learnErrors(filtRs.2, multithread=TRUE, nreads = 2e+6)

# show error rates for each possible transition (eg. A->C, A->G, …) 
plotErrors(errF, nominalQ=TRUE) # looks pretty good
plotErrors(errR, nominalQ=TRUE) # looks all good


# plotErrors(errF.2, nominalQ=TRUE)
# plotErrors(errR.2, nominalQ=TRUE)


# Points are the observed error rates for each consensus quality score. 
# The black line shows the estimated error rates after convergence. 
# The red line shows the error rates expected under the nominal definition of the Q-value.

# **** If using this workflow on your own data ****:
# Parameter learning is computationally intensive, so by default the learnErrors function 
#    uses only a subset of the data (the first 1M reads). If the plotted error model does 
#    not look like a good fit, try increasing the nreads parameter to see if the fit improves.



# from this post about not converging after 10 rounds or error estimation https://github.com/benjjneb/dada2/issues/77
dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)
# 
# dada2:::checkConvergence(errF.2)
# dada2:::checkConvergence(errR.2)










# ****************************************************** #
#### DEREPLICATION, SAMPLE INFERENCE, MERGE PAIRED READS #### 
# ****************************************************** #
# Combines all identical sequences into unique sequences with abundances
# Reduces computation time by eliminating redundant comparisons
# *** UNIQUE TO DADA2: retains summary of quality info for each unique sequence
#     This is an average of positional qualities for duplicated reads to inform the error model

# Since there are so many samples, better to process them in streaming fashion (for loop)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  # cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[ grep(sam, filtFs) ], verbose = T)
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  
  derepR <- derepFastq(filtRs[ grep(sam, filtRs) ], verbose = T)
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
# rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab.rds")


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# only keep 
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(298,310)]






# mergers.2 <- vector("list", length(sample.names))
# names(mergers.2) <- sample.names
# 
# for(sam in sample.names) {
#   # cat("Processing:", sam, "\n")
#   derepF.2 <- derepFastq(filtFs.2[ grep(sam, filtFs.2) ], verbose = T)
#   ddF.2 <- dada(derepF.2, err=errF.2, multithread=TRUE)
#   
#   derepR.2 <- derepFastq(filtRs.2[ grep(sam, filtRs.2) ], verbose = T)
#   ddR.2 <- dada(derepR.2, err=errR.2, multithread=TRUE)
#   
#   merger.2 <- mergePairs(ddF.2, derepF.2, ddR.2, derepR.2)
#   mergers.2[[sam]] <- merger.2
# }
# # rm(derepF); rm(derepR)
# 
# # Construct sequence table and remove chimeras
# seqtab.2 <- makeSequenceTable(mergers.2)
# saveRDS(seqtab.2, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab.2.rds")
# 
# # Inspect distribution of sequence lengths
# table(nchar(getSequences(seqtab.2)))


# save.image(file = "~/Downloads/dada2_SLL1.RData")

# **** NOTE: According to the tutorial on working with big data (https://benjjneb.github.io/dada2/bigdata.html)
#    Large projects can span multiple sequencing runs, and because different runs can have different error profiles, 
#      it is recommended to learn the error rates for each run individually. Typically this means running the Sample 
#      Inference script once for each run or lane, and then merging those runs together into a full-study sequence table. 
#      If your study is contained on one run, that part of this script can be ignored.

# so now will do the same but grouping samples by the plate on which they were sequenced:
by_plate <- read.delim("/users/tg/jwillis/SLL/Part_1/annick/samples_by_plate_1319.csv")

mergers.plate <- vector("list", ncol(by_plate))
names(mergers.plate) <- sapply(3:17, function(x) sprintf("Plate_%s", x))

for (plate in 3:17) {
  print(sprintf("Plate_%s", plate))
  
  samps <- by_plate[ , sprintf("Plate_%s", plate) ]
  filtFs.plate <- filtFs[ grep(samps, filtFs) ]
  filtRs.plate <- filtRs[ grep(samps, filtRs) ]
  
  derepF.plate <- derepFastq(filtFs.plate, verbose = T)
  names(derepF.plate) <- samps
  ddF.plate <- dada(derepF.plate, err=errF, multithread=TRUE)
  
  derepR.plate <- derepFastq(filtRs.plate, verbose = T)
  names(derepR.plate) <- samps
  ddR.plate <- dada(derepR.plate, err=errR, multithread=TRUE)
  
  merger.plate <- mergePairs(ddF.plate, derepF.plate, ddR.plate, derepR.plate)
  mergers.plate[[sprintf("Plate_%s", plate)]] <- merger.plate
}






# ******** #
# Try for only the forward reads:
# mergers.Fonly <- vector("list", length(sample.names))
# names(mergers.Fonly) <- sample.names
# 
# for(sam in sample.names) {
#   # cat("Processing:", sam, "\n")
#   derepF <- derepFastq(filtFs[ grep(sam, filtFs) ], verbose = T)
#   mergers.Fonly[[sam]] <- dada(derepF, err=errF, multithread=TRUE)
# }
# # rm(derepF); rm(derepR)
# 
# # Construct sequence table and remove chimeras
# seqtab.Fonly <- makeSequenceTable(mergers.Fonly)
# saveRDS(seqtab.Fonly, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab.Fonly.rds")
# 
# 
# # Inspect distribution of sequence lengths
# table(nchar(getSequences(seqtab.Fonly)))





# ******** #
# Try for only the forward reads:
# mergers.Fonly.2 <- vector("list", length(sample.names))
# names(mergers.Fonly.2) <- sample.names
# 
# for(sam in sample.names) {
#   # cat("Processing:", sam, "\n")
#   derepF <- derepFastq(filtFs.2[ grep(sam, filtFs.2) ], verbose = T)
#   mergers.Fonly.2[[sam]] <- dada(derepF, err=errF.2, multithread=TRUE)
# }
# # rm(derepF); rm(derepR)
# 
# # Construct sequence table and remove chimeras
# seqtab.Fonly.2 <- makeSequenceTable(mergers.Fonly.2)
# saveRDS(seqtab.Fonly.2, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab.Fonly.2.rds")
# 
# 
# # Inspect distribution of sequence lengths
# table(nchar(getSequences(seqtab.Fonly.2)))







# ****************************************************** #
#### MERGE SEQUENCE RUNS (if necessary), REMOVE CHIMERAS, ASSIGN TAXONOMY #### 
# ****************************************************** #
# This takes the seqtab (or the merged seqtabs from the different sequencing runs) and removes chimeras
# Then assign taxonomy as usual

# Merge multiple runs (if necessary)
# st1 <- readRDS("/users/tg/jwillis/SLL/Part_1/DADA2/seqtab1.rds")
# st2 <- readRDS("/users/tg/jwillis/SLL/Part_1/DADA2/seqtab2.rds")
# st3 <- readRDS("/users/tg/jwillis/SLL/Part_1/DADA2/seqtab3.rds")
# st.all <- mergeSequenceTables(st1, st2, st3)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
dim(seqtab)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Assign taxonomy
path.tut <- "/users/tg/jwillis/SLL/DADA2_files"
tax <- assignTaxonomy(seqtab.nochim, sprintf("%s/silva_nr_v132_train_set.fa.gz",path.tut), multithread=TRUE)

# Add species
tax_spec <- addSpecies(tax, sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))

# remove "." from any taxa names because it will break part of the pipeline later that relies on splitting names by "."
tax[ grepl('\\.',tax) ] <- gsub('\\.','_',tax[grepl('\\.',tax)])
tax_spec[ grepl('\\.',tax_spec) ] <- gsub('\\.','_',tax_spec[grepl('\\.',tax_spec)])

# Write to disk
saveRDS(seqtab.nochim, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab_final.rds")
saveRDS(tax, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_final.rds")
saveRDS(tax_spec, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_species_final.rds")




# inspect taxonomic assignments
tax.print <- tax # Removing sequence rownames for display only
rownames(tax.print) <- NULL
head(tax.print)

tax_spec.print <- tax_spec # Removing sequence rownames for display only
rownames(tax_spec.print) <- NULL
head(tax_spec.print)

summary(is.na(tax_spec.print[,"Phylum"]))
summary(is.na(tax_spec.print[,"Class"]))
summary(is.na(tax_spec.print[,"Order"]))
summary(is.na(tax_spec.print[,"Family"]))
summary(is.na(tax_spec.print[,"Genus"]))
summary(is.na(tax_spec.print[,"Species"]))







# # Remove chimeras
# seqtab.nochim.2 <- removeBimeraDenovo(seqtab.2, method="consensus", multithread=TRUE)
# dim(seqtab.2)
# dim(seqtab.nochim.2)
# sum(seqtab.nochim.2)/sum(seqtab.2)
# 
# 
# # Assign taxonomy
# path.tut <- "/users/tg/jwillis/SLL/mothur_tutorial/MiSeq_SOP_dada2"
# tax.2 <- assignTaxonomy(seqtab.nochim.2, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)
# 
# # Add species
# halfRows <- round(nrow(tax.2) / 2)
# tax_spec.2.p1 <- addSpecies(tax.2[ 1:halfRows, ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.2.p2 <- addSpecies(tax.2[ (halfRows+1):nrow(tax.2), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# # had to split into 2 tables because too large to load into memory all at once
# tax_spec.2 <- rbind(tax_spec.2.p1, tax_spec.2.p2)
# 
# 
# # Write to disk
# saveRDS(seqtab.nochim.2, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab_final.2.rds")
# saveRDS(tax.2, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_final.2.rds")
# saveRDS(tax_spec.2, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_species_final.2.rds")
# 
# 
# # inspect taxonomic assignments
# tax.print.2 <- tax.2 # Removing sequence rownames for display only
# rownames(tax.print.2) <- NULL
# head(tax.print.2)
# 
# tax_spec.print.2 <- tax_spec.2 # Removing sequence rownames for display only
# rownames(tax_spec.print.2) <- NULL
# head(tax_spec.print.2)
# 
# summary(is.na(tax_spec.print.2[,"Phylum"]))
# summary(is.na(tax_spec.print.2[,"Class"]))
# summary(is.na(tax_spec.print.2[,"Order"]))
# summary(is.na(tax_spec.print.2[,"Family"]))
# summary(is.na(tax_spec.print.2[,"Genus"]))
# summary(is.na(tax_spec.print.2[,"Species"]))








# Remove chimeras
# seqtab.nochim.Fonly <- removeBimeraDenovo(seqtab.Fonly, method="consensus", multithread=TRUE)
# dim(seqtab.Fonly)
# dim(seqtab.nochim.Fonly)
# sum(seqtab.nochim.Fonly)/sum(seqtab.Fonly)
# 
# 
# # Assign taxonomy
# path.tut <- "/users/tg/jwillis/SLL/mothur_tutorial/MiSeq_SOP_dada2"
# tax.Fonly <- assignTaxonomy(seqtab.nochim.Fonly, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)
# 
# # Add species
# q1Rows.Fonly <- round(nrow(tax.Fonly) / 4)
# halfRows.Fonly <- round(nrow(tax.Fonly) / 2)
# q3Rows.Fonly <- round(nrow(tax.Fonly) * 3 / 4)
# tax_spec.Fonly.p1 <- addSpecies(tax.Fonly[ 1:q1Rows.Fonly, ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.Fonly.p2 <- addSpecies(tax.Fonly[ (q1Rows.Fonly+1):halfRows.Fonly, ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.Fonly.p3 <- addSpecies(tax.Fonly[ (halfRows.Fonly+1):q3Rows.Fonly, ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.Fonly.p4 <- addSpecies(tax.Fonly[ (q3Rows.Fonly+1):nrow(tax.Fonly), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# # had to split into 2 tables because too large to load into memory all at once
# tax_spec.Fonly <- rbind(tax_spec.Fonly.p1, tax_spec.Fonly.p2, tax_spec.Fonly.p3, tax_spec.Fonly.p4)
# 
# # Write to disk
# saveRDS(seqtab.nochim.Fonly, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab_final.Fonly.rds")
# saveRDS(tax.Fonly, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_final.Fonly.rds")
# saveRDS(tax_spec.Fonly, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_species_final.Fonly.rds")
# 
# 
# 
# 
# # inspect taxonomic assignments
# tax.print.Fonly <- tax.Fonly # Removing sequence rownames for display only
# rownames(tax.print.Fonly) <- NULL
# head(tax.print.Fonly)
# 
# tax_spec.print.Fonly <- tax_spec.Fonly # Removing sequence rownames for display only
# rownames(tax_spec.print.Fonly) <- NULL
# head(tax_spec.print.Fonly)
# 
# summary(is.na(tax_spec.print.Fonly[,"Phylum"]))
# summary(is.na(tax_spec.print.Fonly[,"Class"]))
# summary(is.na(tax_spec.print.Fonly[,"Order"]))
# summary(is.na(tax_spec.print.Fonly[,"Family"]))
# summary(is.na(tax_spec.print.Fonly[,"Genus"]))
# summary(is.na(tax_spec.print.Fonly[,"Species"]))







# # Remove chimeras
# seqtab.nochim.Fonly.2 <- removeBimeraDenovo(seqtab.Fonly.2, method="consensus", multithread=TRUE)
# dim(seqtab.Fonly.2)
# dim(seqtab.nochim.Fonly.2)
# sum(seqtab.nochim.Fonly.2)/sum(seqtab.Fonly.2)
# 
# 
# # Assign taxonomy
# path.tut <- "/users/tg/jwillis/SLL/mothur_tutorial/MiSeq_SOP_dada2"
# tax.Fonly.2 <- assignTaxonomy(seqtab.nochim.Fonly.2, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)
# 
# # Add species
# halfRows.Fonly.2 <- round(nrow(tax.Fonly.2) / 2)
# tax_spec.Fonly.2.p1 <- addSpecies(tax.Fonly.2[ 1:halfRows.Fonly.2, ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.Fonly.2.p2 <- addSpecies(tax.Fonly.2[ (halfRows.Fonly.2+1):nrow(tax.Fonly.2), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# # had to split into 2 tables because too large to load into memory all at once
# tax_spec.Fonly.2 <- rbind(tax_spec.Fonly.2.p1, tax_spec.Fonly.2.p2)
# 
# 
# # Write to disk
# saveRDS(seqtab.nochim.Fonly.2, "/users/tg/jwillis/SLL/Part_1/DADA2/seqtab_final.Fonly.2.rds")
# saveRDS(tax.Fonly.2, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_final.Fonly.2.rds")
# saveRDS(tax_spec.Fonly.2, "/users/tg/jwillis/SLL/Part_1/DADA2/tax_species_final.Fonly.2.rds")
# 
# 
# 
# 
# # inspect taxonomic assignments
# tax.print.Fonly.2 <- tax.Fonly.2 # Removing sequence rownames for display only
# rownames(tax.print.Fonly.2) <- NULL
# head(tax.print.Fonly.2)
# 
# tax_spec.print.Fonly.2 <- tax_spec.Fonly.2 # Removing sequence rownames for display only
# rownames(tax_spec.print.Fonly.2) <- NULL
# head(tax_spec.print.Fonly.2)
# 
# summary(is.na(tax_spec.print.Fonly.2[,"Phylum"]))
# summary(is.na(tax_spec.print.Fonly.2[,"Class"]))
# summary(is.na(tax_spec.print.Fonly.2[,"Order"]))
# summary(is.na(tax_spec.print.Fonly.2[,"Family"]))
# summary(is.na(tax_spec.print.Fonly.2[,"Genus"]))
# summary(is.na(tax_spec.print.Fonly.2[,"Species"]))





# ****************************************************** #
#### TRACK READS THROUGH THE PIPELINE #### 
# ****************************************************** #
# check number of reads that remained through each step in pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim), rowSums(seqtab.nochim)/out[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "tabled", "nonchim", "%_of_start")
rownames(track) <- sample.names
head(track)
mean(track[,"%_of_start"]) # 0.655









# ****************************************************** #
#### GET OTU TABLE (genus level) #### 
# ****************************************************** #

# change NAs to unclassified
tax.fixNames <- tax
tax.fixNames[is.na(tax.fixNames)] <- "unclassified"
# combine/adjust genus names when appropriate
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Butyrivibrio_2"), "Genus"] <- "Butyrivibrio"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Clostridium_sensu_stricto_1"), "Genus"] <- "Clostridium"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Corynebacterium_1"), "Genus"] <- "Corynebacterium"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Escherichia/Shigella"), "Genus"] <- "Escherichia"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Prevotella_2","Prevotella_6","Prevotella_7","Prevotella_9"), "Genus"] <- "Prevotella"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Ruminococcus_2"), "Genus"] <- "Ruminococcus"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Selenomonas_3","Selenomonas_4"), "Genus"] <- "Selenomonas"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Spirochaeta_2"), "Genus"] <- "Spirochaeta"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Treponema_2"), "Genus"] <- "Treponema"

# add other levels to unclassified genera to have unique values
tax.fixNames <- t(apply(tax.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                     paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                     all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.fixNames) <- colnames(tax)
# get otu table
otu.table <- seqtab.nochim
colnames(otu.table) <- unname(tax.fixNames[colnames(seqtab.nochim) , "Genus"])
non.bact <- tax.fixNames[ tax.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota","Archaea"), "Genus"]
otu.table <- otu.table[ , ! colnames(otu.table) %in% non.bact ]
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table <- t(otu.table %*% sapply(unique(colnames(otu.table)),"==", colnames(otu.table)))
otu.table.rel <- apply(otu.table, 2, function(x) 100 * x/sum(x))
otu.table.rel[is.nan(otu.table.rel)] <- 0

write.csv(otu.table, sprintf("%s/Part_1/DADA2/otu_table.csv", home_dir))
write.csv(otu.table.rel, sprintf("%s/Part_1/DADA2/otu_table_rel.csv", home_dir))


# get tax table with rownames as genus values
tax.table <- unique(tax.fixNames)
rownames(tax.table) <- unname(tax.table[,"Genus"])
tax.table <- tax.table[ ! rownames(tax.table) %in% non.bact, ]
tax.table <- tax.table[ sort(rownames(tax.table)), ]
write.csv(tax.table, sprintf("%s/Part_1/DADA2/tax_table.csv", home_dir))

   
# *********** #


# change NAs to unclassified
# tax.2.fixNames <- tax.2
# tax.2.fixNames[is.na(tax.2.fixNames)] <- "unclassified"
# tax.2.fixNames <- t(apply(tax.2.fixNames, 1, function(x) {
#   all.levels <- as.character(x)
#   if ("unclassified" %in% all.levels) {
#     new.row <- sapply(2:ncol(tax.2.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
#                                                                 paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
#                                                                 all.levels[tl]))
#     new.row <- c(all.levels[1], new.row)
#   } else {
#     all.levels
#   }
#   
# }))
# colnames(tax.2.fixNames) <- colnames(tax.2)
# # get otu table
# otu.table.2 <- seqtab.nochim.2
# colnames(otu.table.2) <- unname(tax.2.fixNames[colnames(seqtab.nochim.2) , "Genus"])
# non.bact.2 <- tax.2.fixNames[ tax.2.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
# otu.table.2 <- otu.table.2[ , ! colnames(otu.table.2) %in% non.bact.2]
# # see this post for how to merge all the columns with the same column names:
# #   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
# otu.table.2 <- t(otu.table.2 %*% sapply(unique(colnames(otu.table.2)),"==", colnames(otu.table.2)))
# otu.table.2.rel <- apply(otu.table.2, 2, function(x) 100 * x/sum(x))




# *********** #


# change NAs to unclassified
# tax.Fonly.fixNames <- tax.Fonly
# tax.Fonly.fixNames[is.na(tax.Fonly.fixNames)] <- "unclassified"
# tax.Fonly.fixNames <- t(apply(tax.Fonly.fixNames, 1, function(x) {
#   all.levels <- as.character(x)
#   if ("unclassified" %in% all.levels) {
#     new.row <- sapply(2:ncol(tax.Fonly.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
#                                                                 paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
#                                                                 all.levels[tl]))
#     new.row <- c(all.levels[1], new.row)
#   } else {
#     all.levels
#   }
#   
# }))
# colnames(tax.Fonly.fixNames) <- colnames(tax.Fonly)
# # get otu table
# otu.table.Fonly <- seqtab.nochim.Fonly
# colnames(otu.table.Fonly) <- unname(tax.Fonly.fixNames[colnames(seqtab.nochim.Fonly) , "Genus"])
# non.bact.Fonly <- tax.Fonly.fixNames[ tax.Fonly.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
# otu.table.Fonly <- otu.table.Fonly[ , ! colnames(otu.table.Fonly) %in% non.bact.Fonly]
# # see this post for how to merge all the columns with the same column names:
# #   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
# otu.table.Fonly <- t(otu.table.Fonly %*% sapply(unique(colnames(otu.table.Fonly)),"==", colnames(otu.table.Fonly)))
# otu.table.Fonly.rel <- apply(otu.table.Fonly, 2, function(x) 100 * x/sum(x))


# *********** #

# change NAs to unclassified
# tax.Fonly.2.fixNames <- tax.Fonly.2
# tax.Fonly.2.fixNames[is.na(tax.Fonly.2.fixNames)] <- "unclassified"
# tax.Fonly.2.fixNames <- t(apply(tax.Fonly.2.fixNames, 1, function(x) {
#   all.levels <- as.character(x)
#   if ("unclassified" %in% all.levels) {
#     new.row <- sapply(2:ncol(tax.Fonly.2.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
#                                                                       paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
#                                                                       all.levels[tl]))
#     new.row <- c(all.levels[1], new.row)
#   } else {
#     all.levels
#   }
#   
# }))
# colnames(tax.Fonly.2.fixNames) <- colnames(tax.Fonly.2)
# # get otu table
# otu.table.Fonly.2 <- seqtab.nochim.Fonly.2
# colnames(otu.table.Fonly.2) <- unname(tax.Fonly.2.fixNames[colnames(seqtab.nochim.Fonly.2) , "Genus"])
# non.bact.Fonly.2 <- tax.Fonly.2.fixNames[ tax.Fonly.2.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
# otu.table.Fonly.2 <- otu.table.Fonly.2[ , ! colnames(otu.table.Fonly.2) %in% non.bact.Fonly.2]
# # see this post for how to merge all the columns with the same column names:
# #   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
# otu.table.Fonly.2 <- t(otu.table.Fonly.2 %*% sapply(unique(colnames(otu.table.Fonly.2)),"==", colnames(otu.table.Fonly.2)))
# otu.table.Fonly.2.rel <- apply(otu.table.Fonly.2, 2, function(x) 100 * x/sum(x))


# *********** #







# ****************************************************** #
#### GET OTU TABLE (species level) #### 
# ****************************************************** #

# change NAs to unclassified
tax.fixNames.spec <- tax_spec
tax.fixNames.spec[is.na(tax.fixNames.spec)] <- "unclassified"
# combine/adjust genus names when appropriate
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Butyrivibrio_2"), "Genus"] <- "Butyrivibrio"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Clostridium_sensu_stricto_1"), "Genus"] <- "Clostridium"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Corynebacterium_1"), "Genus"] <- "Corynebacterium"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Escherichia/Shigella"), "Genus"] <- "Escherichia"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Prevotella_2","Prevotella_6","Prevotella_7","Prevotella_9"), "Genus"] <- "Prevotella"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Ruminococcus_2"), "Genus"] <- "Ruminococcus"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Selenomonas_3","Selenomonas_4"), "Genus"] <- "Selenomonas"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Spirochaeta_2"), "Genus"] <- "Spirochaeta"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Treponema_2"), "Genus"] <- "Treponema"
# combine/adjust species names when appropriate
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Butyrivibrio_2"), 
                   "Species"] <- "unclassified.Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Butyrivibrio"
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae_1.Clostridium_sensu_stricto_1"), 
              "Species"] <- "unclassified.Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae_1.Clostridium"
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Actinobacteria.Actinobacteria.Corynebacteriales.Corynebacteriaceae.Corynebacterium_1"), 
              "Species"] <- "unclassified.Bacteria.Actinobacteria.Actinobacteria.Corynebacteriales.Corynebacteriaceae.Corynebacterium"
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae.Escherichia/Shigella"), 
              "Species"] <- "unclassified.Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae.Escherichia"
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella_2",
                                            "unclassified.Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella_6",
                                            "unclassified.Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella_7",
                                            "unclassified.Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella_9"), 
              "Species"] <- "unclassified.Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella"
# tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("Ruminococcus_2"), "Species"] <- "Ruminococcus"
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Firmicutes.Negativicutes.Selenomonadales.Veillonellaceae.Selenomonas_3",
                                              "unclassified.Bacteria.Firmicutes.Negativicutes.Selenomonadales.Veillonellaceae.Selenomonas_4"), 
              "Species"] <- "unclassified.Bacteria.Firmicutes.Negativicutes.Selenomonadales.Veillonellaceae.Selenomonas"
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Spirochaetae.Spirochaetes.Spirochaetales.Spirochaetaceae.Spirochaeta_2"), 
              "Species"] <- "unclassified.Bacteria.Spirochaetae.Spirochaetes.Spirochaetales.Spirochaetaceae.Spirochaeta"
tax.fixNames.spec[ tax.fixNames.spec[,"Species"] %in% c("unclassified.Bacteria.Spirochaetae.Spirochaetes.Spirochaetales.Spirochaetaceae.Treponema_2"), 
              "Species"] <- "unclassified.Bacteria.Spirochaetae.Spirochaetes.Spirochaetales.Spirochaetaceae.Treponema"


# make species value scientific name by combining with genus name if not unclassified
tax.fixNames.spec[,"Species"] <- sapply(rownames(tax.fixNames.spec), function(ro) {
  if (tax.fixNames.spec[ro,"Species"] == "unclassified") {
    tax.fixNames.spec[ro,"Species"]
  } else {
    paste(tax.fixNames.spec[ro,"Genus"], tax.fixNames.spec[ro,"Species"], sep = " ")
  }
})

# add other levels to unclassified genera to have unique values
tax.fixNames.spec <- t(apply(tax.fixNames.spec, 1, function(x) {
  all.levels.spec <- as.character(x)
  if ("unclassified" %in% all.levels.spec) {
    new.row.spec <- sapply(2:ncol(tax.fixNames.spec), function(tl) ifelse(all.levels.spec[tl]=="unclassified", 
                                                                          paste(c("unclassified", all.levels.spec[1:(tl-1)]), collapse = '.'),
                                                                          all.levels.spec[tl]))
    new.row.spec <- c(all.levels.spec[1], new.row.spec)
  } else {
    all.levels.spec
  }
  
}))
colnames(tax.fixNames.spec) <- colnames(tax_spec)
# get otu.spec table
otu.spec.table <- seqtab.nochim
colnames(otu.spec.table) <- unname(tax.fixNames.spec[colnames(seqtab.nochim) , "Species"])
non.bact.spec <- tax.fixNames.spec[ tax.fixNames.spec[,"Kingdom"] %in% c(NA, "Eukaryota","Archaea"), "Species"]
otu.spec.table <- otu.spec.table[ , ! colnames(otu.spec.table) %in% non.bact.spec]
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.spec.table <- t(otu.spec.table %*% sapply(unique(colnames(otu.spec.table)),"==", colnames(otu.spec.table)))
otu.spec.table.rel <- apply(otu.spec.table, 2, function(x) 100 * x/sum(x))
otu.spec.table.rel[is.nan(otu.spec.table.rel)] <- 0

write.csv(otu.spec.table, sprintf("%s/Part_1/DADA2/otu_spec_table.csv", home_dir))
write.csv(otu.spec.table.rel, sprintf("%s/Part_1/DADA2/otu_spec_table_rel.csv", home_dir))


# get tax table with rownames as species values
tax.spec.table <- unique(tax.fixNames.spec)
rownames(tax.spec.table) <- unname(tax.spec.table[,"Species"])
tax.spec.table <- tax.spec.table[ ! rownames(tax.spec.table) %in% non.bact.spec, ]
tax.spec.table <- tax.spec.table[ sort(rownames(tax.spec.table)), ]
write.csv(tax.spec.table, sprintf("%s/Part_1/DADA2/tax_spec_table.csv", home_dir))



# *********** #





# save.image(file = "~/Downloads/dada2_SLL1.RData")
# # save.image(file = "~/Downloads/dada2_all.RData")






# ****************************************************** #
#### EVALUATE ACCURACY #### 
# ****************************************************** #
# In the mock community here, mixture of 20 known strains (ref seqs in HMP_MOCK.v35.fasta)
# Now compare seq variants produced here to expected composition of mock community

mocks <- sample.names[ ! grepl('SLL-', sample.names) ]
mock.ref <- getSequences(file.path(path.tut, "HMP_MOCK.v35.fasta"))
# mock.genera <- c("Acinetobacter","Actinomyces","Bacillus","Bacteroides","Clostridium_sensu_stricto_1","Deinococcus","Enterococcus",
#                  "Escherichia/Shigella","Helicobacter","Lactobacillus","Listeria","Neisseria","Propionibacterium","Pseudomonas",
#                  "Rhodobacter","Staphylococcus","Staphylococcus","Streptococcus","Streptococcus","Streptococcus")
# # changed Escherichia => Escherichia/Shigella and Clostridium => Clostridium_sensu_stricto_1
# # because that is how the names appear in the silva database and thus how they are assigned here
mock.genera <- c("Acinetobacter","Actinomyces","Bacillus","Bacteroides","Clostridium","Deinococcus","Enterococcus",
                 "Escherichia","Helicobacter","Lactobacillus","Listeria","Neisseria","Propionibacterium","Pseudomonas",
                 "Rhodobacter","Staphylococcus","Staphylococcus","Streptococcus","Streptococcus","Streptococcus")

# values found here: http://downloads.ihmpdcc.org/data/HMMC/HMPRP_sT1-Mock.pdf
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




unqs.mock.list <- list()
for (m in mocks) {
  unqs.mock <- otu.table[,m]
  m.unqs <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
  unqs.mock.list[[m]] <- m.unqs
  
  cat("DADA2 inferred", length(m.unqs), "sample sequences present in the Mock community:", m, ".\n")
  
  if (length(m.unqs) > 0) {
    # match.ref <- sum(sapply(names(m.unqs), function(x) any(grepl(x, mock.ref))))
    match.ref <- sum(names(m.unqs) %in% mock.genera)
    cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")
  } else{
    cat("\n")
  }
}

# Residual error rate here is 0%. Mothur found 34 OTUs...so far fewer FPs here


# Assign taxonomy to mock samples
snm <- seqtab.nochim[mocks, ]
snm <- snm[ , colSums(snm) > 0]
tax.mocks <- assignTaxonomy(snm, sprintf("%s/silva_nr_v132_train_set.fa.gz",path.tut), multithread=TRUE)

# Add species
tax_spec.mocks <- addSpecies(tax.mocks, sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))

tax.mocks[ grepl('\\.',tax.mocks) ] <- gsub('\\.','_',tax.mocks[grepl('\\.',tax.mocks)])
tax_spec.mocks[ grepl('\\.',tax_spec.mocks) ] <- gsub('\\.','_',tax_spec.mocks[grepl('\\.',tax_spec.mocks)])





# unqs.mock.list.2 <- list()
# for (m in mocks) {
#   unqs.mock <- otu.table.2[,m]
#   m.unqs <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
#   unqs.mock.list.2[[m]] <- m.unqs
#   
#   cat("DADA2 inferred", length(m.unqs), "sample sequences present in the Mock community:", m, ".\n")
#   
#   if (length(m.unqs) > 0) {
#     # match.ref <- sum(sapply(names(m.unqs), function(x) any(grepl(x, mock.ref))))
#     match.ref <- sum(names(m.unqs) %in% mock.genera)
#     cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")
#   } else{
#     cat("\n")
#   }
# }
# 
# # m <- "HM-782D-plate1"
# m <- "HM-782D-plate4"
# # m <- "HM-783D-plate6"
# unqs.mock <- otu.table.2[,m]
# m1.full <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
# m1.full.rel <- sapply(m1.full, function(x) round(100 * x/sum(m1.full), 2))
# m1.full.both <- cbind(m1.full, m1.full.rel)
# colnames(m1.full.both) <- c("Counts", "Rel_abund")
# m1.full.both <- m1.full.both [ sort(rownames(m1.full.both)), ]
# 
# 
# # Assign taxonomy to mock samples
# snm.2 <- seqtab.nochim.2[mocks, ]
# snm.2 <- snm.2[ , colSums(snm.2) > 0]
# tax.mocks.2 <- assignTaxonomy(snm.2, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)
# 
# # Add species
# tax_spec.mocks.2.p1 <- addSpecies(tax.mocks.2[ 1:(nrow(tax.mocks.2) / 2), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.mocks.2.p2 <- addSpecies(tax.mocks.2[ ((nrow(tax.mocks.2) / 2)+1):nrow(tax.mocks.2), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# # had to split into 2 tables because too large to load into memory all at once
# tax_spec.mocks.2 <- rbind(tax_spec.mocks.2.p1, tax_spec.mocks.2.p2)











# unqs.mock.list.Fonly <- list()
# for (m in mocks) {
#   unqs.mock <- otu.table.Fonly[,m]
#   m.unqs <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
#   unqs.mock.list.Fonly[[m]] <- m.unqs
#   
#   cat("DADA2 inferred", length(m.unqs), "sample sequences present in the Mock community:", m, ".\n")
#   
#   if (length(m.unqs) > 0) {
#     # match.ref <- sum(sapply(names(m.unqs), function(x) any(grepl(x, mock.ref))))
#     match.ref <- sum(names(m.unqs) %in% mock.genera)
#     cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")
#   } else{
#     cat("\n")
#   }
# }
# 
# # Assign taxonomy to mock samples
# snm.Fonly <- seqtab.nochim.Fonly[mocks, ]
# snm.Fonly <- snm.Fonly[ , colSums(snm.Fonly) > 0]
# tax.mocks.Fonly <- assignTaxonomy(snm.Fonly, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)
# 
# # Add species
# tax_spec.mocks.Fonly.p1 <- addSpecies(tax.mocks.Fonly[ 1:(nrow(tax.mocks.Fonly) / 2), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.mocks.Fonly.p2 <- addSpecies(tax.mocks.Fonly[ ((nrow(tax.mocks.Fonly) / 2)+1):nrow(tax.mocks.Fonly), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# # had to split into 2 tables because too large to load into memory all at once
# tax_spec.mocks.Fonly <- rbind(tax_spec.mocks.Fonly.p1, tax_spec.mocks.Fonly.p2)






# unqs.mock.list.Fonly.2 <- list()
# for (m in mocks) {
#   unqs.mock <- otu.table.Fonly.2[,m]
#   m.unqs <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
#   unqs.mock.list.Fonly.2[[m]] <- m.unqs
#   
#   cat("DADA2 inferred", length(m.unqs), "sample sequences present in the Mock community:", m, ".\n")
#   
#   if (length(m.unqs) > 0) {
#     # match.ref <- sum(sapply(names(m.unqs), function(x) any(grepl(x, mock.ref))))
#     match.ref <- sum(names(m.unqs) %in% mock.genera)
#     cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")
#   } else{
#     cat("\n")
#   }
# }
# 
# # m <- "HM-782D-plate1"
# m <- "HM-782D-plate4"
# # m <- "HM-783D-plate6"
# unqs.mock <- otu.table.Fonly.2[,m]
# m1.forward <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
# m1.forward.rel <- sapply(m1.forward, function(x) round(100 * x/sum(m1.forward), 2))
# m1.forward.both <- cbind(m1.forward, m1.forward.rel)
# colnames(m1.forward.both) <- c("Counts", "Rel_abund")
# m1.forward.both <- m1.forward.both [ sort(rownames(m1.forward.both)), ]
# 
# # Assign taxonomy to mock samples
# snm.Fonly.2 <- seqtab.nochim.Fonly.2[mocks, ]
# snm.Fonly.2 <- snm.Fonly.2[ , colSums(snm.Fonly.2) > 0]
# tax.mocks.Fonly.2 <- assignTaxonomy(snm.Fonly.2, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)
# 
# # Add species
# tax_spec.mocks.Fonly.2.p1 <- addSpecies(tax.mocks.Fonly.2[ 1:(nrow(tax.mocks.Fonly.2) / 2), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# tax_spec.mocks.Fonly.2.p2 <- addSpecies(tax.mocks.Fonly.2[ ((nrow(tax.mocks.Fonly.2) / 2)+1):nrow(tax.mocks.Fonly.2), ], sprintf("%s/silva_species_assignment_v128.fa.gz",path.tut))
# # had to split into 2 tables because too large to load into memory all at once
# tax_spec.mocks.Fonly.2 <- rbind(tax_spec.mocks.Fonly.2.p1, tax_spec.mocks.Fonly.2.p2)











# ********************************* #
# function to get scores for each type of mock community
get_mock_comp_score <- function(m, ot, otr, pipeline) {
  
  # use sum of counts in given mock community as a weight for score
  mock.counts <- sum(ot[ , m ])
  # use relative abundance values to compare to known percentages in mock communities
  mock.rel <- otr[ , m ]
  # add known percentage for those genera that do not appear in a given sample
  not_present <- unique(mock.genera[ ! mock.genera %in% names(mock.rel) ])
  
  # if (startsWith(m, "HM-782")) {
  if ((pipeline=="mothur" & grepl("HM782", m)) | (pipeline=="dada2" & startsWith(m, "HM-782")) ) {
    # for the even mock community
    sco <- sum(sapply(names(mock.rel), function(x) {
      ifelse(x %in% mock.genera,
             # if one of the genera that should be present, get percent difference
             abs(mock.rel[x] - mock.genera.copies[x,"Even"]) / mock.genera.copies[x,"Even"],
             # else it is a false positive, just add its abundance
             mock.rel[x])
    }))
    sco <- sco + sum(mock.genera.copies[not_present,"Even"])
    
  # } else if (startsWith(m, "HM-783")) {
  } else if ((pipeline=="mothur" & grepl("HM783", m)) | (pipeline=="dada2" & startsWith(m, "HM-783")) ) {
    # for the staggered mock community
    sco <- sum(sapply(names(mock.rel), function(x) {
      ifelse(x %in% mock.genera,
             # if one of the genera that should be present, get percent difference
             abs(mock.rel[x] - mock.genera.copies[x,"Staggered"]) / mock.genera.copies[x,"Staggered"],
             # else it is a false positive, just add its abundance
             mock.rel[x])
    }))
    sco <- sco + sum(mock.genera.copies[not_present,"Staggered"])
    
  } else {
    # for the control communities
    sco <- sum(mock.rel)
  }
  print(c(m, sco * mock.counts))
  return(sco * mock.counts)
}
# ********************************* #


# ****************** #
# ****** get the data.tax table from the SLL_phyloseq.R script
mocks.mothur <- colnames(data.tax)[grepl("HM", colnames(data.tax)) | 
                                     grepl("NTC", colnames(data.tax)) | 
                                     grepl("H2O", colnames(data.tax))]

dt <- data.tax[, mocks.mothur]
dtr <- apply(dt, 2, function(x) 100*(x/sum(x)))
# ****************** #




total_counts <- sum( otu.table )

# score is a comparison to the known values in the mock community samples
#   calculates percent difference for the 17 genera that should appear in these samples
#   adds full value of any mock genus that does not appear in the sample (as a penalty)
#   then also adds value for any genus that should not be in the sample (as a penalty)
# So the lower the score the better
score.dada2 <- sum(sapply(mocks, get_mock_comp_score, otu.table, otu.table.rel, "dada2" )) / total_counts # 0.1547631
score.mothur <- sum(sapply(mocks.mothur, get_mock_comp_score, dt, dtr, "mothur" )) / sum( dt ) # 42.93582





# ********************************* #
# function to get scores for each type of mock community
get_proportion_non_mocks <- function(m, ot, pipeline) {
  # 
  mock.counts <- ot[ , m ]

  # if (startsWith(m, "HM-78")) {
  if ((pipeline=="mothur" & grepl("HM78", m)) | (pipeline=="dada2" & startsWith(m, "HM-78")) ) {
    
    return(list(total=sum(mock.counts),
                non.mock=sum(mock.counts[! names(mock.counts) %in% mock.genera]),
                prop.non.mock=sum(mock.counts[! names(mock.counts) %in% mock.genera]) / sum(mock.counts)))
    # return(list(sum(mock.counts),
    #             sum(mock.counts[! names(mock.counts) %in% mock.genera]),
    #             sum(mock.counts[! names(mock.counts) %in% mock.genera]) / sum(mock.counts)))
  } else {
    
    return(list(total=sum(mock.counts), 
                non.mock=sum(mock.counts),
                prop.non.mock=1))
  }
}
# ********************************* #



mocks.mothur <- colnames(data.tax)[grepl("HM", colnames(data.tax)) | 
                                     grepl("NTC", colnames(data.tax)) | 
                                     grepl("H2O", colnames(data.tax))]

dt <- data.tax[, mocks.mothur]



HM.782 <- sapply(mocks.mothur[grepl("HM782", mocks.mothur)], get_proportion_non_mocks, dt, "mothur")
HM.783 <- sapply(mocks.mothur[grepl("HM783", mocks.mothur)], get_proportion_non_mocks, dt, "mothur")
HM.all <- sapply(mocks.mothur[grepl("HM78", mocks.mothur)], get_proportion_non_mocks, dt, "mothur")
NTC <- sapply(mocks.mothur[! grepl("HM78", mocks.mothur)], get_proportion_non_mocks, dt, "mothur")

# HM.782 <- sapply(mocks[startsWith(mocks, 'HM-782')], get_proportion_non_mocks, otu.table, "dada2")
# HM.783 <- sapply(mocks[startsWith(mocks, 'HM-783')], get_proportion_non_mocks, otu.table, "dada2")
# HM.all <- sapply(mocks[startsWith(mocks, 'HM-78')], get_proportion_non_mocks, otu.table, "dada2")
# NTC <- sapply(mocks[! startsWith(mocks, 'HM-78')], get_proportion_non_mocks, otu.table, "dada2")

HM.782.means <- rowMeans(apply(HM.782, 2, as.numeric)); names(HM.782.means) <- c("total","non.mock","prop.non.mock")
HM.783.means <- rowMeans(apply(HM.783, 2, as.numeric)); names(HM.783.means) <- c("total","non.mock","prop.non.mock")
HM.all.means <- rowMeans(apply(HM.all, 2, as.numeric)); names(HM.all.means) <- c("total","non.mock","prop.non.mock")
NTC.means <- rowMeans(apply(NTC, 2, as.numeric)); names(NTC.means) <- c("total","non.mock","prop.non.mock")

props.mocks <- t(data.frame(HM.782.means=HM.782.means, HM.783.means=HM.783.means, 
                            HM.all.means=HM.all.means, NTC.means=NTC.means))







# ********************************* #
# function to get counts or abundances of non-mock genera
get_non_mocks_genera <- function(m, ot, pipeline) {
  # 
  mock.counts <- ot[ , m ]
  
  # if (startsWith(m, "HM-78")) {
  if ((pipeline=="mothur" & grepl("HM78", m)) | (pipeline=="dada2" & startsWith(m, "HM-78")) ) {
    
    return(mock.counts[! names(mock.counts) %in% mock.genera][mock.counts[! names(mock.counts) %in% mock.genera] != 0])
    
  } else {
    
    return(mock.counts[mock.counts != 0])
    
  }
}
# ********************************* #







prop_or_count <- "prop"
# prop_or_count <- "count"

if (prop_or_count == "prop") {
  HM.782.nmg.mo <- sapply(mocks.mothur[grepl("HM782", mocks.mothur)], get_non_mocks_genera, as.matrix(dtr), "mothur")
  HM.783.nmg.mo <- sapply(mocks.mothur[grepl("HM783", mocks.mothur)], get_non_mocks_genera, as.matrix(dtr), "mothur")
  HM.all.nmg.mo <- sapply(mocks.mothur[grepl("HM78", mocks.mothur)], get_non_mocks_genera, as.matrix(dtr), "mothur")
  NTC.nmg.mo <- sapply(mocks.mothur[! grepl("HM78", mocks.mothur)], get_non_mocks_genera, as.matrix(dtr), "mothur")
  
  HM.782.nmg.da <- sapply(mocks[startsWith(mocks, 'HM-782')], get_non_mocks_genera, as.matrix(otu.table.rel), "dada2")
  HM.783.nmg.da <- sapply(mocks[startsWith(mocks, 'HM-783')], get_non_mocks_genera, as.matrix(otu.table.rel), "dada2")
  HM.all.nmg.da <- sapply(mocks[startsWith(mocks, 'HM-78')], get_non_mocks_genera, as.matrix(otu.table.rel), "dada2")
  NTC.nmg.da <- sapply(mocks[ ! startsWith(mocks, 'HM-78')], get_non_mocks_genera, as.matrix(otu.table.rel), "dada2")
  
} else if (prop_or_count == "count") {
  HM.782.nmg <- sapply(mocks.mothur[grepl("HM782", mocks.mothur)], get_non_mocks_genera, as.matrix(dt), "mothur")
  HM.783.nmg <- sapply(mocks.mothur[grepl("HM783", mocks.mothur)], get_non_mocks_genera, as.matrix(dt), "mothur")
  HM.all.nmg <- sapply(mocks.mothur[grepl("HM78", mocks.mothur)], get_non_mocks_genera, as.matrix(dt), "mothur")
  NTC.nmg <- sapply(mocks.mothur[! grepl("HM78", mocks.mothur)], get_non_mocks_genera, as.matrix(dt), "mothur")
  
  HM.782.nmg.da <- sapply(mocks[startsWith(mocks, 'HM-782')], get_non_mocks_genera, as.matrix(otu.table), "dada2")
  HM.783.nmg.da <- sapply(mocks[startsWith(mocks, 'HM-783')], get_non_mocks_genera, as.matrix(otu.table), "dada2")
  HM.all.nmg.da <- sapply(mocks[startsWith(mocks, 'HM-78')], get_non_mocks_genera, as.matrix(otu.table), "dada2")
  NTC.nmg.da <- sapply(mocks[ ! startsWith(mocks, 'HM-78')], get_non_mocks_genera, as.matrix(otu.table), "dada2")
  
}



get_non_mocks_table <- function(gen) {
  m <- matrix(NA, nrow = length(unique(unlist(lapply(gen, function(x) names(x))))), ncol = length(gen))
  rownames(m) <- unique(unlist(lapply(gen, function(x) names(x))))
  colnames(m) <- names(gen)
  
  for (ro in rownames(m)) {
    for (co in colnames(m)) {
      m[ro, co] <- gen[[co]][ro]
    }
  }
  m[is.na(m)] <- 0
  return(m)
}


HM.782.meanProp.nmg.mo <- sort(rowMeans(get_non_mocks_table(HM.782.nmg.mo)))
HM.783.meanProp.nmg.mo <- sort(rowMeans(get_non_mocks_table(HM.783.nmg.mo)))
HM.all.meanProp.nmg.mo <- sort(rowMeans(get_non_mocks_table(HM.all.nmg.mo)))
NTC.meanProp.nmg.mo <- sort(rowMeans(get_non_mocks_table(NTC.nmg.mo)))

HM.782.meanProp.nmg.da <- sort(rowMeans(get_non_mocks_table(HM.782.nmg.da)))
HM.783.meanProp.nmg.da <- sort(rowMeans(get_non_mocks_table(HM.783.nmg.da)))
HM.all.meanProp.nmg.da <- sort(rowMeans(get_non_mocks_table(HM.all.nmg.da)))
NTC.meanProp.nmg.da <- sort(rowMeans(get_non_mocks_table(NTC.nmg.da)))

# ******************************************* #
prop.da <- matrix(NA, nrow = length(unique(c(names(HM.782.meanProp.nmg.da),names(HM.783.meanProp.nmg.da),
                                             names(HM.all.meanProp.nmg.da),names(NTC.meanProp.nmg.da)))) + 1, ncol = 4)
rownames(prop.da) <- unique(c(mock.genera[mock.genera %in% names(NTC.meanProp.nmg.da)], # included in order to put these at top of table
                              names(HM.782.meanProp.nmg.da),names(HM.783.meanProp.nmg.da),
                              names(HM.all.meanProp.nmg.da),names(NTC.meanProp.nmg.da), 
                              sprintf("Total mean %s", prop_or_count)))
colnames(prop.da) <- c("HM.782.meanProp.nmg.da","HM.783.meanProp.nmg.da","HM.all.meanProp.nmg.da","NTC.meanProp.nmg.da")

for (ro in rownames(prop.da)) {
  if (ro == sprintf("Total mean %s", prop_or_count)) {
    prop.da[ro, "HM.782.meanProp.nmg.da"] <- sum(HM.782.meanProp.nmg.da)
    prop.da[ro, "HM.783.meanProp.nmg.da"] <- sum(HM.783.meanProp.nmg.da)
    prop.da[ro, "HM.all.meanProp.nmg.da"] <- sum(HM.all.meanProp.nmg.da)
    prop.da[ro, "NTC.meanProp.nmg.da"] <- sum(NTC.meanProp.nmg.da)
  } else {
    for (co in colnames(prop.da)) {
      if (co == "HM.782.meanProp.nmg.da") {
        prop.da[ro, co] <- HM.782.meanProp.nmg.da[ro]
      } else if (co == "HM.783.meanProp.nmg.da") {
        prop.da[ro, co] <- HM.783.meanProp.nmg.da[ro]
      } else if (co == "HM.all.meanProp.nmg.da") {
        prop.da[ro, co] <- HM.all.meanProp.nmg.da[ro]
      } else if (co == "NTC.meanProp.nmg.da") {
        prop.da[ro, co] <- NTC.meanProp.nmg.da[ro]
      }
    }
  }
}

write.csv(prop.da, sprintf("/users/tg/jwillis/SLL/Part_1/SLL1_paper/extras/%s.da.csv", prop_or_count))
# ******************************************* #

prop.mo <- matrix(NA, nrow = length(unique(c(names(HM.782.meanProp.nmg.mo),names(HM.783.meanProp.nmg.mo),
                                             names(HM.all.meanProp.nmg.mo),names(NTC.meanProp.nmg.mo)))) + 1, ncol = 4)
rownames(prop.mo) <- unique(c(mock.genera[mock.genera %in% names(NTC.meanProp.nmg.mo)], # included in order to put these at top of table
                              names(HM.782.meanProp.nmg.mo),names(HM.783.meanProp.nmg.mo),
                              names(HM.all.meanProp.nmg.mo),names(NTC.meanProp.nmg.mo), 
                              sprintf("Total mean %s", prop_or_count)))
colnames(prop.mo) <- c("HM.782.meanProp.nmg.mo","HM.783.meanProp.nmg.mo","HM.all.meanProp.nmg.mo","NTC.meanProp.nmg.mo")

for (ro in rownames(prop.mo)) {
  if (ro == sprintf("Total mean %s", prop_or_count)) {
    prop.mo[ro, "HM.782.meanProp.nmg.mo"] <- sum(HM.782.meanProp.nmg.mo)
    prop.mo[ro, "HM.783.meanProp.nmg.mo"] <- sum(HM.783.meanProp.nmg.mo)
    prop.mo[ro, "HM.all.meanProp.nmg.mo"] <- sum(HM.all.meanProp.nmg.mo)
    prop.mo[ro, "NTC.meanProp.nmg.mo"] <- sum(NTC.meanProp.nmg.mo)
  } else {
    for (co in colnames(prop.mo)) {
      if (co == "HM.782.meanProp.nmg.mo") {
        prop.mo[ro, co] <- HM.782.meanProp.nmg.mo[ro]
      } else if (co == "HM.783.meanProp.nmg.mo") {
        prop.mo[ro, co] <- HM.783.meanProp.nmg.mo[ro]
      } else if (co == "HM.all.meanProp.nmg.mo") {
        prop.mo[ro, co] <- HM.all.meanProp.nmg.mo[ro]
      } else if (co == "NTC.meanProp.nmg.mo") {
        prop.mo[ro, co] <- NTC.meanProp.nmg.mo[ro]
      }
    }
  }
}

write.csv(prop.mo, sprintf("/users/tg/jwillis/SLL/Part_1/SLL1_paper/extras/%s.mo.csv", prop_or_count))
# ******************************************* #











unclassified.db <- readRDS("/users/tg/jwillis/SLL/unclassified_identifier_db.rds")



# ****************************************************** #
#### UPDATE unclassified identifier database #### 
# ****************************************************** #
# the numbers in the identifier will be arbitrary, but will remain the same for each unique taxonomy to maintain
# consistency of naming, even between different datasets (OCD, SLL, etc.)
tls <- c("Phylum","Class","Order","Family","Genus","Species")
for (tl in tls) {
  # new table subset to be added to unclassified identifier database
  table.subset <- tax.spec.table[ startsWith(tax.spec.table[ , tl], "unclassified"), 1:(match(tl, tls)+1) ]
  
  # only update table if there are any unclassified taxa at given level
  if ( length(table.subset) > 0 ) {
    
    # combine the existing unclassified identifier database with new table subset
    # but if only 1 row is present in the subset, will appear as a character vector instead of matrix, must force matrix class
    if ( class(table.subset) == "character" ) {
      tmp_table <- rbind( unclassified.db[[ tl ]], t(as.matrix(table.subset)) )
    } else {
      tmp_table <- rbind( unclassified.db[[ tl ]], table.subset )
    }
    
    # reduce the expanded unclassified values at all levels
    tmp_table[ startsWith( tmp_table, "unclassified") ] <- "unclassified"
    
    # keep only those rows which are unique for the columns up to but excluding the current taxa level 
    # (which is how the unclassified values are being identified)
    if ( is.null( rownames(unique(tmp_table[ , 1:match(tl, tls)])) ) ) {
      # this case occurs for the Phylum level because it checks just unique values in Kingdom,
      # so instead of returning another table (as in the else below), 
      # it returns an unnamed character vector, so cannot identify rows
      tmp_table <- unique(tmp_table)
    } else {
      tmp_table <- tmp_table[rownames(unique(tmp_table[ , 1:match(tl, tls)])),]
    }
    
    # update rownames to numbers, 
    # will keep current order so that those that were already present will be labeled with the same number
    rownames(tmp_table) <- 1:nrow(tmp_table)
    # Use numbered rows plus letter of tax level to label the deepest level unclassified taxa
    tletter <- substr(tl, 1, 1)
    tmp_table[ , tl] <- sapply(rownames(tmp_table), function(x) sprintf("unclassified.%s%s",tletter,x))
    rownames(tmp_table) <- tmp_table[ , tl]
    
    # finally, update the table for the given taxonomic level in the unclassified identifier database
    unclassified.db[[ tl ]] <- tmp_table
    
  }
}
# saveRDS(unclassified.db, "/users/tg/jwillis/SLL/unclassified_identifier_db.rds")












# ****************************************************** #
#### UPDATE OTU and TAX TABLES with new unclassified identifiers #### 
# ****************************************************** #

t0 <- Sys.time()

otus <- list()
otus_rel <- list()
taxTables <- list()


# ************************************ #
# function for converting taxa names to appropriate unclassified identifiers
get_unclass_IDs <- function(taxa, level) {
  # get taxa table up to the indicated level
  tmp_tax <- unique(tax.spec.table[ , 1:(match(level, tls)+1) ])
  
  uID <- sapply(taxa, 
                function(x) {
                  if (startsWith(tmp_tax[x, level], "unclassified")) {
                    # get vector of taxa up to one level above current based on the extended
                    taxonomy <- tmp_tax[ x, 1:match(level, tls) ]
                    # remove upper taxa levels from all unclassified names, 
                    # in order to match to the unclassified identifier database in next step
                    taxonomy <- sapply(taxonomy, function(y) strsplit(y, '\\.')[[1]][1])
                    # get correct unclassified identifier by checking where taxonomy matches a row in the table for the given level
                    new.unclass <- rownames( unclassified.db[[ level ]] )[ sapply(rownames(unclassified.db[[ level ]]), 
                                                                                  function(y) sum(unclassified.db[[ level ]][y, 1:match(level, tls)] == taxonomy) == length(taxonomy) ) ]
                    return( new.unclass )
                  } else {
                    return( tmp_tax[x, level] )
                  }
                })
  if (level != "Species" & length( uID[ uID %in% uID[duplicated(uID)] ] ) > 0) {
    # this convoluted bit must exist because in this and other data, there are 2 instances of 
    # families called "Family_XI", but from 2 different orders, causes problems down the line
    # so I must give further identifiers to any instance of taxa names that are duplicated at a given
    # level, but actually have different taxonomy at higher levels
    dups <- uID[uID %in% uID[duplicated(uID)]]
    lev.up.from.dups <- tmp_tax[ names(dups), ncol(tmp_tax)-1 ]
    uID[ uID %in% uID[duplicated(uID)] ] <- sprintf("%s.%s", dups, lev.up.from.dups)
  }
  return( unname(uID) )
  
}

# ************************************ #
# Start by including the Species tax and otu tables, update the rownames with unclassified identifiers
taxTables[["Species"]] <- tax.spec.table
rownames(taxTables[["Species"]]) <- get_unclass_IDs( rownames(taxTables[["Species"]]), "Species")
taxTables[["Species"]] <- unique(taxTables[["Species"]])
otus[["Species"]] <- otu.spec.table
rownames(otus[["Species"]]) <- get_unclass_IDs( rownames(otus[["Species"]]), "Species")
otus_rel[["Species"]] <- otu.spec.table.rel
rownames(otus_rel[["Species"]]) <- get_unclass_IDs( rownames(otus_rel[["Species"]]), "Species")








# ************************************ #
# then get otu_tables at each level by summing the values from same taxa at that level within the species otu table
for (tl in c("Genus","Family","Order","Class","Phylum")) {
  print(tl)
  
  ### update tax table at given level first
  taxTables[[ tl ]] <- unique(tax.spec.table[ , 1:(match(tl, tls)+1) ])
  rownames(taxTables[[ tl ]]) <- get_unclass_IDs( rownames(taxTables[[ tl ]]), tl)
  
  ### get counts of all taxa at given tl from Species counts table that are within each value of the given tl
  taxa <- taxTables[[ tl ]][ , tl]
  otus[[ tl ]] <- apply(otus[["Species"]], 2,
                        function(x) sapply(taxa,
                                           function(y) sum(x[ names( taxTables[["Species"]][,tl][ y == taxTables[["Species"]][,tl] ] )]
                                           )
                        )
  )
  # example just to verify by hand that values are correct (they are now:
  # otus[["Genus"]][ 101:105, 1:5 ]
  # sum(otus[["Species"]][,4][ names( taxTables[["Species"]][,"Genus"][ taxTables[["Species"]][,"Genus"]==taxTables[["Genus"]][,"Genus"]["F0058"] ] )])
  # otus[["Species"]][names( taxTables[["Species"]][,"Genus"][ taxTables[["Species"]][,"Genus"]==taxTables[["Genus"]][,"Genus"]["F0058"] ] ), 4]
  
  # get normalized values -- relative abundances
  otus_rel[[ tl ]] <- apply( otus[[ tl ]], 2, function(x) 100 * x/sum(x))
  
}









# ************************************ #
# Make species names proper scientific names (by making them <Genus species>)

rownames(otus[["Species"]]) <- unname(sapply(rownames(otus[["Species"]]), 
                                             function(x) {
                                               if (startsWith(x, "unclassified")) {
                                                 if (startsWith(taxTables[["Species"]][x,"Genus"], "unclassified")) {
                                                   paste("unclassified", x, sep = ' ')
                                                 } else {
                                                   paste(taxTables[["Species"]][x,"Genus"], x, sep = ' ')
                                                 }
                                               } else {
                                                 x
                                               }
                                               
                                             }))

rownames(otus_rel[["Species"]]) <- unname(sapply(rownames(otus_rel[["Species"]]), 
                                                 function(x) {
                                                   if (startsWith(x, "unclassified")) {
                                                     if (startsWith(taxTables[["Species"]][x,"Genus"], "unclassified")) {
                                                       paste("unclassified", x, sep = ' ')
                                                     } else {
                                                       paste(taxTables[["Species"]][x,"Genus"], x, sep = ' ')
                                                     }
                                                   } else {
                                                     x
                                                   }
                                                   
                                                 }))

rownames(taxTables[["Species"]]) <- unname(sapply(rownames(taxTables[["Species"]]), 
                                                  function(x) {
                                                    if (startsWith(x, "unclassified")) {
                                                      if (startsWith(taxTables[["Species"]][x,"Genus"], "unclassified")) {
                                                        paste("unclassified", x, sep = ' ')
                                                      } else {
                                                        paste(taxTables[["Species"]][x,"Genus"], x, sep = ' ')
                                                      }
                                                    } else {
                                                      x
                                                    }
                                                    
                                                  }))

taxTables[["Species"]][ , "Species"] <- rownames(taxTables[["Species"]])



# save tax and otu table objects
saveRDS(taxTables, "/users/tg/jwillis/SLL/Part_1/DADA2/taxTables.rds")
saveRDS(otus, "/users/tg/jwillis/SLL/Part_1/DADA2/otus.rds")
saveRDS(otus_rel, "/users/tg/jwillis/SLL/Part_1/DADA2/otus_rel.rds")

print(Sys.time() - t0)
# ************************************ #






# save.image(file = "~/Downloads/dada2_SLL1.RData")









# ****************************************************** #
#### CREATE PHYLOGENETIC TREE #### 
# ****************************************************** #

# From: https://f1000research.com/articles/5-1492/v2
# "Phylogenetic relatedness is commonly used to inform downstream analyses, especially the 
#  calculation of phylogeny-aware distances between microbial communities. The DADA2 sequence 
#  inference method is reference-free, so we must construct the phylogenetic tree relating the 
#  inferred sequence variants de novo. We begin by performing a multiple-alignment using the 
#  DECIPHER R package11."
# library(DECIPHER)
# 
# seqs <- getSequences(seqtab.nochim)
# names(seqs) <- seqs # This propagates to the tip labels of the tree
# alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
# 
# # Determining distance matrix based on shared 7-mers:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 3811.71 secs
# # 
# # Clustering into groups by similarity:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 4527.83 secs
# # 
# # Aligning Sequences:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 678.05 secs
# # 
# # Iteration 1 of 2:
# #   
# #   Determining distance matrix based on alignment:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 574.37 secs
# # 
# # Reclustering into groups by similarity:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 4375.18 secs
# # 
# # Realigning Sequences:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 183.79 secs
# # 
# # Iteration 2 of 2:
# #   
# #   Determining distance matrix based on alignment:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 561.08 secs
# # 
# # Reclustering into groups by similarity:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 4446.97 secs
# # 
# # Realigning Sequences:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 57.36 secs
# # 
# # Refining the alignment:
# #   |========================================================================================================================| 100%
# # 
# # Time difference of 126.27 secs



# randomly sampling 1 sequence from all sequences assigned to each species
#   so will end with vector the length of number of species
#   in theory, the sequences assigned to a given species should be similar enough
#     that they would align together in the total alignment (19626 sequences here)
#     but can now do this much more efficiently (only 598 now)
non.euk <- names(sapply(unique(tax.fixNames.spec[,"Species"]), 
                        function(x) x %in% rownames(otu.spec.table))[ sapply(unique(tax.fixNames.spec[,"Species"]), 
                                                                             function(x) x %in% rownames(otu.spec.table))])
species.seqs <- sapply(sort(non.euk), 
                       function(x) sample(rownames(as.data.frame(tax.fixNames.spec)[tax.fixNames.spec[,"Species"]==x,]),1))
# get proper ID for unclassified species
names(species.seqs) <- unname(unlist(get_unclass_IDs(names(species.seqs), "Species")))

# make proper scientific names by adding the Genus
# spec.only is a vector of species names, for which the name attributes are the full scientific names
spec.only <- sapply(rownames(taxTables[["Species"]]), function(x) strsplit(x,' ')[[1]][2])
names(species.seqs) <- unname(sapply(names(species.seqs), function(x) 
  if (startsWith(x, "unclassified")) {
    names(spec.only[spec.only==x])
  } else {
    x
  }))




library(DECIPHER)
alignment <- AlignSeqs(DNAStringSet(species.seqs), anchor=NA)



# "The phangorn R package is then used to construct a phylogenetic tree. 
#  Here we first construct a neighbor-joining tree, and then fit a GTR+G+I 
#  (Generalized time-reversible with Gamma rate variation) maximum likelihood 
#  tree using the neighbor-joining tree as a starting point."
library(phangorn)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

saveRDS(fitGTR, file = "/users/tg/jwillis/SLL/Part_1/R_objects/SLL1.fitGTR.rds")


















# ****************************************************** #
#### HANDOFF TO PHYLOSEQ #### 
# ****************************************************** #
library(phyloseq)

# make sample_data table
samples.out <- rownames(seqtab.nochim)

data.que <- read.delim(sprintf("%s/Part_1/Surveys/SLL_survey_data_repaired.csv", home_dir), header=T)

# ignore samples of duplicated identifiers bc cannot determine which is correct or what the other should be changed to
dubs <- as.character(data.que[duplicated(data.que[,"Qsamples"]),"Qsamples"])
x <- 1
both_dubs <- c()
for (i in 1:length(rownames(data.que))) {
  if (data.que[i,"Qsamples"] %in% dubs) {
    both_dubs[x] <- as.numeric(rownames(data.que)[i])
    x <- x+1
  }
}
data.que <- data.que[-both_dubs,]
# rownames(data.que) <- gsub('ll','LL', gsub("-","",as.vector(data.que[,"Qsamples"])))
rownames(data.que) <- gsub('ll','LL', as.vector(data.que[,"Qsamples"]))

#collect list of questions and remove question from top of data.que
full_questions <- data.que[1,]
data.que <- data.que[-1,]
data.que <- apply(data.que, 2, factor)
data.que <- as.data.frame(data.que) #because apply converts it to a matrix
colnames(data.que)[102] <- "Socioeconomic"
data.que <- data.que[ rownames(data.que) %in% gsub('-', '', samples.out), ]
rownames(data.que) <- sapply(rownames(data.que), function(x) sprintf("%s-%s-%s", 
                                                                     substr(x, 1, 3),
                                                                     substr(x, 4, 7),
                                                                     substr(x, 8, 11)))

# construct phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(data.que), 
                   tax_table(tax_spec))
# ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

# visualize a-diversity
plot_richness(ps, x="Q2", measures=c("Shannon", "Simpson"), color="Q27") + theme_bw()


# ordinate
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="Q2", title="Bray NMDS")


# bar plot:
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU) * 100)
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Q2", fill="Genus") + facet_wrap(~Q15, scales="free_x")

# These bar plots dont show anything vastly different between Early and Late,
#   plenty more that can be check with phyloseq still...










cursam <- unique(sapply(strsplit(basename(list.files("/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads/")), "_"), `[`, 1))
missing.samps <- sample.names2[!sample.names2 %in% cursam]

missing.files <- fnFs[sapply(strsplit(basename(fnFs), "_"), `[`, 3) %in% missing.samps]
missing.files.r <- fnRs[sapply(strsplit(basename(fnRs), "_"), `[`, 3) %in% missing.samps]

plotQualityProfile(missing.files[1:4])
plotQualityProfile(missing.files.r[1:4])







