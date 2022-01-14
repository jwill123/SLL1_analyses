



# *********** #

# for the pre-assembled F and R reads using Pear v0.9.6
# ran from the folder /users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2
#  $ pear -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -o joined/HM-782D-plate1_PEAR
#  $ gzip joined/HM-782D-plate1_PEAR.assembled.fastq
filt.mock.PEAR <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_PEAR.assembled.fastq.gz"
plotQualityProfile(filt.mock.PEAR) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks not excellent in the middle, but will move forward
filt.mock4.PEAR <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate4_PEAR.assembled.fastq.gz"
plotQualityProfile(filt.mock4.PEAR) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks not excellent in the middle, but will move forward
filt.mock6.PEAR <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-783D-plate6_PEAR.assembled.fastq.gz"
plotQualityProfile(filt.mock6.PEAR) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks not excellent in the middle, but will move forward


# *********** #

# for the pre-assembled F and R reads using pandaseq v2.11 with -A simple_bayesian (option for algorithm)
#  $ pandaseq -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -F -A simple_bayesian \
#      -w joined/HM-782D-plate1_pandaseq.simple_bayesian.fastq 2> joined/HM-782D-plate1_pandaseq.simple_bayesian.fastq.log
#  $ gzip joined/HM-782D-plate1_pandaseq.simple_bayesian.fastq
filt.mock.panda.sb <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_pandaseq.simple_bayesian.fastq.gz"
plotQualityProfile(filt.mock.panda.sb) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks terrible in the middle


# for the pre-assembled F and R reads using pandaseq v2.11 with -A ea_util (option for algorithm) (fastqJoin)
#  $ pandaseq -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -F -A ea_util \
#      -w joined/HM-782D-plate1_pandaseq.ea_util.fastq 2> joined/HM-782D-plate1_pandaseq.ea_util.fastq.log
#  $ gzip joined/HM-782D-plate1_pandaseq.ea_util.fastq
filt.mock.panda.ea <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_pandaseq.ea_util.fastq.gz"
plotQualityProfile(filt.mock.panda.ea) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks fine in the middle, can proceed



# for the pre-assembled F and R reads using pandaseq v2.11 with -A flash (option for algorithm)
#  $ pandaseq -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -F -A flash \
#      -w joined/HM-782D-plate1_pandaseq.flash.fastq 2> joined/HM-782D-plate1_pandaseq.flash.fastq.log
#  $ gzip joined/HM-782D-plate1_pandaseq.flash.fastq
filt.mock.panda.flash <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_pandaseq.flash.fastq.gz"
plotQualityProfile(filt.mock.panda.flash) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks terrible in the middle


# for the pre-assembled F and R reads using pandaseq v2.11 with -A pear (option for algorithm)
#  $ pandaseq -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -F -A pear \
#      -w joined/HM-782D-plate1_pandaseq.pear.fastq 2> joined/HM-782D-plate1_pandaseq.pear.fastq.log
#  $ gzip joined/HM-782D-plate1_pandaseq.pear.fastq
filt.mock.panda.pear <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_pandaseq.pear.fastq.gz"
plotQualityProfile(filt.mock.panda.pear) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks terrible in the middle...no chance these will be good



# for the pre-assembled F and R reads using pandaseq v2.11 with -A rdp_mle (option for algorithm)
#  $ pandaseq -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -F -A rdp_mle \
#      -w joined/HM-782D-plate1_pandaseq.rdp_mle.fastq 2> joined/HM-782D-plate1_pandaseq.rdp_mle.fastq.log
#  $ gzip joined/HM-782D-plate1_pandaseq.rdp_mle.fastq
filt.mock.panda.rdp_mle <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_pandaseq.rdp_mle.fastq.gz"
plotQualityProfile(filt.mock.panda.rdp_mle) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks fine in the middle, can move forward



# for the pre-assembled F and R reads using pandaseq v2.11 with -A stitch (option for algorithm)
#  $ pandaseq -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -F -A stitch \
#      -w joined/HM-782D-plate1_pandaseq.stitch.fastq 2> joined/HM-782D-plate1_pandaseq.stitch.fastq.log
#  $ gzip joined/HM-782D-plate1_pandaseq.stitch.fastq
filt.mock.panda.stitch <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_pandaseq.stitch.fastq.gz"
plotQualityProfile(filt.mock.panda.stitch) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks pretty bad in the middle, will move forward anyway



# for the pre-assembled F and R reads using pandaseq v2.11 with -A uparse (option for algorithm)
#  $ pandaseq -f HM-782D-plate1_F_filt.fastq.gz -r HM-782D-plate1_R_filt.fastq.gz -F -A uparse \
#      -w joined/HM-782D-plate1_pandaseq.uparse.fastq 2> joined/HM-782D-plate1_pandaseq.uparse.fastq.log
#  $ gzip joined/HM-782D-plate1_pandaseq.uparse.fastq
filt.mock.panda.uparse <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_pandaseq.uparse.fastq.gz"
plotQualityProfile(filt.mock.panda.uparse) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 40, 5))
# looks weird in the middle, will move forward anyway

# *********** #


# for the pre-assembled F and R reads using fastq-join v1.3.1
#  $ fastq-join HM-782D-plate1_F_filt.fastq.gz HM-782D-plate1_R_filt.fastq.gz -o joined/HM-782D-plate1_fqj.
#  $ mv joined/HM-782D-plate1_fqj.join joined/HM-782D-plate1_fqj.join.fastq
#  $ gzip joined/HM-782D-plate1_fqj.join.fastq
filt.mock.fqj <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate1_fqj.join.fastq.gz"
plotQualityProfile(filt.mock.fqj) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 45, 5))
# looks good, can move forward
filt.mock4.fqj <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-782D-plate4_fqj.join.fastq.gz"
plotQualityProfile(filt.mock4.fqj) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 45, 5))
# looks good, can move forward
filt.mock6.fqj <- "/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/filtered_reads_2/joined/HM-783D-plate6_fqj.join.fastq.gz"
plotQualityProfile(filt.mock6.fqj) + scale_x_continuous(breaks = seq(0, 500, 25)) + scale_y_continuous(breaks = seq(0, 45, 5))
# looks alright, can move forward

# *********** #






# ****************************************************** #
#### LEARN THE ERROR RATES #### 
# ****************************************************** #

# *********** #
err.mock.PEAR <- learnErrors(filt.mock.PEAR, multithread=TRUE)
plotErrors(err.mock.PEAR, nominalQ=TRUE) # all looks good
err.mock4.PEAR <- learnErrors(filt.mock4.PEAR, multithread=TRUE)
plotErrors(err.mock4.PEAR, nominalQ=TRUE) # all looks good
err.mock6.PEAR <- learnErrors(filt.mock6.PEAR, multithread=TRUE)
plotErrors(err.mock6.PEAR, nominalQ=TRUE) # all looks good

# *********** #

err.mock.panda.sb     <- learnErrors(filt.mock.panda.sb, multithread=TRUE)
err.mock.panda.ea     <- learnErrors(filt.mock.panda.ea, multithread=TRUE)
err.mock.panda.flash  <- learnErrors(filt.mock.panda.flash, multithread=TRUE)
err.mock.panda.pear   <- learnErrors(filt.mock.panda.pear, multithread=TRUE)
err.mock.panda.rdp    <- learnErrors(filt.mock.panda.rdp_mle, multithread=TRUE)
err.mock.panda.stitch <- learnErrors(filt.mock.panda.stitch, multithread=TRUE)
err.mock.panda.uparse <- learnErrors(filt.mock.panda.uparse, multithread=TRUE)

plotErrors(err.mock.panda.sb, nominalQ = TRUE) # mostly good
plotErrors(err.mock.panda.ea, nominalQ = TRUE) # most are fine
plotErrors(err.mock.panda.flash, nominalQ = TRUE) # not good
plotErrors(err.mock.panda.pear, nominalQ = TRUE) # some not great
plotErrors(err.mock.panda.rdp, nominalQ = TRUE) # all looks good
plotErrors(err.mock.panda.stitch, nominalQ = TRUE) # mostly good
plotErrors(err.mock.panda.uparse, nominalQ = TRUE) # not bad

# *********** #

err.mock.fqj <- learnErrors(filt.mock.fqj, multithread=TRUE)
plotErrors(err.mock.fqj, nominalQ = TRUE) # all looks good
err.mock4.fqj <- learnErrors(filt.mock4.fqj, multithread=TRUE)
plotErrors(err.mock4.fqj, nominalQ = TRUE) # all looks good
err.mock6.fqj <- learnErrors(filt.mock6.fqj, multithread=TRUE)
plotErrors(err.mock6.fqj, nominalQ = TRUE) # all looks good

# *********** #



# ****************************************************** #
#### DEREPLICATION #### 
# ****************************************************** #

derep.mock.PEAR  <- derepFastq(filt.mock.PEAR, verbose=TRUE)
derep.mock4.PEAR <- derepFastq(filt.mock4.PEAR, verbose=TRUE)
derep.mock6.PEAR <- derepFastq(filt.mock6.PEAR, verbose=TRUE)

# *********** #

derep.mock.panda.sb     <- derepFastq(filt.mock.panda.sb, verbose=TRUE)
derep.mock.panda.ea     <- derepFastq(filt.mock.panda.ea, verbose=TRUE)
derep.mock.panda.flash  <- derepFastq(filt.mock.panda.flash, verbose=TRUE)
derep.mock.panda.pear   <- derepFastq(filt.mock.panda.pear, verbose=TRUE)
derep.mock.panda.rdp    <- derepFastq(filt.mock.panda.rdp_mle, verbose=TRUE)
derep.mock.panda.stitch <- derepFastq(filt.mock.panda.stitch, verbose=TRUE)
derep.mock.panda.uparse <- derepFastq(filt.mock.panda.uparse, verbose=TRUE)

# *********** #

derep.mock.fqj  <- derepFastq(filt.mock.fqj, verbose=TRUE)
derep.mock4.fqj <- derepFastq(filt.mock4.fqj, verbose=TRUE)
derep.mock6.fqj <- derepFastq(filt.mock6.fqj, verbose=TRUE)

# *********** #




# ****************************************************** #
#### SAMPLE INFERENCE #### 
# ****************************************************** #

dada.mock.PEAR  <- dada(derep.mock.PEAR, err=err.mock.PEAR, multithread=TRUE)
dada.mock4.PEAR <- dada(derep.mock4.PEAR, err=err.mock4.PEAR, multithread=TRUE)
dada.mock6.PEAR <- dada(derep.mock6.PEAR, err=err.mock6.PEAR, multithread=TRUE)

# *********** #

dada.mock.panda.sb     <- dada(derep.mock.panda.sb, err=err.mock.panda.sb, multithread=TRUE)
dada.mock.panda.ea     <- dada(derep.mock.panda.ea, err=err.mock.panda.ea, multithread=TRUE)
dada.mock.panda.flash  <- dada(derep.mock.panda.flash, err=err.mock.panda.flash, multithread=TRUE)
dada.mock.panda.pear   <- dada(derep.mock.panda.pear, err=err.mock.panda.pear, multithread=TRUE)
dada.mock.panda.rdp    <- dada(derep.mock.panda.rdp, err=err.mock.panda.rdp, multithread=TRUE)
dada.mock.panda.stitch <- dada(derep.mock.panda.stitch, err=err.mock.panda.stitch, multithread=TRUE)
dada.mock.panda.uparse <- dada(derep.mock.panda.uparse, err=err.mock.panda.uparse, multithread=TRUE)

# *********** #

dada.mock.fqj  <- dada(derep.mock.fqj, err=err.mock.fqj, multithread=TRUE)
dada.mock4.fqj <- dada(derep.mock4.fqj, err=err.mock4.fqj, multithread=TRUE)
dada.mock6.fqj <- dada(derep.mock6.fqj, err=err.mock6.fqj, multithread=TRUE)

# *********** #
# *********** #

dada2:::checkConvergence(err.mock.PEAR)
dada2:::checkConvergence(err.mock4.PEAR)
dada2:::checkConvergence(err.mock6.PEAR)

# *********** #

dada2:::checkConvergence(err.mock.panda.sb)
dada2:::checkConvergence(err.mock.panda.ea)
dada2:::checkConvergence(err.mock.panda.flash)
dada2:::checkConvergence(err.mock.panda.pear)
dada2:::checkConvergence(err.mock.panda.rdp)
dada2:::checkConvergence(err.mock.panda.stitch)
dada2:::checkConvergence(err.mock.panda.uparse)

# *********** #

dada2:::checkConvergence(err.mock.fqj)
dada2:::checkConvergence(err.mock4.fqj)
dada2:::checkConvergence(err.mock6.fqj)

# *********** #




# ****************************************************** #
#### CONSTRUCT SEQUENCE TABLE #### 
# ****************************************************** #

seqtab.mock.PEAR  <- makeSequenceTable(dada.mock.PEAR)
seqtab.mock4.PEAR <- makeSequenceTable(dada.mock4.PEAR)
seqtab.mock6.PEAR <- makeSequenceTable(dada.mock6.PEAR)

table(nchar(getSequences(seqtab.mock.PEAR)))
table(nchar(getSequences(seqtab.mock4.PEAR)))
table(nchar(getSequences(seqtab.mock6.PEAR)))

# *********** #

seqtab.mock.panda.sb     <- makeSequenceTable(dada.mock.panda.sb)
seqtab.mock.panda.ea     <- makeSequenceTable(dada.mock.panda.ea)
seqtab.mock.panda.flash  <- makeSequenceTable(dada.mock.panda.flash)
seqtab.mock.panda.pear   <- makeSequenceTable(dada.mock.panda.pear)
seqtab.mock.panda.rdp    <- makeSequenceTable(dada.mock.panda.rdp)
seqtab.mock.panda.stitch <- makeSequenceTable(dada.mock.panda.stitch)
seqtab.mock.panda.uparse <- makeSequenceTable(dada.mock.panda.uparse)

table(nchar(getSequences(seqtab.mock.panda.sb)))
table(nchar(getSequences(seqtab.mock.panda.ea)))
table(nchar(getSequences(seqtab.mock.panda.flash))) # only 3 sequences...not good 
table(nchar(getSequences(seqtab.mock.panda.pear)))
table(nchar(getSequences(seqtab.mock.panda.rdp)))
table(nchar(getSequences(seqtab.mock.panda.stitch)))
table(nchar(getSequences(seqtab.mock.panda.uparse)))

# *********** #

seqtab.mock.fqj  <- makeSequenceTable(dada.mock.fqj)
seqtab.mock4.fqj <- makeSequenceTable(dada.mock4.fqj)
seqtab.mock6.fqj <- makeSequenceTable(dada.mock6.fqj)

table(nchar(getSequences(seqtab.mock.fqj)))
table(nchar(getSequences(seqtab.mock4.fqj)))
table(nchar(getSequences(seqtab.mock6.fqj)))

# *********** #




# ****************************************************** #
#### REMOVE CHIMERAS #### 
# ****************************************************** #

mock.ref <- getSequences(file.path(path.tut, "HMP_MOCK.v35.fasta"))
mock.genera <- c("Acinetobacter","Actinomyces","Bacillus","Bacteroides","Clostridium","Deinococcus","Enterococcus","Escherichia","Helicobacter",
                 "Lactobacillus","Listeria","Neisseria","Propionibacterium","Pseudomonas","Rhodobacter","Staphylococcus","Staphylococcus",
                 "Streptococcus","Streptococcus","Streptococcus")
mock.genera.copies <- c(10000,1000,100000,1000,100000,1000,1000,1000000,10000,10000,10000,10000,10000,100000,1000000,100000,1000000,
                        100000,1000000,1000)
mock.genera.copies <- sapply(mock.genera.copies, function(x) round(100 * x/sum(mock.genera.copies), 2))
mock.genera.even <- rep(10000, 20)
mock.genera.even <- sapply(mock.genera.even, function(x) round(100 * x/sum(mock.genera.even), 2))
mock.genera.copies <- as.data.frame(cbind(mock.genera, mock.genera.even, mock.genera.copies))
colnames(mock.genera.copies) = c("Genus","Even (%)","Staggered (%)")


# *********** #

seqtab.nochim.mock.PEAR <- removeBimeraDenovo(seqtab.mock.PEAR, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.PEAR)
dim(seqtab.nochim.mock.PEAR)
sum(seqtab.nochim.mock.PEAR)/sum(seqtab.mock.PEAR)
unqs.mock.PEAR <- seqtab.nochim.mock.PEAR[1,]
unqs.mock.PEAR <- sort(unqs.mock.PEAR[unqs.mock.PEAR > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.PEAR), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.PEAR), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")

seqtab.nochim.mock4.PEAR <- removeBimeraDenovo(seqtab.mock4.PEAR, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock4.PEAR)
dim(seqtab.nochim.mock4.PEAR)
sum(seqtab.nochim.mock4.PEAR)/sum(seqtab.mock4.PEAR)
unqs.mock4.PEAR <- seqtab.nochim.mock4.PEAR[1,]
unqs.mock4.PEAR <- sort(unqs.mock4.PEAR[unqs.mock4.PEAR > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock4.PEAR), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock4.PEAR), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")

seqtab.nochim.mock6.PEAR <- removeBimeraDenovo(seqtab.mock6.PEAR, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock6.PEAR)
dim(seqtab.nochim.mock6.PEAR)
sum(seqtab.nochim.mock6.PEAR)/sum(seqtab.mock6.PEAR)
unqs.mock6.PEAR <- seqtab.nochim.mock6.PEAR[1,]
unqs.mock6.PEAR <- sort(unqs.mock6.PEAR[unqs.mock6.PEAR > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock6.PEAR), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock6.PEAR), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")


# *********** #

seqtab.nochim.mock.panda.sb <- removeBimeraDenovo(seqtab.mock.panda.sb, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.panda.sb)
dim(seqtab.nochim.mock.panda.sb)
sum(seqtab.nochim.mock.panda.sb)/sum(seqtab.mock.panda.sb)
unqs.mock.panda.sb <- seqtab.nochim.mock.panda.sb[1,]
unqs.mock.panda.sb <- sort(unqs.mock.panda.sb[unqs.mock.panda.sb > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.panda.sb), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.panda.sb), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")


seqtab.nochim.mock.panda.ea <- removeBimeraDenovo(seqtab.mock.panda.ea, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.panda.ea)
dim(seqtab.nochim.mock.panda.ea)
sum(seqtab.nochim.mock.panda.ea)/sum(seqtab.mock.panda.ea)
unqs.mock.panda.ea <- seqtab.nochim.mock.panda.ea[1,]
unqs.mock.panda.ea <- sort(unqs.mock.panda.ea[unqs.mock.panda.ea > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.panda.ea), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.panda.ea), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")


seqtab.nochim.mock.panda.flash <- removeBimeraDenovo(seqtab.mock.panda.flash, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.panda.flash)
dim(seqtab.nochim.mock.panda.flash)
sum(seqtab.nochim.mock.panda.flash)/sum(seqtab.mock.panda.flash)
unqs.mock.panda.flash <- seqtab.nochim.mock.panda.flash[1,]
unqs.mock.panda.flash <- sort(unqs.mock.panda.flash[unqs.mock.panda.flash > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.panda.flash), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.panda.flash), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")



seqtab.nochim.mock.panda.pear <- removeBimeraDenovo(seqtab.mock.panda.pear, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.panda.pear)
dim(seqtab.nochim.mock.panda.pear)
sum(seqtab.nochim.mock.panda.pear)/sum(seqtab.mock.panda.pear)
unqs.mock.panda.pear <- seqtab.nochim.mock.panda.pear[1,]
unqs.mock.panda.pear <- sort(unqs.mock.panda.pear[unqs.mock.panda.pear > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.panda.pear), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.panda.pear), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")



seqtab.nochim.mock.panda.rdp <- removeBimeraDenovo(seqtab.mock.panda.rdp, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.panda.rdp)
dim(seqtab.nochim.mock.panda.rdp)
sum(seqtab.nochim.mock.panda.rdp)/sum(seqtab.mock.panda.rdp)
unqs.mock.panda.rdp <- seqtab.nochim.mock.panda.rdp[1,]
unqs.mock.panda.rdp <- sort(unqs.mock.panda.rdp[unqs.mock.panda.rdp > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.panda.rdp), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.panda.rdp), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")


seqtab.nochim.mock.panda.stitch <- removeBimeraDenovo(seqtab.mock.panda.stitch, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.panda.stitch)
dim(seqtab.nochim.mock.panda.stitch)
sum(seqtab.nochim.mock.panda.stitch)/sum(seqtab.mock.panda.stitch)
unqs.mock.panda.stitch <- seqtab.nochim.mock.panda.stitch[1,]
unqs.mock.panda.stitch <- sort(unqs.mock.panda.stitch[unqs.mock.panda.stitch > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.panda.stitch), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.panda.stitch), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")


seqtab.nochim.mock.panda.uparse <- removeBimeraDenovo(seqtab.mock.panda.uparse, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.panda.uparse)
dim(seqtab.nochim.mock.panda.uparse)
sum(seqtab.nochim.mock.panda.uparse)/sum(seqtab.mock.panda.uparse)
unqs.mock.panda.uparse <- seqtab.nochim.mock.panda.uparse[1,]
unqs.mock.panda.uparse <- sort(unqs.mock.panda.uparse[unqs.mock.panda.uparse > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.panda.uparse), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.panda.uparse), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")


# *********** #

seqtab.nochim.mock.fqj <- removeBimeraDenovo(seqtab.mock.fqj, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock.fqj)
dim(seqtab.nochim.mock.fqj)
sum(seqtab.nochim.mock.fqj)/sum(seqtab.mock.fqj)
unqs.mock.fqj <- seqtab.nochim.mock.fqj[1,]
unqs.mock.fqj <- sort(unqs.mock.fqj[unqs.mock.fqj > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock.fqj), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock.fqj), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")

seqtab.nochim.mock4.fqj <- removeBimeraDenovo(seqtab.mock4.fqj, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock4.fqj)
dim(seqtab.nochim.mock4.fqj)
sum(seqtab.nochim.mock4.fqj)/sum(seqtab.mock4.fqj)
unqs.mock4.fqj <- seqtab.nochim.mock4.fqj[1,]
unqs.mock4.fqj <- sort(unqs.mock4.fqj[unqs.mock4.fqj > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock4.fqj), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock4.fqj), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")

seqtab.nochim.mock6.fqj <- removeBimeraDenovo(seqtab.mock6.fqj, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.mock6.fqj)
dim(seqtab.nochim.mock6.fqj)
sum(seqtab.nochim.mock6.fqj)/sum(seqtab.mock6.fqj)
unqs.mock6.fqj <- seqtab.nochim.mock6.fqj[1,]
unqs.mock6.fqj <- sort(unqs.mock6.fqj[unqs.mock6.fqj > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock6.fqj), "sample sequences present in the Mock community.\n")
match.ref <- sum(sapply(names(unqs.mock6.fqj), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")

# *********** #


# *********** #



sn.mock.PEAR <- seqtab.nochim.mock.PEAR
sn.mock.PEAR <- sn.mock.PEAR[ , colSums(sn.mock.PEAR) > 0]
tax.mock.PEAR <- assignTaxonomy(sn.mock.PEAR, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)

sn.mock4.PEAR <- seqtab.nochim.mock4.PEAR
sn.mock4.PEAR <- sn.mock4.PEAR[ , colSums(sn.mock4.PEAR) > 0]
tax.mock4.PEAR <- assignTaxonomy(sn.mock4.PEAR, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)

sn.mock6.PEAR <- seqtab.nochim.mock6.PEAR
sn.mock6.PEAR <- sn.mock6.PEAR[ , colSums(sn.mock6.PEAR) > 0]
tax.mock6.PEAR <- assignTaxonomy(sn.mock6.PEAR, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)

sn.mock.fqj <- seqtab.nochim.mock.fqj
sn.mock.fqj <- sn.mock.fqj[ , colSums(sn.mock.fqj) > 0]
tax.mock.fqj <- assignTaxonomy(sn.mock.fqj, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)

sn.mock4.fqj <- seqtab.nochim.mock4.fqj
sn.mock4.fqj <- sn.mock4.fqj[ , colSums(sn.mock4.fqj) > 0]
tax.mock4.fqj <- assignTaxonomy(sn.mock4.fqj, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)

sn.mock6.fqj <- seqtab.nochim.mock6.fqj
sn.mock6.fqj <- sn.mock6.fqj[ , colSums(sn.mock6.fqj) > 0]
tax.mock6.fqj <- assignTaxonomy(sn.mock6.fqj, sprintf("%s/silva_nr_v128_train_set.fa.gz",path.tut), multithread=TRUE)

# *********** #



# *********** #

# mock.PEAR
# change NAs to unclassified
tax.mock.PEAR.fixNames <- tax.mock.PEAR
tax.mock.PEAR.fixNames[is.na(tax.mock.PEAR.fixNames)] <- "unclassified"
tax.mock.PEAR.fixNames <- t(apply(tax.mock.PEAR.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.mock.PEAR.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.mock.PEAR.fixNames) <- colnames(tax.mock.PEAR)
# get otu table
otu.table.mock.PEAR <- seqtab.nochim.mock.PEAR
colnames(otu.table.mock.PEAR) <- unname(tax.mock.PEAR.fixNames[colnames(seqtab.nochim.mock.PEAR) , "Genus"])
non.bact.mock.PEAR <- tax.mock.PEAR.fixNames[ tax.mock.PEAR.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
otu.table.mock.PEAR <- t(as.matrix(otu.table.mock.PEAR[ , ! colnames(otu.table.mock.PEAR) %in% non.bact.mock.PEAR]))
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table.mock.PEAR <- t(otu.table.mock.PEAR %*% sapply(unique(colnames(otu.table.mock.PEAR)),"==", colnames(otu.table.mock.PEAR)))
otu.table.mock.PEAR.rel <- apply(otu.table.mock.PEAR, 2, function(x) round(100 * x/sum(x), 2))
otu.table.mock.PEAR.both <- cbind(otu.table.mock.PEAR, otu.table.mock.PEAR.rel)
colnames(otu.table.mock.PEAR.both) <- c("Counts", "Rel_abund")
otu.table.mock.PEAR.both <- otu.table.mock.PEAR.both[ sort(rownames(otu.table.mock.PEAR.both)), ]


# mock4.PEAR
# change NAs to unclassified
tax.mock4.PEAR.fixNames <- tax.mock4.PEAR
tax.mock4.PEAR.fixNames[is.na(tax.mock4.PEAR.fixNames)] <- "unclassified"
tax.mock4.PEAR.fixNames <- t(apply(tax.mock4.PEAR.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.mock4.PEAR.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                          paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                          all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.mock4.PEAR.fixNames) <- colnames(tax.mock4.PEAR)
# get otu table
otu.table.mock4.PEAR <- seqtab.nochim.mock4.PEAR
colnames(otu.table.mock4.PEAR) <- unname(tax.mock4.PEAR.fixNames[colnames(seqtab.nochim.mock4.PEAR) , "Genus"])
non.bact.mock4.PEAR <- tax.mock4.PEAR.fixNames[ tax.mock4.PEAR.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
otu.table.mock4.PEAR <- t(as.matrix(otu.table.mock4.PEAR[ , ! colnames(otu.table.mock4.PEAR) %in% non.bact.mock4.PEAR]))
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table.mock4.PEAR <- t(otu.table.mock4.PEAR %*% sapply(unique(colnames(otu.table.mock4.PEAR)),"==", colnames(otu.table.mock4.PEAR)))
otu.table.mock4.PEAR.rel <- apply(otu.table.mock4.PEAR, 2, function(x) round(100 * x/sum(x), 2))
otu.table.mock4.PEAR.both <- cbind(otu.table.mock4.PEAR, otu.table.mock4.PEAR.rel)
colnames(otu.table.mock4.PEAR.both) <- c("Counts", "Rel_abund")
otu.table.mock4.PEAR.both <- otu.table.mock4.PEAR.both[ sort(rownames(otu.table.mock4.PEAR.both)), ]







# mock6.PEAR
# change NAs to unclassified
tax.mock6.PEAR.fixNames <- tax.mock6.PEAR
tax.mock6.PEAR.fixNames[is.na(tax.mock6.PEAR.fixNames)] <- "unclassified"
tax.mock6.PEAR.fixNames <- t(apply(tax.mock6.PEAR.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.mock6.PEAR.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                           paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                           all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.mock6.PEAR.fixNames) <- colnames(tax.mock6.PEAR)
# get otu table
otu.table.mock6.PEAR <- seqtab.nochim.mock6.PEAR
colnames(otu.table.mock6.PEAR) <- unname(tax.mock6.PEAR.fixNames[colnames(seqtab.nochim.mock6.PEAR) , "Genus"])
non.bact.mock6.PEAR <- tax.mock6.PEAR.fixNames[ tax.mock6.PEAR.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
otu.table.mock6.PEAR <- t(as.matrix(otu.table.mock6.PEAR[ , ! colnames(otu.table.mock6.PEAR) %in% non.bact.mock6.PEAR]))
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table.mock6.PEAR <- t(otu.table.mock6.PEAR %*% sapply(unique(colnames(otu.table.mock6.PEAR)),"==", colnames(otu.table.mock6.PEAR)))
otu.table.mock6.PEAR.rel <- apply(otu.table.mock6.PEAR, 2, function(x) round(100 * x/sum(x), 2))
otu.table.mock6.PEAR.both <- cbind(otu.table.mock6.PEAR, otu.table.mock6.PEAR.rel)
colnames(otu.table.mock6.PEAR.both) <- c("Counts", "Rel_abund")
otu.table.mock6.PEAR.both <- otu.table.mock6.PEAR.both[ sort(rownames(otu.table.mock6.PEAR.both)), ]





# mock.fqj
# change NAs to unclassified
tax.mock.fqj.fixNames <- tax.mock.fqj
tax.mock.fqj.fixNames[is.na(tax.mock.fqj.fixNames)] <- "unclassified"
tax.mock.fqj.fixNames <- t(apply(tax.mock.fqj.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.mock.fqj.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                          paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                          all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.mock.fqj.fixNames) <- colnames(tax.mock.fqj)
# get otu table
otu.table.mock.fqj <- seqtab.nochim.mock.fqj
colnames(otu.table.mock.fqj) <- unname(tax.mock.fqj.fixNames[colnames(seqtab.nochim.mock.fqj) , "Genus"])
non.bact.mock.fqj <- tax.mock.fqj.fixNames[ tax.mock.fqj.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
otu.table.mock.fqj <- t(as.matrix(otu.table.mock.fqj[ , ! colnames(otu.table.mock.fqj) %in% non.bact.mock.fqj]))
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table.mock.fqj <- t(otu.table.mock.fqj %*% sapply(unique(colnames(otu.table.mock.fqj)),"==", colnames(otu.table.mock.fqj)))
otu.table.mock.fqj.rel <- apply(otu.table.mock.fqj, 2, function(x) round(100 * x/sum(x), 2))
otu.table.mock.fqj.both <- cbind(otu.table.mock.fqj, otu.table.mock.fqj.rel)
colnames(otu.table.mock.fqj.both) <- c("Counts", "Rel_abund")
otu.table.mock.fqj.both <- otu.table.mock.fqj.both[ sort(rownames(otu.table.mock.fqj.both)), ]






# mock4.fqj
# change NAs to unclassified
tax.mock4.fqj.fixNames <- tax.mock4.fqj
tax.mock4.fqj.fixNames[is.na(tax.mock4.fqj.fixNames)] <- "unclassified"
tax.mock4.fqj.fixNames <- t(apply(tax.mock4.fqj.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.mock4.fqj.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                         paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                         all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.mock4.fqj.fixNames) <- colnames(tax.mock4.fqj)
# get otu table
otu.table.mock4.fqj <- seqtab.nochim.mock4.fqj
colnames(otu.table.mock4.fqj) <- unname(tax.mock4.fqj.fixNames[colnames(seqtab.nochim.mock4.fqj) , "Genus"])
non.bact.mock4.fqj <- tax.mock4.fqj.fixNames[ tax.mock4.fqj.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
otu.table.mock4.fqj <- t(as.matrix(otu.table.mock4.fqj[ , ! colnames(otu.table.mock4.fqj) %in% non.bact.mock4.fqj]))
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table.mock4.fqj <- t(otu.table.mock4.fqj %*% sapply(unique(colnames(otu.table.mock4.fqj)),"==", colnames(otu.table.mock4.fqj)))
otu.table.mock4.fqj.rel <- apply(otu.table.mock4.fqj, 2, function(x) round(100 * x/sum(x), 2))
otu.table.mock4.fqj.both <- cbind(otu.table.mock4.fqj, otu.table.mock4.fqj.rel)
colnames(otu.table.mock4.fqj.both) <- c("Counts", "Rel_abund")
otu.table.mock4.fqj.both <- otu.table.mock4.fqj.both[ sort(rownames(otu.table.mock4.fqj.both)), ]






# mock6.fqj
# change NAs to unclassified
tax.mock6.fqj.fixNames <- tax.mock6.fqj
tax.mock6.fqj.fixNames[is.na(tax.mock6.fqj.fixNames)] <- "unclassified"
tax.mock6.fqj.fixNames <- t(apply(tax.mock6.fqj.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.mock6.fqj.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                          paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                          all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.mock6.fqj.fixNames) <- colnames(tax.mock6.fqj)
# get otu table
otu.table.mock6.fqj <- seqtab.nochim.mock6.fqj
colnames(otu.table.mock6.fqj) <- unname(tax.mock6.fqj.fixNames[colnames(seqtab.nochim.mock6.fqj) , "Genus"])
non.bact.mock6.fqj <- tax.mock6.fqj.fixNames[ tax.mock6.fqj.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota"), "Genus"]
otu.table.mock6.fqj <- t(as.matrix(otu.table.mock6.fqj[ , ! colnames(otu.table.mock6.fqj) %in% non.bact.mock6.fqj]))
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table.mock6.fqj <- t(otu.table.mock6.fqj %*% sapply(unique(colnames(otu.table.mock6.fqj)),"==", colnames(otu.table.mock6.fqj)))
otu.table.mock6.fqj.rel <- apply(otu.table.mock6.fqj, 2, function(x) round(100 * x/sum(x), 2))
otu.table.mock6.fqj.both <- cbind(otu.table.mock6.fqj, otu.table.mock6.fqj.rel)
colnames(otu.table.mock6.fqj.both) <- c("Counts", "Rel_abund")
otu.table.mock6.fqj.both <- otu.table.mock6.fqj.both[ sort(rownames(otu.table.mock6.fqj.both)), ]



# # change NAs to unclassified
# tax.mock6.fqj[is.na(tax.mock6.fqj)] <- "unclassified"
# for (tl in colnames(tax.mock6.fqj)) {
#   t.abbr <- substr(tl,1,1) # letter to use as part of identifier for unclassified taxa
#   tax.mock6.fqj[ , tl] <- sapply(1:nrow(tax.mock6.fqj), function(x) ifelse(tax.mock6.fqj[x, tl]=="unclassified",
#                                                                            sprintf("unclassified.%s%s", t.abbr, x),
#                                                                            tax.mock6.fqj[x, tl]))
# }
# otu.table.mock6.fqj <- seqtab.nochim.mock6.fqj
# colnames(otu.table.mock6.fqj) <- unname(tax.mock6.fqj[colnames(seqtab.nochim.mock6.fqj) , "Genus"])
# otu.table.mock6.fqj <- otu.table.mock6.fqj %*% sapply(unique(colnames(otu.table.mock6.fqj)),"==", colnames(otu.table.mock6.fqj))
# otu.table.mock6.fqj.rel <- apply(otu.table.mock6.fqj, 1, function(x) 100 * x/sum(x))
# otu.table.mock6.fqj.both <- cbind(t(otu.table.mock6.fqj), otu.table.mock6.fqj.rel)
# colnames(otu.table.mock6.fqj.both) <- c("Counts.uparse", "Rel_abund.uparse")
# rownames(otu.table.mock6.fqj.both)[ ! rownames(otu.table.mock6.fqj.both) %in% mock.genera ]
# mock.genera[ ! mock.genera %in% rownames(otu.table.mock6.fqj.both) ]


