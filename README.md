# SLL1_analyses
Code for analyzing the data of the 1st edition of Saca La Lengua

The primary scripts for the analyses are:
- **src/dada2_pipeline.R** for read filtering and taxonomy assignment
- **src/SLL_phyloseq.R** for statistical analyses and figure generation.


The raw fastq files can be found in the Sequence Read Archive (SRA) with the accession number **PRJNA427101** and can be found here: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA427101


The additional data necessary for the analyses are:
- **data/SLL.OTU.counts.csv** - a table of taxon counts per sample
- **data/SLL.questionnaire_responses.csv** - a table of the metadata used for the analyses
- **data/SLL.taxonomy_table.csv** - the taxonomy table
