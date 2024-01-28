# RNA-Seq-Analysis

Developed independently as part of coursework for Bioinformatics Programming and System Management at the University of Edinburgh.

Pipeline taking *Trypanosoma congolense* RNA-Seq data in compressed FASTQ format.
Runs a quality check on the data, aligns the read pairs to a *Trypanosoma congolense* genome (BAM format), generates count data of the number of reads that align to gene-coding regions of the genome, generates TXT files with average expression level per gene, after which the fold-change is calculated between chosen group-wise comaprisons.

Step by step instructions indicated in "Help_manual.pdf". Script files run Bash and are intended to be run on Linux.
