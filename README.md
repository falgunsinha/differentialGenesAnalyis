# Differential gene expression analyis

Perform differential gene expression analysis on RNA-seq data.
The data set contains 6 RNA-seq samples. Samples 0-2 belong to Condition A and samples 3-5 to Condition B.
• Condition A: log2 fold-change >= 2
• Condition B: log2 fold-change <= -2

## Steps:
- fastqc: for quality control of the raw and quality filtered data.
- cutadapt: for quality filtering of the raw reads.
- STAR: to map reads to the genome.
- featureCounts: to obtain reads counts on the gene level.
- DESeq2: to detect differentially expressed genes.

**Note**: the RNA-seq data is single-end. We need to change program options accordingly, i.e. for featureCounts, use the option
-s 1.

## Files description:
- featureCounts_output.txt : Contains the featureCounts results of 6 samples.
- deseq2_up.txt : Contains the subset of DESeq2 results for genes with a log2 fold-change >= 2.
- deseq2_down.txt : Contains the subset of DESeq2 results for genes with a log2 fold-change <= -2.
- workflow.smk : snakemake workflow used to generate the results.
- protocol.ipynb : A description of steps performed and commands executed.
