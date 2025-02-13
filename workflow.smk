samples = ["sample_0", "sample_1", "sample_2", "sample_3", "sample_4", "sample_5"]

rule all:
    input:
        "bioproject/alignment/star_completed.flag",
        "bioproject/counts/featureCounts_output.txt",
        "bioproject/deseq2/deseq2_up.txt",
        "bioproject/deseq2/deseq2_down.txt"

rule run_fastqc:
    input:
        expand("bioproject/fastqc/{sample}_fastqc.zip", sample=samples),
        expand("bioproject/fastqc/{sample}_fastqc.html", sample=samples)

rule run_cutadapt:
    input:
        expand("bioproject/filtered/{sample}_filtered.fastq", sample=samples)

rule run_starindex:
    input:
        "bioproject/alignment/STAR_index"

rule run_starmapping:
    input:
        expand("bioproject/alignment/{sample}_Aligned.out.bam", sample=samples)

rule fastqc:
    input:
        "bioproject/rawdata/{sample}.fastq"
    output:
        zip="bioproject/fastqc/{sample}_fastqc.zip",
        html="bioproject/fastqc/{sample}_fastqc.html"
    shell:
        "fastqc {input} --outdir=bioproject/fastqc"

rule cutadapt:
    input:
        "bioproject/rawdata/{sample}.fastq",
        adapter="bioproject/reference/illumina_adapter.fa"
    output:
        "bioproject/filtered/{sample}_filtered.fastq"
    params:
        options="--minimum-length 25 --quality-cutoff 30 --cut 10"
    shell:
        "cutadapt {params.options} -a file:{input.adapter} -o {output} {input[0]}"

rule starindex:
    input:
        genome_fasta="bioproject/reference/genome.fa",
        annotation="bioproject/reference/annotation.gtf"
    output:
        directory("bioproject/alignment/STAR_index")
    shell:
        """
        STAR --runThreadN 8 \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome_fasta} \
             --sjdbGTFfile {input.annotation} \
             --sjdbOverhang 100 \
             --genomeSAindexNbases 11
        """

rule starmapping:
    input:
        reads="bioproject/filtered/{sample}_filtered.fastq",
        index="bioproject/alignment/STAR_index"
    output:
        bam="bioproject/alignment/{sample}_Aligned.out.bam",
        sam="bioproject/alignment/{sample}_Aligned.out.sam"
    params:
        prefix="bioproject/alignment/{sample}_"
    shell:
        """
        STAR --runThreadN 8 \
             --genomeDir {input.index} \
             --readFilesIn {input.reads} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM Unsorted \
             --outSAMattributes All

        samtools view -h -o {output.sam} {output.bam}
        """

rule featurecounts:
    input:
        expand("bioproject/alignment/{sample}_Aligned.out.bam", sample=samples)
    output:
        original="bioproject/counts/featureCounts_org.txt",
        modified="bioproject/counts/featureCounts_output.txt"
    params:
        annotation="bioproject/reference/annotation.gtf"
    shell:
        """
        featureCounts -T 8 \
                      -a {params.annotation} \
                      -o {output.original} \
                      -s 1 {input}

        echo -e "Geneid\tsample_0\tsample_1\tsample_2\tsample_3\tsample_4\tsample_5" > {output.modified}
        awk 'NR>2 {{print $1, $7, $8, $9, $10, $11, $12}}' OFS="\t" {output.original} >> {output.modified}
        """

rule deseq2:
    input:
        counts="bioproject/counts/featureCounts_org.txt"
    output:
        up="bioproject/deseq2/deseq2_up.txt",
        down="bioproject/deseq2/deseq2_down.txt"
    shell:
        """
        Rscript -e 'library(DESeq2);
        countfile <- "{input.counts}";
        up <- "{output.up}";
        down <- "{output.down}";
        df <- read.table(countfile, header = TRUE, row.names = 1, check.names = FALSE);
        df <- df[, -(1:5)];
        df1 <- colnames(df);
        condition <- factor(c("A", "A", "A", "B", "B", "B"));
        col <- data.frame(row.names = df1, condition = condition);
        df2 <- DESeqDataSetFromMatrix(countData = df, colData = col, design = ~condition);
        df2 <- DESeq(df2, fitType = "local");
        df3 <- results(df2, contrast = c("condition", "B", "A"));
        result_up <- df3[which(df3$log2FoldChange >= 2 & df3$padj < 0.05), ];
        result_down <- df3[which(df3$log2FoldChange <= -2 & df3$padj < 0.05), ];
        result_up <- cbind(Geneid = rownames(result_up), as.data.frame(result_up));
        result_down <- cbind(Geneid = rownames(result_down), as.data.frame(result_down));
        write.table(result_up, file = up, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE);
        write.table(result_down, file = down, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE);'
        """