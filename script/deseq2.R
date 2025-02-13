library(DESeq2)

#File path
countfile <- "bioproject/counts/featureCounts_org.txt"
up <- "bioproject/deseq2/deseq2_up.txt"
down <- "bioproject/deseq2/deseq2_down.txt"

#DataFrame
df <- read.table(countfile, header = TRUE, row.names = 1, check.names = FALSE)
df <- df[, -(1:5)]

df1 <- colnames(df)
condition <- factor(c("A", "A", "A", "B", "B", "B"))
col <- data.frame(row.names = df1, condition = condition)

#DESeq2 analysis
df2 <- DESeqDataSetFromMatrix(countData = df, colData = col, design = ~condition)
df2 <- DESeq(df2, fitType = "local")

#Results
df3 <- results(df2, contrast = c("condition", "B", "A"))

result_up <- df3[which(df3$log2FoldChange >= 2 & df3$padj < 0.05), ]
result_down <- df3[which(df3$log2FoldChange <= -2 & df3$padj < 0.05), ]

result_up <- cbind(Geneid = rownames(result_up), as.data.frame(result_up))
result_down <- cbind(Geneid = rownames(result_down), as.data.frame(result_down))

write.table(result_up, file = up, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(result_down, file = down, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
