library("optparse")
library("tximport")
library("readr")
library("DESeq2")

option_list <- list(
  make_option(c("-m", "--meta"), action="store", type="character", default=NULL,
             help="design matrix for experiment in csv format"),
  make_option(c("-c", "--contrast"), action="store", type="character", default=NULL,
               help="contrast matrix for experiment in csv format"),
  make_option(c("-o", "--output"), action="store", type="character",
               help="folder to output results"),
  make_option(c("-t", "--transcriptCounts"), action="store", type="character",
               help="transcript counts file for differential expression, should be space between file names,
               should be tvs format")
)

opt = parse_args(OptionParser(option_list=option_list))

# metaData <- read.csv(opt$meta)
metaData <- data.frame(condition=c('control', 'control', 'control',
                                  'mars2', 'mars2', 'mars2'))
rownames(metaData) <- list.files("/home/jake/rna_data/mars2_gastroc_kallisto/")

# contrast <- read.csv(opt$contrast)
contrast <- c("condition", "mars2", "control")
outputDir <- opt$output
#transcriptDirs <- unlist(strsplit(opt$transcriptCounts, split=" "))
transcriptDirs <- paste0("/home/jake/rna_data/mars2_gastroc_kallisto/",
                         list.files("/home/jake/rna_data/mars2_gastroc_kallisto/"))
transcriptFiles <- paste0(transcriptDirs, "/abundance.tsv")
names(transcriptFiles) <- list.files("/home/jake/rna_data/mars2_gastroc_kallisto/")

tx_full <- read_tsv(transcriptFiles[1])
tx_gene_names <- read.table(transcriptFiles[1], sep="|", skip=1)
tx2gene <- data.frame(TXNAME=tx_full$target_id, GENEID=tx_gene_names$V6)

txi <- tximport(transcriptFiles, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)

formula <- as.formula(paste("~ ", paste(colnames(metaData), collapse = "+")))

dds <- DESeqDataSetFromTximport(txi, colData = metaData, design = formula)
dds <- dds[rowSums(counts(dds)) > 0,]
dds <- DESeq(dds)
res <- results(dds, contrast = contrast)
res <- res[order(res$padj),]
 
des_full <- data.frame(res@rownames, res@listData)

des_padj.05 <- subset(des_full, des_full$padj <= .05)

dir.create(output)
write.table(des_full, file = paste0(output, "/gene_full.csv"),
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(des_padj.05, file = paste0(output, "/gene_padj.05.csv"),
            row.names = FALSE, sep = ",", quote = FALSE)