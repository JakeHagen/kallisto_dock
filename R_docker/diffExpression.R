library("optparse")
library("tximport")
library("readr")
library("DESeq2")
library("jsonlite")

option_list <- list(
  make_option(c("-m", "--meta"), action="store", type="character", default=NULL,
             help="design matrix for experiment in csv format"),
  make_option(c("-c", "--contrast"), action="store", type="character", default=NULL,
               help="contrast matrix for experiment in csv format"),
  make_option(c("-o", "--output"), action="store", type="character",
               help="folder to output results"),
  make_option(c("-d", "--transcriptDir"), action="store", type="character",
               help="directory that contains the kallisto output")
               )

opt = parse_args(OptionParser(option_list=option_list))

metaData <- fromJSON(opt$meta)
#read.csv(opt$meta)
#rownames(preMetaData) <- preMetaData[,1]
#print(metaData)
#rownames(metaData) <- metaData[,1]
#print(metaData)
#metaData <- DataFrame(metaData[,-1])
# metaData <- data.frame(condition=c('control', 'control', 'control',
#                                    'mars2', 'mars2', 'mars2'))
print(metaData)
# rownames(metaData) <- list.files("/home/jake/rna_data/mars2_gastroc/kallisto/")

# contrast <- read.csv(opt$contrast)
contrast <- c("condition", "mars2", "control")
outputDir <- opt$output
transcriptDir <- opt$transcriptDir
transcriptFiles <- paste0(transcriptDir, list.files(transcriptDir), "/abundance.tsv")
names(transcriptFiles) <- list.files(transcriptDir)

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

des_padj05 <- subset(des_full, des_full$padj <= .05)

dir.create(outputDir)

write.table(des_padj05, file = paste0(outputDir, "/", contrast[2],
                                       "_", contrast[3], "_", "gene_padj05.csv"),
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(des_full, file = paste0(outputDir, "/", contrast[2],
                                    "_", contrast[3], "_", "gene_full.csv"),
            row.names = FALSE, sep = ",", quote = FALSE)
