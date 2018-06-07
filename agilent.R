############## ANALYSIS GSE80832 DATASET ################

library(biomaRt)
library(limma)
library(HsAgilentDesign026652.db)

# TIPS BEFORE START
# 1. Download the raw data. Agilent files are generated on tsv format (.txt).
# Before any procedure, certify yourself that they are unziped into the same folder.
# 2. Try to find the platform annotation package on Bioconductor. It may be needed.
# 3. It is always good having a file with sample group informations. For this experiment
# we have the "samples.txt".

# Get the SDRF file in Array Express experiment page. 
# This file is a metadata with informations about the experiment.
SDRF <- read.delim("~/Dropbox/Trab_Dalmolin/primeiro_teste/GSE80832/sdrf.txt",
                   check.names = FALSE, stringsAsFactors = FALSE)

# Read the files
raw <- read.maimages(SDRF[,"Array Data File"], 
                     path = "~/Dropbox/Trab_Dalmolin/primeiro_teste/GSE80832/",
                     source = "agilent", green.only = TRUE)

# Normalize data
norm <- backgroundCorrect(raw, method = "normexp")
norm <- normalizeBetweenArrays(norm, method = "quantile")
norm.ave <- avereps(norm, ID = norm$genes$ProbeName)

# Get and process the expression set
eset <- ExpressionSet(assayData = norm.ave$E, annotation = "HsAgilentDesign026652.db")
my_exp <- exprs(eset)
my_exp <- cbind(my_exp, norm.ave$genes[, c(3,4)])
my_exp <- my_exp[my_exp$ControlType == 0,]
probes <- make.names(my_exp$ProbeName, unique = T)
row.names(my_exp) <- probes
my_exp[, 55] <- NULL

########### QUALITY ASSESSMENT ############

# MAplot raw data
par(mfrow = c(2,3))
for(i in seq_along(raw$targets$FileName)) {
    limma::plotMA.EListRaw(raw, array = i)
    lines(c(1,15), c(0,0), col = "blue", lwd = 2)
}

# MAplot normalized data
par(mfrow = c(2,3))
for(i in seq_along(raw$targets$FileName)) {
    limma::plotMA.EList(norm.ave, array = i)
    lines(c(1,15), c(0,0), col = "blue", lwd = 2)
}

### PCA
targets <- factor(SDRF$`Comment [Sample_source_name]`)
# Retirar os NAs
data.exprs <- my_exp[, -55]
data.exprs <- as.data.frame(na.omit(data.exprs))

# Nomear as colunas da tabela de expressao de acordo com os grupos
colnames(data.exprs) <- SDRF$`Comment [Sample_source_name]`

# Obter a matrix transposta
tdata.exprs <- t(data.exprs)

# Obter a matriz design
lev <- levels(targets)
design <- model.matrix(~ 0 + targets)
colnames(design) <- c("mock_0", "mock_18", "mock_30", "mock_8",
                      "mucin_0", "mucin_18", "mucin_30", "mucin_8",
                      "wild_0", "wild_18", "wild_30", "wild_8")

# Estabelecer o esquema de cores
color.code <- rainbow(12)

# calcular a PCA
pca <- prcomp(tdata.exprs)

# Plot 1 - Todos os grupos
par(mfrow = c(1,1))
samples <- read.delim("~/Dropbox/Trab_Dalmolin/primeiro_teste/GSE80832/samples.txt")$Names
plot(x = pca$x[,1:2], pch=19, col=color.code[samples], main = "PCA after normalization")
text(pca$x[,1]-1,pca$x[,2]-1, samples, cex=0.7) # nomear de acordo com as amostras presentes no arquivo 'samples.txt'

# Plot 2 - Cepas 
samples <- read.delim("~/Dropbox/Trab_Dalmolin/primeiro_teste/GSE80832/samples.txt")$CellTypes
plot(x = pca$x[,1:2], pch = 19, col = color.code[samples], main = "PCA after normalization")
text(pca$x[,1]-1,pca$x[,2]-1, samples, cex=0.7)

#plot 3
samples <- as.factor(read.delim("~/Dropbox/Trab_Dalmolin/primeiro_teste/GSE80832/samples.txt")$Hour)
plot(x = pca$x[,1:2], pch=19, col=color.code[samples], main = "PCA after normalization")
text(pca$x[,1]-1,pca$x[,2]-1, samples, cex = 0.7)

### Dendrogram
rownames(tdata.exprs) <- samples
d <- dist(as.matrix(tdata.exprs), method = "euclidean")
clusters <- hclust(d)
par(mar = c(1, 1, 1, 1))
plot(clusters)

############ DIFFERENTIAL EXPRESSION ###############

# Sample names
samples <- read.delim(
    file = "~/Dropbox/Trab_Dalmolin/primeiro_teste/GSE80832/samples.txt")$Names

# Map the probes with biomaRt
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
att <- listAttributes(ensembl)
probes <- make.names(my_exp$ProbeName, unique = T)
ids <- getBM(attributes = c("agilent_wholegenome_4x44k_v2",
                            "entrezgene",
                            "hgnc_symbol"),
             filters = 'agilent_wholegenome_4x44k_v2',
             values = probes,
             mart = ensembl)

entrez <- ids$entrezgene
symbol <- ids$hgnc_symbol

# Alternatively, we can map probes with the platform database package 
library(HsAgilentDesign026652.db)
x <- HsAgilentDesign026652ENTREZID
entrez <- unlist(as.list(x[my_exp$ProbeName]))
x <- HsAgilentDesign026652SYMBOL
symbol <- unlist(as.list(x[my_exp$ProbeName]))

# Assing groups
targets <- factor(SDRF$`Comment [Sample_source_name]`)

# Build model matrix
lev <- levels(targets)
design <- model.matrix(~ 0 + targets)
colnames(design) <- c("mock_0", "mock_18", "mock_30", "mock_8",
                      "mucin_0", "mucin_18", "mucin_30", "mucin_8",
                      "wild_0", "wild_18", "wild_30", "wild_8")

# Fit linear model
fit <- lmFit(my_exp[, -55], design)

# Assign the mapped ids into the fitted model
fit$genes$entrez <- entrez
fit$genes$symbol <- symbol

# Establish contrasts in order to compare groups:
# Determine all possible combinations
groups <- combn(c("mock_0", "mock_18", "mock_30", "mock_8",
                "mucin_0", "mucin_18", "mucin_30", "mucin_8",
                "wild_0", "wild_18", "wild_30", "wild_8"), m = 2)
comb <- apply(groups, 2, paste, collapse = "-")
contrasts <- makeContrasts(contrasts = comb, levels = design)

# Bayes stats
ct.fit <- eBayes(contrasts.fit(fit, contrasts))
res.fit <- decideTests(ct.fit, method = "global", adjust.method = "BH", p.value = 0.001)

# Combine results into a single dataset
SH.limma <- data.frame(logFC = ct.fit$coef, p.value = ct.fit$p.value, 
                       entrez = ct.fit$genes$entrez, symbol = ct.fit$genes$symbol,
                       degenes = unclass(res.fit),
                       stringsAsFactors = FALSE)

# Select only the genes differentially expressed
features <- rowSums(res.fit != 0) > 0
features <- names(features)[features]

# Filter the expression table to obtain the differentially expressed probes
DEexp <- my_exp[features, ]
DElimma <- SH.limma[rownames(DEexp), ]
DElimma <- DElimma[complete.cases(DElimma), ]

# save results
save(DElimma, DEexp, SH.limma, features, 
     file = "~/Dropbox/Trab_Dalmolin/primeiro_teste/GSE80832/Exp_values.RData")






