setwd("D:\\IOB\\Projects\\FED_062723\\WGCNA\\26_11_24\\")
data = read.table("expressionProfile_deicemt_allDatasets_adjustedCount_112824.txt", sep = '\t',header = TRUE)
head(data)
dim(data)
rownames(data) <- data$hgnc_symbol
data <- data[,-1]
head(data)
colnames(data)
rownames(data)


dim(data)
class(data)

library("WGCNA")
datExpr = data.matrix(data)
gsg <- goodSamplesGenes(t(datExpr))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples) 
data <- data[gsg$goodGenes == TRUE,]
outlier = data[gsg$goodGenes == FALSE,]

head(data)
nrow(data)

library('DESeq2')
#dim(data)
condition1 = factor(c(rep("control",29), rep("fed",65)))
countdata1 = as.matrix(data)
coldata1 = data.frame(row.names = colnames(countdata1), condition1)

ddsFull <- DESeqDataSetFromMatrix(countData = countdata1,
                              colData = coldata1,
                              design = ~ condition1)
dds = DESeq(ddsFull)

vsd <- varianceStabilizingTransformation(dds)
wpn_vsd = getVarianceStabilizedData(dds)

expr_normalized = wpn_vsd

dim(expr_normalized)

input_mat = t(expr_normalized)

input_mat[1:5,1:10]

library('WGCNA')
library('flashClust')
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,1))
power = c(c(1:10), seq(from = 12, to = 50, by = 2))
sft = pickSoftThreshold(input_mat,
                  powerVector = power,
                  verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free topolofy Fit, signed R^2", type = "n",
     main = paste("Scale independence - Tenon"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = power, col="red");
abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity - tenon"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power,col="red")

softPower1 = 9
adjacency1 = adjacency(input_mat, power = softPower1, type = "signed")
dissTOM1   = 1-TOMsimilarity(adjacency1, TOMType="signed")
geneTree1  = flashClust(as.dist(dissTOM1), method="average")
par(mfrow=c(1,1))
minModuleSize = 50
dynamicMods1 = cutreeDynamic(dendro = geneTree1, distM = dissTOM1,
                             deepSplit = 0, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods1)
dynamicColors1 = labels2colors(dynamicMods1)
colors1 = table(dynamicColors1)
colors1
plotDendroAndColors(geneTree1, dynamicColors1, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Tenon")
par(mar=c(5.1, 4.1, 4.1, 2.1))
traitdata1 = read.csv("RNASeq_traits.txt", sep = '\t',header = TRUE)
patient1 = rownames(input_mat)
traitRows1 = match(patient1, traitdata1$Sample_ID)
datTraits1 = traitdata1[traitRows1, -1]
names(datTraits1)
nGenes1 = ncol(input_mat)
nSamples1 = nrow(input_mat)
MEs0 = moduleEigengenes(input_mat, dynamicColors1)$eigengenes
MEs1 = orderMEs(MEs0)
MEs1
moduleTraitCor1 = cor(MEs1, datTraits1, use = "p")
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamples1)
textMatrix1 = paste(signif(moduleTraitCor1, 2), "\n(",signif(moduleTraitPvalue1, 1), ")", sep = "")
dim(textMatrix1) = dim(moduleTraitCor1)

labeledHeatmap(Matrix = moduleTraitCor1,
              xLabels = names(datTraits1),
              yLabels = names(MEs1),
              ySymbols = NULL,
              colorLabels = FALSE,
              colors = blueWhiteRed(50),
              textMatrix = textMatrix1,
              setStdMargins = FALSE,
              cex.text = 1,
              zlim = c(-1,1), 
              main = paste("Module- Trait relationships - Tenon"), ylab = "Gene Expression-Based Modules")

par(mfrow = c(1,1))
fed = as.data.frame(datTraits1$fed)
names(fed) = "fed"
modNames1 = substring(names(MEs1), 3)

geneModuleMembership1 = as.data.frame(cor(input_mat, MEs1, use = "p"))
MMPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership1), nSamples1))
names(geneModuleMembership1) = paste("MM", modNames1, sep="")
names(MMPvalue1) = paste("p.MM", modNames1, sep="")
geneTraitSignificance1 = as.data.frame(cor(input_mat, fed, use = "p"))
GSPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance1), nSamples1))
names(geneTraitSignificance1) = paste("GS.", names(fed), sep="")
names(GSPvalue1) = paste("p.GS.", names(fed), sep="")

module = "blue"
column = match(module, modNames1)  
moduleGenes = dynamicColors1==module

verboseScatterplot(abs(geneModuleMembership1[moduleGenes, column]),
                   abs(geneTraitSignificance1[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for fed",
                   cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, col = module)

#Import to cytoscape
TOM1 = TOMsimilarityFromExpr(input_mat, power = 9)
modules = c("blue")
genes = colnames(input_mat)
inModule = is.finite(match(dynamicColors1, modules));
modGenes = genes[inModule];
modTOM = TOM1[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("ControlVsFED_CytoscapeInput-edges-", paste(modules, collapse="-"), "_0.txt", sep=""),
                               nodeFile = paste("ControlVsFED_CytoscapeInput-nodes-", paste(modules, collapse="-"), "_0.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               nodeAttr = dynamicColors1[inModule])

