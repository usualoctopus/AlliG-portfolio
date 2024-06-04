---
layout: post
title: "WGCNA Script for Gene Network Analysis"
date: 2024-06-04
categories: RNA Sequencing Analysis
---

# This is an R markdown containing a full script, visualizations, and explanations for performing Weighted Gene Co-expression Network Analysis (WGCNA) to understand the topology of gene networks across a data set. (Zhang, Horvath; 2005)

## WGCNA can be interpreted as a biologically motivated data reduction scheme that allows for dependency between the resulting components. WGCNA identifies gene modules from RNA sequencing data. Module eigengenes (MEs) are defined as the first principal component (PC).

## Here, I used miRNA sequencing data. Though the n is low with respect to what is recommended (>15), this script nevertheless provides an introduction to and explanation for WGCNA analysis.

***
```
# Load necessary packages and sample files

library(factoextra)
library(tidyverse)
library(devtools)
library(cluster)
library(ggrepel)
library(kableExtra)
library(gplots)
library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE);

# Input samples used for analysis and rename for easier viewing
counts <- read.csv("Partek_miRNA-Reseq_Other_counts_reformat_omit21.csv")
cpm <- read.csv("ALL_miRNAreseq_filtered_CPMs.csv")

# Samples to be used for UT control
ut.control.names <- c(
  "Control_413",
  "Control_392",
  "Control_494",
  "Control_49",
  "Control_43",
  "Control_52",
  "Control_323",
  "Control_449",
  "Control_483",
  "Control_237",
  "Control_492",
  "Control_141",
  "Control_417"
)

ut.control <- rep('control_ut', length(ut.control.names))

num <- 1:length(ut.control.names)
prefix <- rep("Control_")
new.ut.control.names <- paste(prefix,num, sep = "")
new.ut.control.names


# Samples to be used for UT PMDD
# PMDD 21 was omitted based on z-score value 2.8 (standard deviations from the mean) 
# PMDD 98 was omitted for low read count, poor alignment

ut.pmdd.names <- c(
  "PMDD_328",
  "PMDD_406",
  "PMDD_408",
  "PMDD_484",
  "PMDD_491",
  #"PMDD_21",
  #"PMDD_98",
  "PMDD_13",
  "PMDD_486"
)


num <- 1:length(ut.pmdd.names )
prefix <- rep("PMDD_")
new.ut.pmdd.names <- paste(prefix,num, sep = "")
new.ut.pmdd.names

ut.pmdd <- rep('pmdd_ut', length(ut.pmdd.names))


labels <- c(ut.control, ut.pmdd)
labels.dx <- c(rep('control', length(c(ut.control))), rep('pmdd', length(c(ut.pmdd))))


names.new <- c("X", new.ut.control.names, new.ut.pmdd.names)
names <- c("X", ut.control.names, ut.pmdd.names)
names

colnames(cpm)
colnames(cpm)<- c(names.new)
head(cpm)



```

### Data Preprocessing and Exploratory Analysis: Filtering/Normalizing, PCA plot, and Gene Variance Heatmap
```
## Reduce noise: filter out low count genes
# calculate mean expression for each gene
sums <- rowMeans(cpm[,-1]) 
summary(sums) 

# normalize (log-transform), reduce the effect of extreme values, allows for threshold-based filtering
hist(log(sums+1), breaks = 50) # include +1 so don't do log of 0)
summary(log(sums+1))
```

![miRNA_Histogram_logCPM](_posts/r-scripts/assets/images/1ef1563c-5bb7-46d8-9cc0-c7f6748e27ca.png)

The above illustrates the distribution of the unfiltered log-transformed mean expression values of genes across all samples (x-axis) and the number of genes that fall within each bin of log-transformed mean expression values (y-axis, Frequency) before any filtering is applied. This helps decide an appropriate filtering threshold.

```
# retain only genes with a log-transformed mean expression >1.7
cpm.filtered <- cpm[log(sums+1) > 1.7,] 
hist(log(t(cpm.filtered[,-1]+1)))
dim(cpm.filtered)
```

![miRNA_Histogram_logCPMfiltered](assets/images/miRNA_Histogram_logCPMfiltered.png)

This shows the filtered log-transformed gene expression values of genes retained after filtering (>1.7 mean expression value) by the number of expression values across all genes and samples that fall within each bin. Filtering appears to clean up the distribution of gene expression by reducing low-expression genes. Note the heavy tail.

```
# pcaplot
pca <- prcomp(t(cpm.filtered[,-1]))
pca_out <- as.data.frame(pca$x)
pca_out$group <- labels
pca_out$dx <- labels.dx

plot<-ggplot(pca_out,aes(x=PC1,y=PC2,color=dx, size=dx, label=row.names(pca_out) ))
plot<-plot+geom_point()+xlab('PC1 (40.2%)')+ylab('PC2 (23.1%)')
plot
eigs <- pca$sdev^2
# # obtain % variance explained by each PC
# eigs[1] / sum(eigs)
# eigs[2] / sum(eigs)
# # obtain full component summary
summary(pca)
```

#### __PCA Plot, all samples__
![miRNA_PCA.png](assets/images/miRNA_PCA.png)

```
# make a heat map to visualize variation in gene expression
hmcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
# <rowVars() from genefilter package calculates the variances of each row of an array>
x <- cpm.filtered[,-1]
y <- log(x+1)
rv <- rowVars(y)
idx <- order(-rv)[1:50]
idx

# Color-blind friendly (control = green, pmdd = orange)
# Set up the margins: mar = c(bottom, left, top, right)
par(mar = c(10, 4, 6, 2)) # Increase bottom margin for sample names
cols1 <- palette(brewer.pal(6, "Set2"))[as.factor(labels.dx)]
rownames(x) <- cpm.filtered[,1]
head(cbind(colnames(x), cols1))
heatmap.2(as.matrix(log(x[idx,]+1)), labCol=names[-1],
          trace="none", 
          ColSideColors=cols1, 
          margins = c(10, 10),
          col=hmcol)

head(cpm)
names(cpm)

head(cpm)
names(cpm)
```

#### __Heatmap of the top 50 most variable miRNA transcripts__
![miRNA_Top50variance_heatmap.png](assets/images/miRNA_Top50variance_heatmap.png)
This heatmap highlights the top 50 most variable genes after filtering low count & normalizing miR transcripts. Columns indicate samples and rows indicates genes. The colored column bar at the top (green = control, orange = PMDD) illustrates sample diagnosis. This heatmap reveals generally homogenous sample clusters, indicating that diagnostic condition appears to have a consistent effect on gene expression.


#### Check count data for missing entries, entries with weights below a threshold, and zero-variance genes. gsg returns a list of samples and genes that pass criteria on maximum number of missing or low weight values. If necessary, filtering is iterated (unnecessary for this dataset).

```
# load in all count data (but no gene list)
datExpr0 <- as.data.frame(t(cpm[,-1])) 
rownames(datExpr0)
colnames(datExpr0) <- cpm$X
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK # if this returns TRUE, all genes and samples have passed the test
# str(gsg) shows a shortened list: 1) goodGenes 2) goodSamples 3) all0K

# transposed columns/rows
dim(datExpr0)
dim(cpm)

# Hierarchical clustering dendrogram to visualize outliers.
# hclust() performs HC, which is calculated using dist() as the distance metric (here we use the default: Euclidean distance.
# Method = average; i.e., average linkage method (UPGMA- Unweighted Pair Group Method with Arithmetic Mean).
sampleTree = hclust(dist(datExpr0), method = "average");
baseHeight = 108000
sizeGrWindow(4,3)
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# # if you want to remove outliers
# abline(h = 108000, col = "red");
# clust = cutreeStatic(sampleTree, cutHeight = 108000, minSize = 5)
# table(clust)
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]

datExpr <- datExpr0

# excluded samples
# rownames(datExpr0)[!as.logical(clust)]

# genes and samples that will be used for the rest of the analysis 
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

nGenes
nSamples

```
#### Load in clinical data and create a sample dendrogram with trait heatmap.

```
traitData <- as.data.frame(cbind(names.new[-1], labels.dx))
colnames(traitData) <- c('Sample', 'Diagnosis') 
head(traitData)
rownames(datExpr)
str(traitData)

rownames(datExpr)
traitData$Sample

head(traitData)
rownames(traitData) <- traitData$Sample
traitData <- traitData[-1]
head(traitData)
# head(traitData)
#              Diagnosis
# UT_control_1   control
# UT_control_2   control
# UT_control_3   control
# UT_control_4   control
# UT_control_5   control

traitData$Diagnosis <- as.numeric((as.factor(traitData$Diagnosis)))

collectGarbage()
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData, signed = TRUE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")
```

#### Heirarchical sample clustering with corresponding diagnoses
![SampleDendro_TraitHeatmap.png](assets/images/output_4_6.png)
The Height metric on the y-axis corresponds to the Euclidean distance, which is calculated using Unweighted Pair Group Method with Arithmetic Mean (UPGMA). That is, the distance between two clusters is defined as the average distance between all pairs of samples of samples in the two clusters. The samples/sample clusters at the lowest height are the most similar to each other, whereas those at the heighest are the least similar to each other.
Note: I inexplicably renamed the samples.

### Begin WGCNA analysis: construct a signed hybrid network

#### Determine soft-thresholding power and fit the scale-free topology
```
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers,
                        networkType="signed hybrid", verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

```

![SoftPowerThreshold12.png](assets/images/SoftPowerThreshold12.png)
The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis). Constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency (Zhang, Horvath; 2005). The soft thresholding power is chosen based on the criterion of approximate scale-free topology. For these data a power threshold of 4 was chosen, it is the lowest power for which the scale-free topology fit index reaches 0.90.
For these data, a power threshold of 12 is best.

#### Construct Gene Clustering Dendrogram with Module Colors (Merged Dynamic), TOM-based dissimilarity
```
# calculate the adjacencies using power 12
softPower = 12;
adjacency = adjacency(datExpr, power = softPower, type="signed hybrid");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
minModuleSize = 10;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram with colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Hybrid",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Hybrid", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
```

![miRNA_DynamicHybrid_dendro.png](assets/images/miRNA_DynamicHybrid_dendro.png)
To minimize effects of noise and spurious associations, the adjacency is transformed into a signed Topological Overlap dissimilarity Matrix. Module eigengenes are calculated and clustered. Dynamic tree cut identifies modules whose expression profiles are very similar. It's ideal to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them based on their correlation. 

The cut height is set to 0.25, corresponding to a correlation of 0.75, to merge. 
Minimum module size was set to 10.
cutreeDynamic was set to the Hybrid method.
 - 29 modules before clustering, merging based on eigengene
 - 15 modules after clustering, merging

A hierarchical clustering tree (dendrogram) detects clusters of gene expression profiles. In the dendrogram, each leaf corresponds to a gene. Branches of the dendrogram group together are densely interconnected, highly co-expressed genes. Module identification amounts to the identification of individual branches (”cutting branches off the dendrogram”). There are several branch cutting methods, here we used the dynamic hybrid algorithm. This bottom-up algorithm consists of two steps. First, branches that satisfy criteria for being clusters are detected:

   1) Defined cluster minimum (here it was set to 10).
    
   2) Objects too far from a cluster are excluded even if they belong to the same branch.
    
   3) Each cluster separated from its surroundings by a gap.
    
   4) The _core_ (i.e., the lowest-merged objects in the cluster) of each cluster should be tightly connected.

Second, previously unassigned objects are tested for sufficient proxitimity to clusters detected in the first step; if the closest cluster is close enough in the precise sense we define, the gene is assigned to that cluster. This second step could be considered a modified k-medoid partitioning (PAM), hence a hybrid of hierarchical clustering and modified PAM.

#### Relate modules to traits
```
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traitData, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

sizeGrWindow(10,6)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# look at which module has the strongest relationship with diagnosis
moduleTraitPvalue[order(moduleTraitPvalue[,1]),1]
```

![miRNA_ModTraitRelationship.png](assets/images/miRNA_ModTraitRelationship.png)
Each row corresponds to a module eigengene, each column to diagnosis. Each cell contains the corresponding correlation coefficient (R value, top) with respect to diagnosis, and associated p-value (bottom, parentheses). The lightcyan and darkgrey modules are significantly different (p<0.05) compared to diagnosis.

#### Determine Gene Significance (GS) and Module Membership (MM)
#### Calculate module eigenvalue, MM p-value, and GS p-value for each module. Intramodule analysis identifies genes with high GS and MM by diagnosis. Scatterplots can be made for each module. 
Here, quantify associations of individual genes with the trait of interest (diagnosis) by defining GS as (the absolute value of) the correlation between the gene and the trait. For each module, also define a quantitative measure of MM as the correlation of the module eigengene and the gene expression profile. This allows quantification the similarity of all genes on the array to every module.
```
diagnosis <- as.data.frame(traitData$Diagnosis)
names(diagnosis) <- 'diagnosis'
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance.dx = as.data.frame(cor(datExpr, diagnosis, use = "p"));

GSPvalue.dx = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance.dx), nSamples));

names(geneTraitSignificance.dx) = paste("GS.", names(diagnosis), sep="");

names(GSPvalue.dx) = paste("p.GS.", names(diagnosis), sep="")

# All miRNAs, each module eigenvalue
#write.csv(geneModuleMembership, "miRNA_module_values_all.csv")
# All miRNAs, each module pvalue
#write.csv(MMPvalue, "miRNA_module_pvalues_all.csv")
# All miRNA GS diagnosis pvalue
#write.csv(GSPvalue.dx, "miRNA_GSdx_pvalues.csv")

# Intramodular analysis: identify genes with high GS and MM
# identify genes that have a high significance for dx as well as high MM in interesting modules.

# diagnosis modules: lightcyan, darkgrey (GS and MM in Absolute Value)
module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance.dx[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for pmdd diagnosis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# significant diagnosis modules: lightcyan, darkgrey (GS and MM NOT in absolute value)
module = "darkgrey"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mar = c(4.5,6,4,2));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance.dx[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "GS for PMDD diagnosis",
                   main = paste("Module membership vs. gene significance\n"), cex = 1.5,
                   cex.main = 1.5, cex.lab = 1.45, cex.axis = 1.25, bg = module,
                   col="black", abline.color = "black", pch = 21, lwd=1, abline = TRUE)

```

![miRNA_lightcyan_MMvGS.png](assets/images/miRNA_lightcyan_MMvGS.png)
lightcyan (55 genes)
p-value = 0.019

![miRNA_darkgrey_MMvGS.png](assets/images/miRNA_darkgrey_MMvGS.png)
darkgrey (19 genes)
p-value = 0.048

![miRNA_lightcyan_MMvGS_NoAbsValue.png](assets/images/miRNA_lightcyan_MMvGS_NoAbsValue.png)
lightcyan (55 genes)
p-value = 0.019

![miRNA_darkgrey_MMvGS_NoAbsValue.png](assets/images/miRNA_darkgrey_MMvGS_NoAbsValue.png)
darkgrey (19 genes)
p-value = 0.048

GS and MM quantify associations of individual genes with our trait of interest (diagnosis). GS is defined as the (Pearson) correlation coefficient between the gene and the trait. Here, the top two plots are shown using the absolute value of the GS, as is standard. The bottom two plots show the same two plots without absolute value GS, which I find easier to interpret. For each module, quantitative measure of MM is defined as the correlation between the module eigengene and the gene expression profile. This allows for quantification the similarity of all genes on the array to every module.

#### Summary Output of Network Analysis
```
# return all 'probes' (genes) in analysis
# names(datExpr)
# return all genes belongning to a specific module
# names(datExpr)[moduleColors=="lightcyan"]

# list all genes in all modules, with signficance
moduleMap <- as.data.frame(moduleColors)
rownames(moduleMap) <- colnames(datExpr)
head(moduleMap)
moduleMap <- cbind(moduleMap, geneTraitSignificance.dx, abs(geneTraitSignificance.dx))

# determine how many genes are in a module
dim(moduleMap[moduleMap$moduleColors=='lightcyan',])

#write.csv(moduleMap, "ModuleMap_clust5_MM10.csv")

# find a gene in the modules
moduleMap['hsa-miR-503-5p',]

# find genes with highest membership
module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#geneModuleMembership[moduleGenes, column]

module.info.lightcyan <- cbind(moduleMap[moduleMap$moduleColors==module,], abs(geneModuleMembership[moduleGenes, column]))
#write.csv(module.info.lightcyan, "lightcyan.csv")
```
#### Data Visualizations
```
# Calculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, networkType = "signed hybrid", 
                            TOMType = "signed", power = 12);
dissTOM <- 1-TOM

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
#save(TOM, dissTOM, plotTOM, file = "TOM_dissTOM_plotTOM-auto.RData")

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
```
![Network_heatmap_allgenes.png](assets/images/Network_heatmap_allgenes.png)

Each row and column of the heatmap correspond to a single gene. The heatmap depicts topological overlaps (though can depict adjacencies), with light colors denoting low adjacency (overlap) and darker colors higher adjacency (overlap). In addition, the gene dendrograms and module colors are plotted along the top and left side of the heatmap.

```
#### Construct an eigengene Network dendrogram and heatmap with relationship to dx

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# Isolate dx from the clinical traits
dx = as.data.frame(traitData$Diagnosis);
names(dx) = "dx"
# Add the Controls to existing module eigengenes
MET = orderMEs(cbind(MEs, dx))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                      xLabelsAngle = 90)

# Split the dendrogram and heatmap
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot overwrites the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
```

![miRNA_dendro_heatmap_dx.png](assets/images/miRNA_dendro_heatmap_dx.png)
The dendrogram shows hierarchical clustering of the eigengenes in which the dissimilarity of eigengenes Ei, Ej is given by 1 − cor(Ei, Ej). The heatmap below shows the eigengene adjacency Aij = (1 + cor(Ei, Ej))/2. The purpose of these is to identify meta-modules, which are define as tight custers of modules (e.g., modules with a correlation of eigengenes of at least 0.5).

```
# Classic Multi-Dimensional Scaling plot (MDS)
cmd1 = cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(moduleColors), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

```

![MDS_plot.png](assets/images/MDS_plot.png)

Multidimensional scaling takes a set of dissimilarities and returns a set of points such that the distances between the points are approximately equal to the dissimilarities (ordination). Modules here correspond to “fingers”, where intramodular hubs are in the finger tips.


```
#### Construct adjacency matrix, compute Intramodular Connectivity (IM), explore scale-free topology, create module heatmaps with corresponding ME bargraphs

# define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(datExpr,use="p"))^12
# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# Histogram of k and a scale-free topology plot of adjacency matrix
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
```

![miRNA_kHist.png](assets/images/miRNA_kHist.png)
The left panel shows a histogram of network connectivities. The right panel is a log-log plot of the same
histogram; the approximate linear relationship (high R^2 value) shows approximate scale free topology.

```
## Compute INTRAMODULAR CONNECTIVITY
## Once you've computed your modules and identified modules of interest, we can identify the module hub genes.
# Select a module: "lightcyan", "darkgrey"
module <- "lightcyan";
xname = module
# Select module probes
probes <- names(datExpr)
inModule <- (moduleColors==module);
modProbes <- probes[inModule];
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule];
dimnames(modTOM) <- list(modProbes, modProbes)
# Since miRNA modules are small, top hub genes is limited to 10
nTop = 10;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)

length(modProbes)

# Make connectivity table
IMConn.table <- as.data.frame(cbind(modProbes, IMConn))
rownames(IMConn.table) <- IMConn.table$modProbes
IMConn.table[order(-as.numeric(IMConn.table$IMConn)),]

# Make a histagram of connectivity. Confirm that this follows a scale-free, power law distribution.
hist(as.numeric(IMConn.table$IMConn), breaks = 20, main = paste("Histogram of" , xname),
     xlab = "Intramodular Connectivity")
```

![miRNA_lightcyan_hist.png](assets/images/miRNA_lightcyan_hist.png)
![miRNA_darkgrey_hist.png](assets/images/miRNA_darkgrey_hist.png)

```
topIMConn <- IMConn.table[order(-as.numeric(IMConn.table$IMConn)),][1:30,]
# write.csv(topIMConn, "miRNA_Top_IMConn_lightcyan.csv")


# Adjacency can be used to define a separate measure of similarity, the TOM
dissTOM=TOMdist(ADJ1)
collectGarbage()

# Average linkage hierachical clustering with adjacency-based dissimilarity
# Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1
# use of average linkage hierachical clustering.
hierADJ=hclust(as.dist(dissADJ), method="average")

#  Representing modules by eigengenes and relating eigengenes to one another
datME=moduleEigengenes(datExpr,col=as.character(moduleColors))$eigengenes
signif(cor(datME, use="p"), 2)
#write.csv(datME, "datME_miRNA_all.csv")

# define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
# between the module eigengenes, and use it to cluster the eigengene
dissimME=(1-t(cor(datME_dx, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average")

# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
```

![HierClus_dendo_modEigs.png](assets/images/HierClus_dendo_modEigs.png)

Hierarchical clustering dendrograms of module eigengenes show a summary of each module by its eigengene, i.e., first principal component. The relatedness of the modules is reflected in low merge heights.  

```
# create a pairwise scatter plots of the samples (arrays) along the module eigengenes
sizeGrWindow(8,9)
plotMEpairs(datME_dx)
```
![miRNA_modGene_relationship.png](assets/images/miRNA_modGene_relationship.png)
Above is a pairwise scatterplot of module eigengenes of interest, the distribution of their values, and their pairwise correlation. On the diagonal are histograms of sample values for each eigengene. Above the diagonal is a pairwise scatterplot of the module eigengenes, below the diagonal is a pairwise correlation of the eigengene. The module eigengenes (first PCs) of the two different modules appears to be correlated (r = 0.58). 

```
# Create gene expression heatmaps with corresponding module eigengene expression for significant modules

col = as.character(moduleColors)

sizeGrWindow(9,8)
which.module="lightcyan";
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
plotMat(t(scale(datExpr[,col==which.module ]) ),
        nrgcols=30,rlabels=F,
        clabels=names[-1],rcols=which.module,
        main=which.module,cex.main=2)
par(mar=c(4, 4.2, 0.2, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")

sizeGrWindow(9,8)
which.module="darkgrey";
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
plotMat(t(scale(datExpr[,col==which.module ]) ),
        nrgcols=30,rlabels=F,
        clabels=names[-1],rcols=which.module,
        main=which.module,cex.main=2)
par(mar=c(4, 4.2, 0.2, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")

datME0 = as.data.frame(as.matrix(datME))
```
![miRNA_lightcyan_heatmap_barplot.png](assets/images/miRNA_lightcyan_heatmap_barplot.png)
![miRNA_darkgrey_heatmap_barplot.png](assets/images/miRNA_darkgrey_heatmap_barplot.png)

The top row shows the heatmap of the module genes (rows) across the samples (columns). The lower row shows the corresponding module eigengene (ME) expression values (y-axis) versus the same samples. The ME can be considered the most representative gene expression profile of the module.

```
# A measure of module significance (as the average gene significance in the module) with respect to Dx

GS1=as.numeric(cor(dx,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next, module significance is defined as average gene significance.
ModuleSignificance = tapply(GeneSignificance, moduleColors, mean, na.rm=T)
# Plot module significance
sizeGrWindow(8,7)
plotModuleSignificance(GeneSignificance,moduleColors, cex.names = 0.7, ylim = c(0, 0.5))


collectGarbage()
```

![GS_Dx_across_mods_SignedHybrid.png](assets/images/GS_Dx_across_mods_SignedHybrid.png)

The top row shows the heatmap of the module genes (rows) across the samples (columns). The lower row shows the corresponding module eigengene (ME) expression values (y-axis) versus the same samples. The ME can be considered the most representative gene expression profile of the module.

```
#### Correlate module eigengenes with Dx

dx = traitData$Diagnosis;
datME=moduleEigengenes(datExpr,col=as.character(moduleColors))$eigengenes
signif(cor(dx, datME, use="p"), 2)
#     MEblue MEdarkgreen MEdarkgrey MEdarkorange MEdarkred MEdarkturquoise MEgreenyellow MEgrey
# dx   -0.4        0.33      -0.45          0.3      0.26           -0.19         -0.44  -0.32
#       MElightcyan MEmagenta MEorange MEpink MEred MEroyalblue MEwhite
# dx       -0.52      0.27    -0.12    -0.42    -0.35    -0.028    0.44

# Pearson's product-moment correlation
cor.test(dx, datME$MElightcyan)
# data:  dx and datME$MElightcyan
# t = -2.5742, df = 18, p-value = 0.01911
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.78179228 -0.09891708
# sample estimates:
# cor 
# -0.5187315 

cor.test(dx, datME$MEdarkgrey)
# data:  dx and datME$MEdarkgrey
# t = -2.1346, df = 18, p-value = 0.0468
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.743998392 -0.008653582
# sample estimates:
#   cor 
# -0.4494532


#### Intramodular connectivity, module membership, and screening for intramodular hub genes

# calculate intramodular connectivity for each gene
ADJ1=abs(cor(datExpr,use="p"))^12
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1)
#write.table(Alldegrees1, "miRNA_IM_connectivity_AllGenes.txt", sep="\t")

# Identify gene in each module with the highest connectivity, looking at all genes in the expression file
hubGene_TopPerMod = chooseTopHubInEachModule(datExpr, moduleColors,
                                             type = "unsigned")
# write.csv(hubGene_TopPerMod, "hubGene_TopPerMod_miRNA.csv")

# plot gene significance against intramodular connectivity
colorlevels=unique(moduleColors)
# for MExDx
colorlevels=colorlevels[-c(1,3:9,11:12,15)]
colorlevels
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
```

![miRNA_darkgrey_lightcyan_ConnectivityvGS.png](assets/images/miRNA_darkgrey_lightcyan_ConnectivityvGS.png)
Intramodular connectivity, or "degree", was calculated for each gene; i.e., the whole network connectivity and within-module connectivity were calculated. The above plots show approximately how significant the intramodular hub genes are to the module.



```R

```
