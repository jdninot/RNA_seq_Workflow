################################################################################
############################## RNA-seq workflow ################################
################################################################################

################################################################################
############################## Packages ########################################
################################################################################

library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(corrplot)
library(PoiClaClu)
library(pheatmap)
library(DESeq2)
library(vsn)
library(ggbeeswarm)
library(glmpca)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(IHW)
library(gridExtra)

################################################################################
############################## About data ######################################
################################################################################

# Title: Characterization of the SARS-CoV-2 Host Response in Primary Human
# Airway Epithelial Cells from Aged Individuals

# Overall design: mRNA profiling of primary human bronchial epithelial
# cells (HBECs) from normal, nonsmoker, aged individuals (>67 years old)
# either uninfected or infected with SARS-CoV-2 as a function of time
# post-infection (24, 48, 72 and 96 hr)

# Reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175779

################################################################################
############################## Reading data ####################################
################################################################################

df = read.table("GSE175779.txt", header = TRUE, row.names = 1)
str(df)

################################################################################
############################## Cleaning data ###################################
################################################################################

colnames(df) = gsub("SCOV2", "infected", colnames(df))
colnames(df) = gsub("hr", "", colnames(df))

condsplit = function(x){paste(str_split(x, "_")[[1]][c(3, 2, 4)], collapse = "")}
colnames(df) = unlist(lapply(colnames(df),condsplit))

colnames(df) = gsub("d", "d_", colnames(df))
colnames(df) = gsub("S", "_S", colnames(df))
colnames(df) = gsub("NANANA", "ENSEMBL", colnames(df))

df = df[, sort(colnames(df))]
df = df[, !colnames(df) %in% c("uninfected_0_S1",
                                "uninfected_0_S10",
                                "uninfected_0_S19",
                                "uninfected_0_S28")]

str(df)

################################################################################
############################## Data quality control ############################
################################################################################

QSplot = function(data, nameplot = "Other", namedata = "No info"){
  
  data_mtrx = as.matrix(data)
  data["SYM"] = rownames(data)
  data_melt = melt(data, id = "SYM")
  nb.cols = dim(data_mtrx)[2]
  
  if(nameplot == "boxplot"){
    
    ggplot(data_melt, aes(variable, value, color = variable)) + 
      geom_boxplot() + 
      theme_minimal() + 
      scale_color_manual(values = colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)) +
      theme(axis.text.x = element_text(angle = 90), 
            legend.position = "none", 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      labs(title = paste("Box plot |", namedata), x = "Sample", y = "Count")
    
  }else if(nameplot == "corrplot"){
    
    corrplot(cor(data_mtrx), method = 'color', tl.cex = 0.6, tl.col = "black",
             main = paste("Correlogram |",namedata), mar = c(0, 2, 2, 2))
    
  }else if(nameplot == "heatmap"){
    
    poisd = PoissonDistance(t(data_mtrx))
    samplePoisDistMatrix = as.matrix(poisd$dd)
    rownames(samplePoisDistMatrix) = colnames(data_mtrx)
    colnames(samplePoisDistMatrix) = colnames(data_mtrx)
    
    pheatmap(samplePoisDistMatrix,
             clustering_distance_rows = poisd$dd,
             clustering_distance_cols = poisd$dd,
             col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255)),
             fontsize_row = 7,fontsize_col = 7,
             main = paste("Heatmap |", namedata))
    
  }else if(nameplot == "meansdplot"){
    
    msd = meanSdPlot(data_mtrx, ranks = TRUE, plot = FALSE)
    msd$gg + ggtitle(paste("Plot row standard deviations versus row means |", namedata)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(hjust = 0.5)) 
    
  }else{
    
    infoqs = list(paste("Summary about empty data",namedata),
                  data.frame(Empty = apply(data, 2, function(x){sum(is.na(x)*1)})),
                  paste("Total NA data:", sum(apply(data, 2, function(x){sum(is.na(x)*1)}))),
                  "You can indicate a graph: boxplot, corrplot, heatmap, meansdplot")
    print(infoqs)
    
  }
}

############################## With all the data ###############################

QSplot(df, namedata = "GSE175779")
QSplot(df, nameplot = "boxplot", namedata = "GSE175779")
QSplot(df, nameplot = "corrplot", namedata = "GSE175779")
QSplot(df, nameplot = "heatmap", namedata = "GSE175779")
QSplot(df, nameplot = "meansdplot", namedata = "GSE175779")

############################# Filtering rule ###################################

# keep = rowSums(df >= 10) >= 3 / rowSums(df) > 1
keep = rowSums(df >= 10) >= 3 & rowSums(df) < 4e5
df2 = df[keep, ]

dim(df)
dim(df2)

QSplot(df2, namedata = "GSE175779")
QSplot(df2, nameplot = "boxplot", namedata = "GSE175779")
QSplot(df2, nameplot = "corrplot", namedata = "GSE175779")
QSplot(df2, nameplot = "heatmap", namedata = "GSE175779")
QSplot(df2, nameplot = "meansdplot", namedata = "GSE175779")

######################## Variance filtering rule ###############################

keep = 0 < apply(df2, 1, var)
df3 = df2[keep, ]
 
dim(df)
dim(df2)
dim(df3)

QSplot(df3, namedata = "GSE175779")
QSplot(df3, nameplot = "boxplot", namedata = "GSE175779")
QSplot(df3, nameplot = "corrplot", namedata = "GSE175779")
QSplot(df3, nameplot = "heatmap", namedata = "GSE175779")
QSplot(df3, nameplot = "meansdplot", namedata = "GSE175779")

####################### Filter samples with low correlation ####################

df4 = df3[, apply(cor(as.matrix(df3)), 1, min) > 0.70]

dim(df)
dim(df2)
dim(df3)
dim(df4)

QSplot(df4, namedata = "GSE175779")
QSplot(df4, nameplot = "boxplot", namedata = "GSE175779")
QSplot(df4, nameplot = "corrplot", namedata = "GSE175779")
QSplot(df4, nameplot = "heatmap", namedata = "GSE175779")
QSplot(df4, nameplot = "meansdplot", namedata = "GSE175779")

################################################################################
####################### Differential expression analysis #######################
################################################################################

####################### Table of conditions and factors ########################

countdata = as.matrix(df4)

condsplit_condition = function(x){paste(str_split(x, "_")[[1]][1], collapse = "")}
condsplit_time = function(x){paste(str_split(x, "_")[[1]][2], collapse = "")}

condition =  factor(unlist(lapply(colnames(countdata),condsplit_condition)))
time =  factor(unlist(lapply(colnames(countdata),condsplit_time)))
coldata = data.frame(condition, time)
row.names(coldata) = colnames(countdata)

coldata$condition = relevel(coldata$condition, ref = 'uninfected')

############################### DESeq Matrix ###################################

deseq_mtrx = DESeqDataSetFromMatrix(countData = countdata,
                                    colData = coldata,
                                    design = ~ time*condition)

deseq_mtrx

######################### Summary DESEq 2 Model ################################

# The estimation of size factors (controlling for differences in the sequencing 
# depth of the samples), the estimation of dispersion values for each gene, and 
# fitting a generalized linear model.

deseq_model = DESeq(deseq_mtrx, test="LRT", reduced = ~ condition + time)
deseq_reslt = results(deseq_model)
mcols(deseq_reslt, use.names = TRUE)
summary(deseq_reslt)
table(deseq_reslt$padj < 0.1)

deseq_reslt
resultsNames(deseq_model)

############################## Plots summary ###################################

DESeq_plot = function(deseqmodel, dqsresfcS, nameplot = "Other", namedata = "No info", ylim = c(-1, 1)){
  
  datalfc = data.frame(
    baseMean = dqsresfcS[,1],
    log2FoldChange = dqsresfcS[,2],    
    lfcSE = dqsresfcS[,3],     
    pvalue = dqsresfcS[,4],      
    padj = dqsresfcS[,5],
    names = rownames(dqsresfcS))
  
  datalfc$expression = ifelse(is.na(datalfc$padj), FALSE, datalfc$padj < 0.1)
  datalfc$expression[datalfc$expression == TRUE & datalfc$log2FoldChange < 0] = "Low"
  datalfc$expression[datalfc$expression == TRUE & datalfc$log2FoldChange > 0] = "High"
  datalfc$expression[datalfc$expression == FALSE] = "Neutral"

  if(nameplot == "Sizefactorsplot"){
    
    sF = data.frame(Value = sizeFactors(deseqmodel))
    sF["Sample"] = row.names(sF)
    ggplot(sF, aes(x = Sample, y = Value)) +
      geom_point(aes(size=Value), color = "gray") +
      geom_hline(yintercept=mean(t(sF['Value'])), linetype='dashed') +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 90), 
            legend.position = "none", 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      labs(title = paste("Size Factors plot |", namedata), x = "Variable", y = "Size Factors")
    
  }else if(nameplot == "pvhisplot"){
    
    ggplot(datalfc, aes(x=padj)) + geom_histogram(fill = "#2166AC") +
      labs(y='Frequency') + scale_x_continuous(breaks = seq(0,0.95,0.05)) +
      theme_bw() +
      geom_vline(aes(xintercept=0.1), color = "#B2182B",
                 size = 0.7, show.legend = T) +
      theme(axis.text.x = element_text(angle = 90)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = 'none',
            plot.title = element_text(hjust = 0.5)) +
      xlab("P-adjusted") +
      ggtitle(paste("P-adjusted values |", namedata))
    
  }else if(nameplot == "Dispersionplot"){
    
    plotDispEsts(deseqmodel, main = paste("Dispersion plot |", namedata)) 
    
  }else if(nameplot == "TopGeneplot"){
    
    topGene = rownames(dqsresfcS)[which.min(dqsresfcS$padj)]
    geneCounts = plotCounts(deseqmodel, gene = topGene, returnData = TRUE)
    ggplot(geneCounts, aes(x = condition, y = count, color = condition)) +
      scale_y_log10() +
      geom_boxplot() + 
      geom_beeswarm(cex = 3) +
      theme_minimal() +
      scale_color_manual(values = colorRampPalette(brewer.pal(8, "RdBu"))(2)) +
      theme(axis.text.x = element_text(angle = 90), 
            legend.position = "none", 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      labs(title = paste(topGene, "|", namedata), x = "Condition", y = "Count")
    
  }else if(nameplot == "PCAplot"){
    
    vsd = vst(deseqmodel, blind = FALSE)
    pcaData = plotPCA(vsd, intgroup = c("condition","time"), returnData = TRUE)
    percentVar = round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = time)) +
      geom_point(size =3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      theme_minimal() +
      ggtitle("PCA with VST data") +
      theme(axis.text.x = element_text(angle = 90), 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = colorRampPalette(brewer.pal(8, "RdBu"))(2)) +
      labs(title = paste("PCA |", namedata))
    
  }else if(nameplot == "GMLPCAplot"){
    
    gpca = glmpca(counts(deseqmodel), L=2)
    gpca.dat = gpca$factors
    gpca.dat$condition = deseq_model$condition
    gpca.dat$time = deseq_model$time
    ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition, shape = time)) +
      geom_point(size =3) + coord_fixed() + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90), 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = colorRampPalette(brewer.pal(8, "RdBu"))(2)) +
      labs(title = paste("GML PCA |", namedata))
    
  }else if(nameplot == "Cookplot"){
    
    nb.cols = dim(log10(assays(deseqmodel)[["cooks"]]))[2]
    dist_cook = melt(log10(assays(deseqmodel)[["cooks"]]))
    colnames(dist_cook)
    ggplot(dist_cook, aes(Var2, value, color = Var2)) + 
      geom_boxplot() + 
      theme_minimal() + 
      scale_color_manual(values = colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)) +
      theme(axis.text.x = element_text(angle = 90), 
            legend.position = "none", 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      labs(title = paste("Cook plot |", namedata), x = "Sample", y = "Count")
    
  }
}

DESeq_plot(deseq_model, as.data.frame(deseq_reslt), nameplot = "Sizefactorsplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_reslt), nameplot = "pvhisplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_reslt), nameplot = "Dispersionplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_reslt), nameplot = "Cookplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_reslt), nameplot = "PCAplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_reslt), nameplot = "GMLPCAplot", namedata = "GSE175779")

################################## Time plot ###################################

plots = list()
p_list = head(order(deseq_reslt$padj), 6)
for (i in p_list) {
  df_counts = plotCounts(deseq_model, i, intgroup = c("time","condition"), returnData = TRUE)
  name_gen = row.names(as.data.frame(deseq_reslt)[i,])
  p = ggplot(df_counts, aes(x = time, y = count, color = condition, group = condition)) + 
    geom_point() + 
    theme_minimal() +
    scale_y_log10() +
    stat_summary(fun.y=mean, geom="line") +
    scale_color_manual(values = colorRampPalette(brewer.pal(8, "RdBu"))(2)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), 
          plot.title = element_text(hjust = 0.5)) +
    labs(title = paste("Gene expression plot |", name_gen, "| GSE175779"), x = "Time", y = "Count")
  plots[[i]] = p
}
grid.arrange(grobs=plots[p_list],ncol=2)
resultsNames(deseq_model)

################################## Condition ###################################

deseq_reslt_1 = results(deseq_model, name = "time_72_vs_24", test="Wald")
summary(deseq_reslt_1)
deseq_reslt_1[which.min(deseq_reslt_1$padj),]

deseq_reslt_2 = results(deseq_model, name = "time48.conditioninfected", test="Wald")
summary(deseq_reslt_2)
deseq_reslt_2[which.min(deseq_reslt_2$padj),]

deseq_reslt_3 = results(deseq_model, name = "time72.conditioninfected", test="Wald")
summary(deseq_reslt_3)
deseq_reslt_3[which.min(deseq_reslt_3$padj),]

deseq_reslt_4 = results(deseq_model, name = "time96.conditioninfected", test="Wald")
summary(deseq_reslt_4)
deseq_reslt_4[which.min(deseq_reslt_4$padj),]


betas = coef(deseq_model)
colnames(betas)

topGenes = head(order(deseq_reslt$padj),20)
mat = betas[topGenes, -c(1,2)]
thr = 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, 
         main = paste("Heatmap of log2 fold changes for genes \n with smallest adjusted p value |", "GSE175779"),
         col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255)),
         fontsize_row = 7,fontsize_col = 7)

sessionInfo()

# More details
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html