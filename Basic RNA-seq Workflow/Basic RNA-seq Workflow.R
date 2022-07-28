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

######################### Minimal filtering rule ###############################

keep = rowSums(df) > 1
keep = rowSums(df >= 10) >= 3
df2 = df[keep, ]

dim(df)
dim(df2)

QSplot(df2, namedata = "GSE175779")
QSplot(df2, nameplot = "boxplot", namedata = "GSE175779")
QSplot(df2, nameplot = "corrplot", namedata = "GSE175779")
QSplot(df2, nameplot = "heatmap", namedata = "GSE175779")
QSplot(df2, nameplot = "meansdplot", namedata = "GSE175779")

####################### Filter samples with low correlation ####################

# df2 = df[, apply(cor(df_mtrx), 1, min) > 0.85]
df3 = df2[, !colnames(df) %in% c("uninfected_0_S1",
                                 "uninfected_0_S10",
                                 "uninfected_0_S19",
                                 "uninfected_0_S28",
                                 "uninfected_96_S14", 
                                 "uninfected_96_S23",
                                 "infected_96_S18",
                                 "infected_96_S27")]

dim(df)
dim(df2)
dim(df3)

QSplot(df3, namedata = "GSE175779")
QSplot(df3, nameplot = "boxplot", namedata = "GSE175779")
QSplot(df3, nameplot = "corrplot", namedata = "GSE175779")
QSplot(df3, nameplot = "heatmap", namedata = "GSE175779")
QSplot(df3, nameplot = "meansdplot", namedata = "GSE175779")

################################################################################
####################### Differential expression analysis #######################
################################################################################

####################### Table of conditions and factors ########################

countdata = as.matrix(df3)

condsplit_condition = function(x){paste(str_split(x, "_")[[1]][1], collapse = "")}
condsplit_time = function(x){paste(str_split(x, "_")[[1]][2], collapse = "")}

condition =  factor(unlist(lapply(colnames(countdata),condsplit_condition)))
time =  factor(unlist(lapply(colnames(countdata),condsplit_time)))
coldata = data.frame(condition, time)
row.names(coldata) = colnames(countdata)

coldata$condition = relevel(coldata$condition, ref = 'uninfected')
str(coldata)

deseq_mtrx = DESeqDataSetFromMatrix(countData = countdata,
                                    colData = coldata,
                                    design = ~ condition)

deseq_mtrx

######################### Summary DESEq 2 Model ################################

# The estimation of size factors (controlling for differences in the sequencing 
# depth of the samples), the estimation of dispersion values for each gene, and 
# fitting a generalized linear model.

deseq_model = DESeq(deseq_mtrx, test="LRT", reduced=~1)
deseq_reslt = results(deseq_model)

summary(deseq_reslt)

deseq_reslt
mcols(deseq_reslt, use.names = TRUE)

######################## Summary DESEq 2 Model with Shrink #####################

resultsNames(deseq_model)
deseq_lfcS = lfcShrink(deseq_model, 
                       coef="condition_infected_vs_uninfected", 
                       type="apeglm")

summary(deseq_lfcS)

deseq_lfcS
mcols(deseq_lfcS, use.names = TRUE)

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
  LFCy = datalfc$log2FoldChange
  resSig = subset(datalfc, expression!="Neutral")
  
  if(nameplot == "MAplot"){
    
    ggplot(datalfc, aes(x = baseMean, y = pmax(ylim[1], pmin(ylim[2], log2FoldChange)),
                        col=as.factor(expression))) + 
      geom_point(alpha = 0.7,
                 shape = ifelse(LFCy < ylim[1], 6, ifelse(LFCy > ylim[2], 2, 16)),
                 size  = ifelse(LFCy < ylim[1], 2, ifelse(LFCy > ylim[2], 2, 1.5))) + 
      theme_minimal() +
      scale_x_log10() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
      scale_color_manual(values = c('#B2182B', '#2166AC', "gray65")) +
      labs(title = paste("MA plot |", namedata), x = "Mean of normalized counts", y = "Log fold change") + 
      theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = 'none',
            plot.title = element_text(hjust = 0.5)) 
    
  }else if(nameplot == "Volcanoplot"){
    
    ggplot(datalfc, aes(x = pmax(ylim[1], pmin(ylim[2], log2FoldChange)), y = -log10(pvalue),
                        col=as.factor(expression))) +
      geom_point(alpha = 0.7,
                 shape = ifelse(LFCy < ylim[1], 10, ifelse(LFCy > ylim[2], 10, 16)),
                 size  = ifelse(LFCy < ylim[1], 2, ifelse(LFCy > ylim[2], 2, 1.5))) + 
      theme_minimal() +
      scale_color_manual(values = c('#b2182b', '#2166ac', "#979A9A")) +
      geom_hline(yintercept= -log(0.1), linetype = "dashed", col="darkgrey") +
      geom_vline(xintercept= max(resSig$log2FoldChange[resSig$log2FoldChange<0]),
                 linetype = "dashed", col="darkgrey") +
      geom_vline(xintercept= min(resSig$log2FoldChange[resSig$log2FoldChange>0]),
                 linetype = "dashed", col="darkgrey") +
      labs(title = paste("Volcano plot |", namedata), x = "log Fold Change", y = "-log10(p-value)") + 
      theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = 'none',
            plot.title = element_text(hjust = 0.5))
    
  }else if(nameplot == "Sizefactorsplot"){
    
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
    
    ggplot(datalfc, aes(x=padj)) + geom_histogram() +
      labs(y='Frequency') + # scale_x_continuous(breaks = seq(0,0.95,0.05)) +
      theme_bw() +
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
    plotPCA(vsd) + theme_minimal() +
      theme(axis.text.x = element_text(angle = 90), 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = colorRampPalette(brewer.pal(8, "RdBu"))(2)) +
      labs(title = paste("PCA |", namedata))
    
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
    
  }else if(nameplot == "GMLPCAplot"){
    
    gpca = glmpca(counts(deseqmodel), L=2)
    gpca.dat = gpca$factors
    gpca.dat$condition = deseqmodel$condition
    ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition)) +
      geom_point(size =3) + coord_fixed() + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90), 
            plot.margin = unit(c(1,1,1,1), "cm"), 
            plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = colorRampPalette(brewer.pal(8, "RdBu"))(2)) +
      labs(title = paste("GML PCA |", namedata))

  }
}

DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "MAplot", namedata = "GSE175779")  
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "Volcanoplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "Sizefactorsplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "pvhisplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "Dispersionplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "Cookplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "TopGeneplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "PCAplot", namedata = "GSE175779")
DESeq_plot(deseq_model, as.data.frame(deseq_lfcS), nameplot = "GMLPCAplot", namedata = "GSE175779")

################################## Info Gene ###################################

resSig = subset(deseq_lfcS, padj < 0.1)

head(resSig[order(resSig$log2FoldChange), ]) 
head(resSig[order(resSig$log2FoldChange, decreasing = TRUE), ])

############################### Cluster Gene ###################################

pheatmap(counts(deseq_model)[row.names(resSig),][1:20,],
         cluster_cols=FALSE,
         col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255)),
         fontsize_row = 7,fontsize_col = 7,
         main = paste("Heatmap TopGene |", "GSE175779"))

vsd = vst(deseq_model, blind = FALSE)
select = order(rowMeans(counts(deseq_model)),decreasing=TRUE)[1:20]
df_cluster = as.data.frame(colData(deseq_model)[,c("condition","time")])
pheatmap(assay(vsd)[select,],
          cluster_rows=TRUE, show_rownames=TRUE,
          cluster_cols=FALSE, 
          col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255)),
          fontsize_row = 7,fontsize_col = 7,
          main = paste("Heatmap TopGene |", "GSE175779"))

################################################################################
########################## Summary info and Annotation #########################
################################################################################

DEG_GEN = as.data.frame(resSig)
DEG_GEN['ENSEMBL'] = rownames(DEG_GEN)
edb = EnsDb.Hsapiens.v86
annotationdata = as.data.frame(select(edb, 
                                      keys = row.names(DEG_GEN),
                                      column = "SYMBOL",
                                      keytype = "GENEID"))
DEG_GEN_INFO = DEG_GEN %>% left_join(annotationdata, by = c("ENSEMBL" = "GENEID"))
apply(DEG_GEN_INFO,2,function(x){sum(is.na(x)*1)})

head(DEG_GEN_INFO)

##################################### Save #####################################

write.csv(DEG_GEN_INFO, file = "DEG_Deseq.csv", row.names = FALSE)

sessionInfo()

citation("tidyverse")
citation("reshape2")
citation("RColorBrewer")
citation("corrplot")
citation("PoiClaClu")
citation("pheatmap")
citation("DESeq2")
citation("vsn")
citation("ggbeeswarm")
citation("glmpca")
citation("biomaRt")
citation("EnsDb.Hsapiens.v86")


