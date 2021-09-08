# eRuption
#
# This is the main function for producging
# eRuption volcano plots
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T

#Michael P Mariani PhD 2019-2021

#Make a customizable volcano plot using gene expression data output
#by DESeq2 or edgeR

library(data.table)
library(ggplot2)
library(ggrepel)

eRuption <- function(input.file = "",
                     plot.title = "",
                     subset.genes = NULL,
                     select.genes = NULL,
                     edge.r = FALSE,
                     output.dir = FALSE,
                     output.name = NULL){

if(edge.r==FALSE){

#################### DESeq2 ###############################################

input.file <- file.in

de.data <- read.table(file=input.file,
                      header=TRUE,
                      stringsAsFactors = FALSE,
                      sep=",")

setDT(de.data)
colnames(de.data)[1] <- "gene"
print(nrow(de.data))

if(!is.null(subset.genes)){
  de.data <- subset(de.data, gene %in% subset.genes)
}

volcano.plot <- ggplot(data=de.data,aes(x=log2FoldChange,y=log10padj)) +
  geom_point(color=ifelse((abs(de.data$log2FoldChange)<1 | de.data$padj > 0.05),
                          "gray",
                          ifelse(de.data$log2FoldChange<0, 'blue', 'red'))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab("log 2 fold change") +
  ylab("-log10(p-value)") +
  ggtitle(title.in) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(0.05), color="green") +
  geom_vline(xintercept = -1, color="green") +
  geom_vline(xintercept = 1, color="green") +
  geom_label_repel(data=subset(de.data, gene %in% select.genes),
                   #data=head(subset(de.data[order(de.data$padj,
                   #decreasing=FALSE)], log2FoldChange>0), 10),
                   aes(label = gene),
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 200) +
  geom_label_repel(data=subset(de.data, gene %in% select.genes),
                   #data=head(subset(de.data[order(de.data$padj,
                   #decreasing=FALSE)], log2FoldChange<0), 10),
                   aes(label = gene),
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 200)

if(save.plot==TRUE){
  ggsave(filename = paste0(output.dir,"\\",output.name),
         device="pdf",
         plot=volcano.plot,
         height=8,
         width=8)
}

return(volcano.plot)

################### edgeR ############################################

}else{

de.data <- read.table(file = input.file,
                      header=TRUE,
                      sep="\t",
                      stringsAsFactors = FALSE)

if(!is.null(subset.genes)){
  de.data <- subset(de.data, gene %in% subset.genes)
}

setDT(de.data)
colnames(de.data)[1] <- "gene"
print(nrow(de.data))

if(any(is.na(de.data))){de.data=na.omit(de.data)}

volcano.plot <- ggplot(data=de.data,aes(x=log2FoldChange,y=-log10(pvalue))) +
  geom_point(color=ifelse((abs(de.data$log2FoldChange)<1 | de.data$padj > 0.05),
                          "gray",ifelse(de.data$log2FoldChange<0, 'blue', 'red'))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab("log 2 fold change") +
  ylab("-log10(p-value)") +
  ggtitle(title.in) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(0.05), color="green") +
  geom_vline(xintercept = -1, color="green") +
  geom_vline(xintercept = 1, color="green") +
  geom_label_repel(##data=subset(de.data, gene %in% select.genes),
                   subset(de.data[order(de.data$padj, decreasing=FALSE)],
                          log2FoldChange>0 & gene %in% select.genes),
                   aes(label = gene),
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 200) +
  geom_label_repel(##data=subset(de.data, gene %in% select.genes),
                   data=subset(de.data[order(de.data$padj, decreasing=FALSE)],
                               log2FoldChange<0 & gene %in% select.genes),
                   aes(label = gene),
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 200)

if(save.plot==TRUE){
  ggsave(filename = paste0(output.dir,"\\",output.name),
         device="pdf",
         plot=volcano.plot,
         height=8,
         width=8)
}

return(volcano.plot)

}

}
