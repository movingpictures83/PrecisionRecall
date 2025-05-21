#Library R packages
library(RColorBrewer)
library(ggplot2)
library(ggsignif)
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
source("RIO.R")
library(reshape2)
library(dplyr)

input <- function(inputfile) {
   parameters <<- readParameters(inputfile)
}

run <- function() {}

output <- function(outputfile) {
	pfix <- prefix()
#Read plot data
Plot.Fig3 <-read.table(file = paste(pfix, parameters["metadata", 2], sep="/"),header = T,sep="\t", row.names=1)

#Plot
ggplot(data=Plot.Fig3, aes(x=Profilers, y=Distance, fill=Profilers))+
  stat_boxplot(geom = "errorbar", colour="black", width=0.5)+
  geom_point(aes(x=Profilers, y=Distance, color=Colour))+
  geom_boxplot() + theme_bw() + facet_grid(Matrix~Compared_with, scales = "free")+
  geom_signif(comparisons = list(c("Kraken2","mOTUs2"), c("Kraken2","MPA2"),
                                 c("Bracken","mOTUs2"), c("Bracken","MPA2")),
              step_increase = 0.1, map_signif_level = T, test = wilcox.test)+
  scale_fill_manual(values = brewer.pal(7, "Set1"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
 #library R packages

    auprc <- read.table(paste(pfix, parameters["auprc", 2], sep="/"),header = T, sep="\t")

 #plotting
  
    ggplot(data=auprc, aes(x=profilers, y=Precision, fill=profilers))+
      stat_boxplot(geom = "errorbar",width=0.15)+
      geom_boxplot() + theme_bw() +
      scale_fill_manual(values = brewer.pal(8, "Set1")) +
      geom_signif(comparisons = list(c("Bracken","mOTUs2"),c("Bracken","MPA2"),c("Kraken2","mOTUs2"),c("Kraken2","MPA2")),step_increase = 0.1,map_signif_level = T,test = wilcox.test)

    ggplot(data=auprc, aes(x=profilers, y=Recall, fill=profilers))+
      stat_boxplot(geom = "errorbar",width=0.15)+
      geom_boxplot() + theme_bw() +
      scale_fill_manual(values = brewer.pal(8, "Set1")) +
      geom_signif(comparisons = list(c("Bracken","mOTUs2"),c("Bracken","MPA2"),c("Kraken2","mOTUs2"),c("Kraken2","MPA2")),step_increase = 0.1,map_signif_level = T,test = wilcox.test)

    ggplot(data=auprc, aes(x=profilers, y=F1_score, fill=profilers))+
      stat_boxplot(geom = "errorbar",width=0.15)+
      geom_boxplot() + theme_bw() +
      scale_fill_manual(values = brewer.pal(8, "Set1")) +
      geom_signif(comparisons = list(c("Bracken","mOTUs2"),c("Bracken","MPA2"),c("Kraken2","mOTUs2"),c("Kraken2","MPA2")),step_increase = 0.1,map_signif_level = T,test = wilcox.test)
    
    ggplot(data=auprc, aes(x=profilers, y=AUPRC, fill=profilers))+
      stat_boxplot(geom = "errorbar",width=0.15)+
      geom_boxplot() + theme_bw() +
      scale_fill_manual(values = brewer.pal(8, "Set1")) +
      geom_signif(comparisons = list(c("Bracken","mOTUs2"),c("Bracken","MPA2"),c("Kraken2","mOTUs2"),c("Kraken2","MPA2")),step_increase = 0.1,map_signif_level = T,test = wilcox.test)
 
 # For the plotting script for Fig.S3 d-f, please refer to the FigS3d_f.ipynb#'-------------------------------
# install and load necessary libraries for data analyses
#-------------------------------
#p <- c("reshape2","ggplot2", "reshape2", "dplyr")
#usePackage <- function(p) {
#  if (!is.element(p, installed.packages()[,1]))
#    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
#  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
#}
#invisible(lapply(p, usePackage))



d<-read.table(paste(pfix, parameters["rawanalysis", 2], sep="/"), header=TRUE, sep='\t')

d_m<-melt(d)
write.csv(d_m$value, paste(outputfile, "csv", sep="."))

tmp1 <- do.call(rbind, strsplit(as.character(d_m[,1]), '_'))
tmp2 <- do.call(rbind, strsplit(as.character(d_m[,2]), '_'))

d_u<-data.frame(tmp1, tmp2, d_m)
head(d_u)
colnames(d_u)<-c("SampleID", "Profiler", "Type", "Metric", "AbundanceThreshold", "X", "variable", "value")
tail(d_u)

summary(d_u)

precision_plot <- d_u %>% filter(Metric=='precision') %>% 
ggplot(aes(x=AbundanceThreshold, y=value)) + 
geom_point(aes(color=Profiler)) + geom_smooth(aes(group=Profiler, color=Profiler)) + 
#coord_flip()+
ylab("Precision")+
facet_wrap(~Type) + theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename=paste(outputfile,"precision","pdf",sep="."), plot=precision_plot, width=30, height=4)

recall_plot <- d_u %>% filter(Metric=='recall') %>% 
ggplot(aes(x=AbundanceThreshold, y=value)) + 
geom_point(aes(color=Profiler)) + geom_smooth(aes(group=Profiler, color=Profiler)) + 
ylab("Recall")+
#coord_flip()+
facet_wrap(~Type)+ theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename=paste(outputfile,"recallthreshold","pdf",sep="."), plot=recall_plot, width=30, height=4)

f1_plot <- d_u %>% filter(Metric=='f1') %>% 
ggplot(aes(x=AbundanceThreshold, y=value)) + 
geom_point(aes(color=Profiler)) + geom_smooth(aes(group=Profiler, color=Profiler)) + 
ylab("F1 score")+
#coord_flip()+
facet_wrap(~Type)+ theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f1_plot
ggsave(filename=paste(outputfile,"F1threshold","pdf",sep="."), plot=f1_plot, width=30, height=2)

}

