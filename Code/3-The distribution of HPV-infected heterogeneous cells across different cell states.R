##HPV结果的配色修改
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/细胞stat比例")
hpv_proportions<-read.table("stat_percent训练.txt",header=T)
library(tidyr)
hpv_proportions$sample<-paste(hpv_proportions$cell_type,hpv_proportions$State,sep="_")

long_data <- pivot_longer(hpv_proportions[,5:7], cols = c(prop_minus, prop_plus), names_to = "HPV", values_to = "percent")
head(long_data)
long_data$HPV[which(long_data$HPV=="prop_plus")]<-"HPV+"
long_data$HPV[which(long_data$HPV=="prop_minus")]<-"HPV-"
long_data$sample<-factor(long_data$sample)
colors<-c(
  "HPV+" = "#F8807B",
  "HPV-" = "#FFCB4E"
)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot)
library(ggpubr);library(beeswarm);library(presto);
library(GSEABase)
mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1.2)),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = rel(1.5), color = 'black'),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.5),color = 'black'),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5))
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/细胞stat比例/补充")
p = ggplot(long_data, aes( x = sample, weight = percent, fill = HPV))+
  geom_bar( position = "stack")+scale_fill_manual( values = colors)+labs(y = "Percent")+  theme_classic()+coord_flip()
p
ggsave("sample_percent.pdf",width=4,height=8)

p = ggplot(long_data, aes( x = sample, weight = percent, fill = HPV))+
  geom_bar( position = "stack")+scale_fill_manual( values = colors)+labs(y = "Percent")+  mytheme
p
ggsave("sample_percent_x.pdf",width=17,height=4.5)



setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/细胞stat比例/验证集")
hpv_proportions<-read.table("stat_percent.txt",header=T)
hpv_proportions$sample<-paste(hpv_proportions$cell_type,hpv_proportions$State,sep="_")

long_data <- pivot_longer(hpv_proportions[,5:7], cols = c(prop_minus, prop_plus), names_to = "HPV", values_to = "percent")
head(long_data)
long_data$HPV[which(long_data$HPV=="prop_plus")]<-"HPV+"
long_data$HPV[which(long_data$HPV=="prop_minus")]<-"HPV-"
long_data$sample<-factor(long_data$sample)
library(tidyr)


colors<-c(
  "HPV+" = "#F8807B",
  "HPV-" = "#FFCB4E"
)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot)
library(ggpubr);library(beeswarm);library(presto);
library(GSEABase)
mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1.2)),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = rel(1.5), color = 'black'),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.5),color = 'black'),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5))
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/细胞stat比例/补充/验证")
p = ggplot(long_data, aes( x = sample, weight = percent, fill = HPV))+
  geom_bar( position = "stack")+scale_fill_manual( values = colors)+labs(y = "Percent")+  theme_classic()+coord_flip()
p
ggsave("sample_percent.pdf",width=4,height=8)

p = ggplot(long_data, aes( x = sample, weight = percent, fill = HPV))+
  geom_bar( position = "stack")+scale_fill_manual( values = colors)+labs(y = "Percent")+  mytheme
p
ggsave("sample_percent_x.pdf",width=17,height=4.5)

