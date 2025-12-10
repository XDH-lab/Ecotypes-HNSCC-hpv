##immune
rm(list = ls())
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/富集/富集结果7")
road='/lkn_lab/znh/12月课题实验/实验/训练集/结果7/'
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 6, spaceLegend = 0.8) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
cell_types <- list.files(path = "/lkn_lab/znh/12月课题实验/实验/验证集/bulk/TCGA/exp", full.names = FALSE)
cell_types<-cell_types[-which(cell_types=="Ecotypes")]

cells<-cell_types

#immune
immuneGene<- readLines("/lkn_lab/znh/12月课题实验/实验/训练集/富集/Immune.gmt")
resGene <- strsplit(immuneGene, "\t")
names(resGene) <- vapply(resGene, function(y) y[1], character(1))
resGene <- lapply(resGene, "[", -c(1:2))
ResMart<-c()
for (i in 1:length(cells)) {
  gene_info<-read.table(paste0(road,cells[i],'/gene_info.txt'),header = T,as.is = T,sep = '\t')
  states<-unique(gene_info$State)
  ResMartMid3<-c()
  for (j in 1:length(states)) {
    states_gene<-gene_info[which(gene_info$State==states[j]),]$Gene
    ResMartMid2<-c()
    for (s in 1:17) {
      HallMarkerGenes<-resGene[[s]]
      n <- states_gene 
      n_num <- length(n)
      m <- HallMarkerGenes 
      m_num <- length(m)
      N_num <- 21876
      k <- intersect(n,m) 
      k_num <- length(k)
      pvalues <- 1-phyper(k_num-1, m_num, N_num-m_num, n_num)
      
      ResMartMid<-c(cells[i],states[j],names(resGene)[[s]],k_num,pvalues,paste(k,collapse=","))
      ResMart<-rbind(ResMart,ResMartMid)
      
    }
    
  }
  
}
colnames(ResMart)<-c('CellType','State','immunePathway','NumEnrich','pValue','Genes')
dim(ResMart)

###
p<-p.adjust(ResMart[,5],n=nrow(ResMart),method = "fdr")###bonferroni
index1<-which(ResMart[,5]<0.05)
length(index1)
index2<-which(p<0.05)
length(index2)

ResMart_adjust<-ResMart[index2,]
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/富集/另外图片")
write.table(ResMart_adjust,"ResMart_adjust_immune_adjust.txt",quote = F,sep = "\t",row.names = F)



rm(list = ls())
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/富集/另外图片")
#brand设置成cell和stat，通路和type d的合并
library(dplyr)
library(tidyr)
states_pathway_immue<-read.csv("/lkn_lab/znh/12月课题实验/实验/训练集/富集/另外图片/ResMart_adjust_immune_adjust.txt",stringsAsFactors = F,header = T,sep = "\t")
col_stat<-c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")
celltype_colors  <- c(
  "Epithelial_cells" = "#699ECA",
  "Dendritic_cells" = "#FF8C00",
  "Fibroblasts" = "#F898CB",
  "Mono/Mac" = "#4DAF4A",
  "CD8_T_cells" = "#D65190",
  "Neutrophils" = "#731A73",
  "CD4_T_cells" = "#FFCB5B",
  "Endothelial_cells" = "#E87B1E",
  "B_cells" = "#0076B9",
  "NK_cells" = "#3D506A",
  "Mast_cells" = "#0099B2",
  "Plasma_cells" = "#F7B799"
)
col_path<-"#E7298A"


# 构造 alluvium 格式
library(dplyr)
library(tidyr)

sankey_data <- states_pathway_immue %>%
  mutate(flow = as.factor(paste(CellType, State, immunePathway, sep = "_"))) %>%
  pivot_longer(cols = c(CellType, State, immunePathway), names_to = "variable", values_to = "value") %>%
  mutate(
    link = case_when(
      variable == "CellType" ~ 1,
      variable == "State" ~ 2,
      variable == "immunePathway" ~ 3
    ),
    variable = factor(variable, levels = c("CellType", "State", "immunePathway"))
  )
# 构造颜色映射
col_sankey <- c(
  celltype_colors[unique(sankey_data$value[sankey_data$variable == "CellType"])],
  col_stat[unique(sankey_data$value[sankey_data$variable == "State"])],
  setNames(rep(col_path, length(unique(sankey_data$value[sankey_data$variable == "immunePathway"]))),
           unique(sankey_data$value[sankey_data$variable == "immunePathway"]))
)
ggplot(sankey_data,
       aes(x = link, y = 1, stratum = value, alluvium = flow, fill = value)) +
  geom_stratum(width = 0.3, color = "white") +
  geom_flow(aes.flow = "forward", width = 0.3, color = "white") +
  scale_fill_manual(values = col_sankey) +
  geom_text(stat = "stratum", aes(label = value), size = 2.5, color = "black") +
  scale_x_continuous(breaks = 1:3, labels = c("CellType", "State", "immunePathway")) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
ggsave("immue.pdf",height=5,width=6)








