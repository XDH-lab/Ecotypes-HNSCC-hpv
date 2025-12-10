rm(list = ls())


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

#Hallmarker
immuneGene<- readLines("/lkn_lab/znh/12月课题实验/实验/训练集/富集/Hallmarker.gmt")
resGene <- strsplit(immuneGene, "\t")
names(resGene) <- vapply(resGene, function(y) y[1], character(1))
resGene <- lapply(resGene, "[", -c(1:2))
ResMart<-c()
#setwd("F:\\2study\\CRC_study\\data\\real_data\\细胞状态富集分析Circos\\Hallmarker")
for (i in 1:length(cells)) {
  gene_info<-read.table(paste0(road,cells[i],'/gene_info.txt'),header = T,as.is = T,sep = '\t')
  states<-unique(gene_info$State)
  ResMartMid3<-c()
  for (j in 1:length(states)) {
    states_gene<-gene_info[which(gene_info$State==states[j]),]$Gene
    ResMartMid2<-c()
    for (s in 1:50) {
      HallMarkerGenes<-resGene[[s]]
      n <- states_gene 
      n_num <- length(n)
      m <- HallMarkerGenes 
      m_num <- length(m)
      N_num <- 21876
      k <- intersect(n,m) 
      k_num <- length(k)
      pvalues <- 1-phyper(k_num-1, m_num, N_num-m_num, n_num)
      #输出
      ResMartMid<-c(cells[i],states[j],names(resGene)[[s]],k_num,pvalues,paste(k,collapse=","))
      ResMart<-rbind(ResMart,ResMartMid)
      
    }
    
  }
  
}
dim(ResMart)
colnames(ResMart)<-c('CellType','State','HallMarkerPathway','NumEnrich','pValue','Genes')
ResMart_lower<-as.data.frame(ResMart)
ResMart_lower$HallMarkerPathway<-tolower(ResMart_lower$HallMarkerPathway)


###
p<-p.adjust(ResMart_lower[,5],n=nrow(ResMart_lower),method = "fdr")###bonferroni
index1<-which(ResMart_lower[,5]<0.05)##
index2<-which(p<0.05)##
ResMart_adjust<-ResMart_lower[index2,]
length(index2)
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/富集/另外整理")
write.table(ResMart_adjust,"ResMart_adjust_HallMarker_lower.txt",quote = F,sep = "\t",row.names = F)


rm(list = ls())
#brand设置成cell和stat，通路和type d的合并
library(dplyr)
library(tidyr)
#states_pathway_immue<-read.csv("ResMart_adjust_immune_adjust.txt",stringsAsFactors = F,header = T,sep = "\t")
states_pathway_HallMarker<-read.csv("/lkn_lab/znh/12月课题实验/实验/训练集/富集/另外整理/ResMart_adjust_HallMarker_lower.txt",stringsAsFactors = F,header = T,sep = "\t")
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/富集/另外整理")

#head(states_pathway_immue)
head(states_pathway_HallMarker)
states_pathway_HallMarker <- states_pathway_HallMarker %>%
  mutate(
    HallMarkerPathway = str_remove(HallMarkerPathway, "^hallmark_"),
    HallMarkerPathway = str_to_title(str_replace_all(HallMarkerPathway, "_", " "))
  )
states_pathway_HallMarker1<-states_pathway_HallMarker
colnames(states_pathway_HallMarker1)[3]<-"Pathway"

combined_data<-states_pathway_HallMarker
head(combined_data)

hallmarker_long <- states_pathway_HallMarker %>% 
  pivot_longer(cols = -c(CellType, State, NumEnrich, pValue, Genes), 
               names_to = "Pathway", 
               values_to = "Type") %>%
  mutate(Type = "HallMarker")
long_data<-hallmarker_long
head(long_data)
circos_data <- long_data %>% 
  group_by(CellType, State, Pathway) %>% 
  summarise(NumEnrich = sum(NumEnrich), .groups = 'drop')
head(circos_data)

circos_inter<-combined_data[,c(1,3,4)]
colnames(circos_inter)<-c("from","to","value")
head(circos_inter)


#细胞细分到stat，通路分2类
combined_data$cell<-paste(combined_data[,1],combined_data[,2])
circos_inter<-circos_inter<-combined_data[,c(7,3,1)]
colnames(circos_inter)<-c("from","to","value")
head(circos_inter)
pathway_type<-combined_data[,c(7,1)]#用这个

colnames(pathway_type)<-c('Node','Type')
##考虑合并
head(combined_data[,c(3,1)])
colnames(pathway_type)[1]<-"Pathway"
pathway_names <- states_pathway_HallMarker1$Pathway
pathway_types <- rep("HallMarker", length(states_pathway_HallMarker$HallMarkerPathway))
pathway_types <- data.frame(Node = pathway_names, Type = pathway_types)
head(pathway_types)

##细胞对应通路+通路类型对应通路，12+2设置颜色
colnames(pathway_type)<-c('Node','Type')
pathway_type<-rbind(pathway_type,pathway_types)
head(pathway_type)
table(pathway_type[,2])

brand = c(structure(pathway_type$Type, names=pathway_type$Node))
length(brand)
brand = brand[!duplicated(names(brand))]#一个通路可能对应不同细胞的stat去重复有问题吧，后面用table(brand)就不应该去重复
brand = brand[order(brand, names(brand))]
brand=factor(brand,levels = unique(brand)[c(1:7,10:14,8:9)])
brand = brand[order(brand)]



mycolor<-c("#E6AB02","#7570B3","#D95F02","#1B9E77","#CAB2D6",
           "#FDBF6F","#F7AEB3","#B2DF8A","#A6CEE3","#999999","#F780BF",
           "#D6372E","#FADD4B")
library(RColorBrewer)
library(randomcoloR)
brand_color = structure(mycolor, names = levels(brand))   
model_color = structure(rep(mycolor,as.numeric(table(brand))), names = names(brand))

library(RColorBrewer)
library(randomcoloR)
brand_color = structure(mycolor, names = levels(brand))
model_color = structure(rep(mycolor,as.numeric(table(brand))), names = names(brand))

##########
pdf("Hallmarker_immune_circos.pdf", width = 5, height = 5)
library(circlize)
circos.clear()

# 调小不同类之间的间距
gap.after = do.call("c", lapply(table(brand), function(i) c(rep(0.3, i - 1), 2)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))

# 绘制 chord diagram
chordDiagram(
  circos_inter,
  order = names(brand),
  grid.col = model_color,
  directional = 1,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = 0.02)
  ),
  annotationTrackHeight = mm_h(c(2, 2))
)
# 添加垂直排列的文字标签（圈外）
circos.track(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  
  circos.text(
    mean(xlim),
    max(ylim) + mm_h(20),  # 圈外
    labels = sector.index,
    col = "black",
    cex = 0.5,
    facing = "reverse.clockwise",
    niceFacing = TRUE
  )
}, bg.border = NA)

# 高亮 sector 并添加分类标签
for (b in unique(brand)) {
  model = names(brand[brand == b])
  highlight.sector(
    sector.index = model,
    track.index = 1,
    col = brand_color[b],
    text = b,
    text.vjust = -1,
    niceFacing = TRUE
  )
}

circos.clear()
dev.off()























circos_inter$from <- substr(circos_inter$from, nchar(circos_inter$from) - 2, nchar(circos_inter$from))
pdf("Hallmarker_immune_circosls.pdf",width = 4,height = 4)
library(circlize)
circos.clear()
gap.after = do.call("c", lapply(table(brand), function(i) c(rep(0.3, i-1), 4)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))

chordDiagram(circos_inter, order = names(brand),
             grid.col = model_color,
             directional = 1, annotationTrack = "grid", preAllocateTracks = list(
               list(track.height = 0.02)
             ),
             annotationTrackHeight = mm_h(c(2, 2))
)
#
circos.track(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), sector.index, col = "black", cex =0.5, 
              facing = "reverse.clockwise", niceFacing = T)
}, bg.border = NA)
# 
for(b in unique(brand)) {
  model = names(brand[brand == b])
  highlight.sector(sector.index = model, track.index = 1, col = brand_color[b],
                   text = b, text.vjust = -1, niceFacing = TRUE
  )
}

circos.clear()
dev.off()















pdf("Hallmarker_immune_circos8.pdf", width = 10, height = 10)

library(circlize)

# 清除之前的图
circos.clear()

# 设置每个 sector 之间的间距
gap.after = do.call("c", lapply(table(brand), function(i) c(rep(0.3, i - 1), 4)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))

# 绘制 chord diagram
chordDiagram(
  circos_inter,
  order = names(brand),
  grid.col = model_color,
  directional = 1,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = 0.02)  # 用于放置标签
  ),
  annotationTrackHeight = mm_h(c(2, 2))
)

# 添加 sector 标签（放在圈外侧）
circos.track(
  track.index = 1,  # 使用 track 1（圈外侧）
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")

    circos.text(
      mean(xlim),
      max(ylim) + mm_h(15),  # 向上偏移一点，放在圈外
      labels = sector.index,
      col = "black",
      cex = 0.6,
      facing = "clockwise",
      niceFacing = TRUE
    )
  },
  bg.border = NA
)

# 高亮 sector 并添加分类标签
for (b in unique(brand)) {
  model = names(brand[brand == b])
  highlight.sector(
    sector.index = model,
    track.index = 1,
    col = brand_color[b],
    text = b,
    text.vjust = -2,  # 向上偏移，避免重叠
    niceFacing = TRUE
  )
}

# 清除图形
circos.clear()
dev.off()

gap.after = do.call("c", lapply(table(brand), function(i) c(rep(1, i - 1), 5)))




##############
pdf("Hallmarker_immune_circos9.pdf", width = 6, height = 6)
library(circlize)
circos.clear()

# 增大间隔
gap.after = do.call("c", lapply(table(brand), function(i) c(rep(1, i - 1), 6)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))

# 绘制 chord diagram
chordDiagram(
  circos_inter,
  order = names(brand),
  grid.col = model_color,
  directional = 1,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = 0.02)
  ),
  annotationTrackHeight = mm_h(c(2, 2))
)

# 添加 sector 标签（圈外）
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")

    circos.text(
      mean(xlim),
      max(ylim) + mm_h(10),  # 文字往外放一大截
      labels = sector.index,
      col = "black",
      cex = 0.6,
      facing = "clockwise",
      niceFacing = TRUE
    )
  },
  bg.border = NA
)

# 高亮 sector 并添加分类标签
for (b in unique(brand)) {
  model = names(brand[brand == b])
  highlight.sector(
    sector.index = model,
    track.index = 1,
    col = brand_color[b],
    text = b,
    text.vjust = -3,  # 也适当上调分类标签
    niceFacing = TRUE
  )
}

circos.clear()
dev.off()