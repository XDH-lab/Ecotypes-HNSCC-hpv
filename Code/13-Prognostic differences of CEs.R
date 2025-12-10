rm(list=ls())

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/bulk合并4/预后/ECO预后")

eco<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/bulk合并4/结果/exp_laster/Ecotypes/ecotype_assignment.txt")
head(eco)
lnc<-eco[,c(1,2,7,8)]

head(lnc)
table(lnc$MaxEcotype)
lnc<-na.omit(lnc)
colnames(lnc)[c(3,4)]<-c("fustat","futime")

#lnc$fustat <- ifelse(lnc$fustat == TRUE, 1, 0)

is.numeric(lnc$futime)
is.numeric(lnc$fustat)

ls<-lnc
head(lnc)


library(survival)
survival_object <- Surv(lnc$futime, lnc$fustat)
fit <- survfit(survival_object ~ lnc$MaxEcotype)
library(ggplot2)

library(survminer)

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/bulk合并4/预后/ECO预后/颜色")
col_sankey<-c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 
              'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono_Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF",
              'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3")
ggsurvplot(fit, data = lnc, conf.int = FALSE,  # 关闭置信区间阴影
                     xlab = "Time", ylab = "Survival probability",
                     pval = TRUE,
					 palette = c("#D6372E", "#FADD4B","#70B460","#E690C1", "#985EA8","#A3A3A3"))
ggsave("surviva_eco.pdf")



######
rm(list=ls())

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果")

eco<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/结果/TCGA-HNSCC最后/Ecotypes/ecotype_assignment.txt")
head(eco)
su<-eco[,c(1,2)]
setwd("/groups/g900008/home/zhangnihui/新课题/实验/验证/TCGA-HNSCC")

# 获取所有细胞类型的文件夹名称
cell_types <- list.files(path = "/lkn_lab/znh/12月课题实验/实验/验证集/bulk/TCGA/exp", full.names = FALSE)
cell_types <- cell_types[-which(cell_types == "Ecotypes")]

# 读取临床数据
data1 <- read.delim("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/clin-mouth.txt")
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/ECO预后")

lnc <- cbind(data1[,1], data1$OS.time, data1$OS) # rownames(lnc)<-data1[,1]
lnc <- as.data.frame(lnc)
colnames(lnc) <- c("ID", "futime", "fustat")
lnc$futime <- as.numeric(lnc$futime)
lnc$fustat <- as.numeric(lnc$fustat)
ls <- lnc
head(lnc)
head(su)
lnc <- lnc %>%
  mutate(ID = str_replace_all(ID, "-", "\\."))
lnc<-merge(lnc, su, by = "ID", all.x = TRUE)
table(lnc$MaxEcotype)
lnc<-na.omit(lnc)

#lnc$fustat <- ifelse(lnc$fustat == TRUE, 1, 0)

is.numeric(lnc$futime)
is.numeric(lnc$fustat)

ls<-lnc
head(lnc)



library(survival)
survival_object <- Surv(lnc$futime, lnc$fustat)
fit <- survfit(survival_object ~ lnc$MaxEcotype)
library(ggplot2)

library(survminer)

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/ECO预后")
col_sankey<-c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 
              'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono_Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF",
              'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3")
ggsurvplot(fit, data = lnc, conf.int = FALSE,  # 关闭置信区间阴影
                     xlab = "Time", ylab = "Survival probability",
                     pval = TRUE,
					 palette = c("#D6372E", "#FADD4B","#70B460","#E690C1", "#985EA8","#A3A3A3"))
ggsave("surviva_eco.pdf")

