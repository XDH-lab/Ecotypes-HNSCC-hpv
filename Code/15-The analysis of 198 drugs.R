#198种药物敏感性得分预测
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
rm(list=ls())
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/药物敏感性")
GDSC2_Expr = readRDS('/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/药物反应/ls/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/药物反应/ls/GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 

 
testExpr=read.delim('/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/TCGA-HNSCC最后.txt')
testExpr[1:5,1:5]
rownames(testExpr)=testExpr$Gene
testExpr=testExpr[,-1]
max(testExpr)
class(testExpr)
testExpr=as.matrix(testExpr)

 
library(preprocessCore)
dir.create('GDSC2')
setwd('GDSC2')
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat   qn
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' ,
              rsq = F,
              cc =F,
              report_pc =F,
              pcr=F
)
