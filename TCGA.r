rm(list = ls())
setwd("D:/Project/TCGA/ESCA")


# Step1 数据整合 --------------------------------------------------------------
library("rjson")#install.packages("rjson") 
library(jsonlite)#install.packages("jsonlite") 
json <- jsonlite::fromJSON("metadata.cart.2024-08-30.json") 

sample_id <- sapply(json$associated_entities, function(x){x[, 1]})
file_sample <- data.frame(sample_id, file_name=json$file_name)  

count_file <- list.files('gdc_download_20240830_085953.118358', pattern = '*gene_counts.tsv', recursive = TRUE)
count_file_name <- strsplit(count_file, split='/')
count_file_name <- sapply(count_file_name, function(x){x[2]})

matrix = data.frame(matrix(nrow=60660, ncol=0))

for (i in 1:length(count_file_name)){
  path = paste0('gdc_download_20240830_085953.118358//', count_file[i])
  data<- read.delim(path, fill = TRUE, header = FALSE, row.names = 1)
  colnames(data)<-data[2, ]
  data <-data[-c(1:6), ]
  data <- data[3] 
  #取出unstranded列（得到COUNT矩阵），
  #若想提取tpm_unstranded则改为data[6]，
  #若想提取fpkm-unstranded则改为data[7]，
  #fpkm-up-unstranded改为data[8]，
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix, data)
}

matrix1 = data.frame(matrix(nrow=60660, ncol=0))

for (i in 1:length(count_file_name)){
  path = paste0('gdc_download_20240830_085953.118358//', count_file[i])
  data<- read.delim(path, fill = TRUE, header = FALSE, row.names = 1)
  colnames(data)<-data[2, ]
  data <-data[-c(1:6), ]
  data <- data[6] 
  #取出unstranded列（得到COUNT矩阵），
  #若想提取tpm_unstranded则改为data[6]，
  #若想提取fpkm-unstranded则改为data[7]，
  #fpkm-up-unstranded改为data[8]，
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix1 <- cbind(matrix1, data)
}

path = paste0('gdc_download_20240830_085953.118358//', count_file[1])
data<- as.matrix(read.delim(path, fill = TRUE, header = FALSE, row.names = 1))
gene_name <-data[-c(1:6), 1]

matrix0 <- cbind(gene_name, matrix)
matrix0 <- aggregate( . ~ gene_name, data=matrix0, max)   
rownames(matrix0) <- matrix0[, 1]
matrix0 <- matrix0[, -1]


matrix2 <- cbind(gene_name, matrix1)
matrix2 <- aggregate( . ~ gene_name, data=matrix2, max)    
rownames(matrix2) <- matrix2[, 1]
matrix2 <- matrix2[, -1]

sample <- colnames(matrix0)

normal <- c()
tumor <- c()

for (i in 1:length(sample)){
  if((substring(colnames(matrix0)[i], 14, 15)>=10)){    #14、15位置大于10的为normal样本
    normal <- append(normal, sample[i])
  } else {
    tumor <- append(tumor, sample[i])
  }
}

tumor_matrix <- matrix0[, tumor]
normal_matrix <- matrix0[, normal]

save(matrix0,matrix2,  file = 'step1_output.Rdata')



# Step2 数据预处理 ------------------------------------------------------------------
rm(list = ls()) #清空当前环境
load('step1_output.Rdata') # 加载前一步的过程文件

# 创建一个数据框 samples，其中列名为数据矩阵 matrix0 的列名（样本名）
samples <- as.data.frame(colnames(matrix0))
colnames(samples) <- 'samples' # 将列名设为 'samples'
samples$group <- NA # 创建一个空列 'group'

# 根据样本名中的特定位置信息判断每个样本所属的组别 ('normal' 或 'tumor')
samples$group <- ifelse(as.integer(substr(samples$samples, 14, 15)) < 10, 'tumor', 'normal')

# 将 'group' 列转换为因子，并指定水平顺序为 'normal' 和 'tumor'
samples$group <- factor(samples$group, levels = c('normal', 'tumor'))

# 按照 'group' 列对 samples 数据框进行排序
samples <- samples[order(samples$group), ]

# 提取排序后的 'group' 列，作为后续操作的组别列表
group_list <- samples$group

# 匹配 samples$samples 中的样本名，获取其在 matrix0 和 matrix2 中的索引位置
sample_indices <- match(samples$samples, colnames(matrix0))

# 根据样本索引提取 matrix0 中的表达计数(Count)数据
exp_counts <- matrix0[, sample_indices]

# 根据样本索引提取 matrix2 中的表达 TPM 数据
exp_TPM <- matrix2[, sample_indices]

# 根据 TPM 数据中每行大于1的计数比例筛选基因
filter_count <- apply(exp_TPM, MARGIN = 1, FUN = function(x){
  sum(x > 1) >= ncol(exp_TPM) * 0.5
}) # 计算每行中大于1的个数是否超过总列数的一半

# 根据筛选条件过滤 exp_counts 和 exp_TPM 中的基因
exp_counts <- exp_counts[filter_count, ]
exp_TPM <- exp_TPM[filter_count, ]


save(exp_counts, exp_TPM, group_list, samples, file = 'step2_output.Rdata')



# Step3 差异分析 ------------------------------------------------------------------
rm(list = ls())
load('step2_output.Rdata')

# 将 exp_counts 和 exp_TPM 中的所有元素转换为数值类型
exp_counts[] <- lapply(exp_counts, as.numeric)
exp_TPM[] <- lapply(exp_TPM, as.numeric)

#DESeq2包
library(DESeq2)
library(tidyverse)
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("DESeq2")

# 设置条件变量
condition = factor(group_list) #levels = c("normal", "tumor")
coldata <- data.frame(row.names = colnames(exp_counts), condition)

# 创建 DESeq2 数据集对象
dds <- DESeqDataSetFromMatrix(countData = exp_counts, colData = coldata, design = ~  condition)

# 进行差异分析
dds <- DESeq(dds)   
result <- as.data.frame(results(dds))
result$genes <- rownames(result)

# 根据 log2FoldChange 和 pvalue 计算基因的变化情况
log2FC_t= 2
change = ifelse(result$padj>0.05,'Stable',
                ifelse(abs(result$log2FoldChange) < log2FC_t,'Stable',
                       ifelse(result$log2FoldChange >= log2FC_t,'Up','Down') ))
table(change)
result <- mutate(result, change)
result <- result[order(result$log2FoldChange, decreasing = TRUE), ]
# 统计各类变化的基因数目
table(result$change)
# 保存差异分析结果
write.csv(result,'result.csv')

# 绘制火山图
#Volcano Plot
library(ggplot2)
fivenum(result$log2FoldChange)
fivenum(-log10(result$padj))
volcano_plot <- ggplot(result,aes(log2FoldChange,
                                  -log10(pvalue)))+
  geom_point(size = 1.5, 
             alpha = 0.8, 
             aes(color = change),
             show.legend = T)+
  scale_color_manual(values = c('#00468BFF','gray','#ED0000FF'))+
  ylim(0, 50)+
  xlim(-10, 10)+
  labs(x = 'Log2FoldChange',y = '-Log10(P.Value Adjusted)')+
  geom_hline(yintercept = -log10(0.05),
             linetype = 2,
             color = 'black',lwd = 0.8)+
  geom_vline(xintercept = c(-2, 2),
             linetype = 2, 
             color = 'black', lwd = 0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(volcano_plot)

ggsave(filename = "volcano_plot.tiff", plot = volcano_plot, 
       device = "tiff", dpi = 300,
       width = 9, height = 8)

# 绘制热图
#Heatmap Plot
library(pheatmap)#install.packages("pheatmap")
library(ggplot2)
deg_genes <- result %>%
  filter(change %in% c('Up', 'Down'))

exp_deg <- exp_TPM[rownames(exp_TPM) %in% deg_genes$genes,]
table(deg_genes$change)

# 将表达数据取对数并进行标准化处理
library(ComplexHeatmap)
exp_deg <- log2(exp_deg + 1)
exp_hot <- apply(exp_deg, 1, scale)
rownames(exp_hot) <- colnames(exp_deg)
exp_hot <- t(exp_hot)

# 设置颜色渐变函数
library(circlize)
col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c("#00468BFF", "white", "#ED0000FF")
)

# 绘制热图
heatmap_plot <- Heatmap(exp_hot, name = "Expression", cluster_rows = T, row_km = 2, 
                        col = col_fun ,show_column_names = F , show_row_names = F, 
                        row_gap = unit(2, "mm"),  # 关闭行之间的间隙
                        column_gap = unit(2, "mm"),  # 关闭列之间的间隙
                        clustering_distance_columns = 'pearson', row_title  = "Genes", 
                        row_title_side = "left" , column_km = 2, column_title = NULL ,
                        top_annotation = HeatmapAnnotation(
                          foo = anno_block(# 设置填充色
                            gp = gpar(fill = c("#00468BFF","#cD0000FF")),
                            # 设置文本标签
                            labels = c("Normal", "Tumor"), 
                            # 文本标签样式
                            labels_gp = gpar(col = "white", fontsize = 10)
                          )))
print(heatmap_plot)

tiff(filename = "heatmap_plot.tiff", width = 9, height = 8, units = "in", res = 300)  # Set resolution as needed
draw(heatmap_plot)
dev.off()

save(deg_genes, samples, exp_counts, exp_TPM, file = 'step3_output.Rdata')
save('volcano_plot','heatmap_plot', file = 'DEGs_plot.Rdata')



# Step4 富集分析 ------------------------------------------------------------------
rm(list = ls()) #清空
load('step3_output.Rdata')
library(clusterProfiler)
library(org.Hs.eg.db)#将Symbol转换为基因符号Entrez ID
library(dplyr)
library(DOSE)
library(ggplot2)
library(tidyr)
gene <- deg_genes$genes

#将Symbol转换为EntrezID
diff_entrez<-bitr(
  gene,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Hs.eg.db'
)
head(diff_entrez)

#GO富集
go_enrich<-clusterProfiler::enrichGO(gene = diff_entrez$ENTREZID,
                                     ont = 'all',#可选'BP','CC','MF' or 'all'
                                     keyType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",#p值矫正方法
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)

go_enrich<-DOSE::setReadable(go_enrich,
                             OrgDb = org.Hs.eg.db,
                             keyType = 'ENTREZID')

#去除冗余的GO term
go_geo<- simplify(go_enrich, cutoff=0.7, by="p.adjust",
                  select_fun=min)
#去除冗余的GO term
go_result<-go_geo@result
head(go_result)

library(ggplot2)
#dot
go_dot <- dotplot(go_geo,
                  x = "GeneRatio",
                  color = "p.adjust",
                  showCategory=10,
                  split='ONTOLOGY',
                  label_format = Inf)+#不换行
  #分面
  facet_grid(ONTOLOGY~.,
             space = 'free_y',#面板大小根据y轴自行调整
             scale='free_y'#子图坐标轴根据y轴自行调整
  )
go_dot
ggsave("go_dot.tiff", plot = go_dot, device = "tiff", 
       width = 9,  height = 9, 
       units = "in",dpi = 300)
#barplot
go_bar <- barplot(go_geo,
                  x = "Count",
                  color = "p.adjust",
                  showCategory=10,
                  split='ONTOLOGY',
                  label_format = Inf)+
  facet_grid(ONTOLOGY~.,
             space = 'free_y',#面板大小根据y轴自行调整
             scale='free_y'#子图坐标轴根据y轴自行调整
  )
go_bar
ggsave("go_bar.tiff", plot = go_bar, device = "tiff", 
       width = 9,  height = 9, 
       units = "in",dpi = 300)


#KEGG富集
KEGG_enrich <- clusterProfiler::enrichKEGG(gene = diff_entrez$ENTREZID,
                                           organism = "hsa", #物种Homo sapiens 
                                           pvalueCutoff = 0.05,#pvalue阈值
                                           qvalueCutoff = 0.05,#qvalue阈值
                                           pAdjustMethod = "BH",#p值矫正方法，one of "holm", 
                                           #"hochberg", "hommel", 
                                           #"bonferroni", "BH", 
                                           #"BY", "fdr", "none"
                                           minGSSize = 10,#富集分析中考虑的最小基因集合大小
                                           maxGSSize = 500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_enrich<-setReadable(KEGG_enrich,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENTREZID')
#提取KEGG富集结果表格                 
KEGG_result<-KEGG_enrich@result
head(KEGG_result)
# write.csv(KEGG_result, file = "KEGG_result.csv")
#dotplot
kegg_dot <- dotplot(KEGG_enrich,
                    x = "GeneRatio",
                    color = "p.adjust",
                    showCategory = 10)#展示top10通路
kegg_dot
ggsave("kegg_dot.tiff", plot = kegg_dot, device = "tiff", 
       width = 9,  height = 9, 
       units = "in",dpi = 300)

#barplot
kegg_bar <- barplot(KEGG_enrich,
                    x = "Count",
                    color = "p.adjust",
                    showCategory = 10)
kegg_bar

ggsave("kegg_bar.tiff", plot = kegg_bar, device = "tiff", 
       width = 9,  height = 9, 
       units = "in",dpi = 300)


###GSEA
symbol <- rownames(deg_genes)
entrez <- bitr(symbol, 
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = 'org.Hs.eg.db'
)
head(symbol)


gene_list <- deg_genes$log2FoldChange
names(gene_list) <- rownames(deg_genes)

gene_list <- gene_list[names(gene_list) %in% entrez[,1]]
names(gene_list) <- entrez[match(names(gene_list), entrez[, 1]), 2]
length(gene_list)
head(gene_list)

gene_list <- sort(gene_list, decreasing = T)
head(gene_list)


#KEGG-GSEA
KEGG_gse <- gseKEGG(geneList = gene_list, 
                    organism = "hsa", 
                    minGSSize = 5, 
                    maxGSSize = 500, 
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH", 
                    verbose = FALSE, 
                    eps = 0)
KEGG_gse <- setReadable(KEGG_gse, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID")
KEGG_gse_result <- KEGG_gse@result


#gsea排序：按照NES排

gsea <- read.csv("D:\\Project\\TCGA\\ESCA\\KEGG_gse_result.csv",row.names = 1)
gsea <- gsea[order(gsea$NES,decreasing = T),]


# write.csv(KEGG_gse_result, file = 'KEGG_gse_result.csv')
#可视化
library(ggridges)
library(ggplot2)
library(enrichplot)

###山峦图
kegg1 <- ridgeplot(KEGG_gse, 
                   showCategory = 10, 
                   fill = "p.adjust", 
                   decreasing = T)
kegg1

###气泡图
kegg2 <- dotplot(KEGG_gse,
                 x = "GeneRatio",
                 color = "p.adjust",
                 showCategory = 10)
kegg2
ggsave("KEGG_gsea.tiff", plot = kegg2, device = "tiff", 
       width = 9,  height = 9, 
       units = "in",dpi = 300)

x1 <- KEGG_result$Description
x2 <- KEGG_gse_result$Description
m <- intersect(x1,x2)
m
#ES图绘制(gseaplot2函数)
#单个基因集可视化：
gsea_kegg <- gseaplot2(KEGG_gse,
                       geneSetID = 4,
                       color = "red",
                       rel_heights = c(1.5, 0.5, 1), #子图高度
                       subplots = 1:3, #显示哪些子图
                       pvalue_table = F, #是否显示pvalue表
                       title = KEGG_gse$Description[3],
                       ES_geom = "line") #"dot"将线转换为点
gsea_kegg

ggsave("gsea_kegg.tiff", plot = gsea_kegg, device = "tiff", 
       width = 9,  height = 6, 
       units = "in",dpi = 300)


# Step5 机器学习 ------------------------------------------------------------------
rm(list = ls()) #清空
load('step3_output.Rdata')
rm(deg_genes,exp_counts)
library(tidyverse)

tnf <- read.table("D:\\Project\\TCGA\\ESCA\\clinic\\TNF.txt")

exp <- log2(exp_TPM + 1)

exp_tnf <- exp[rownames(exp) %in% tnf$V1,]
rm(exp,exp_TPM,tnf)
exp <- as.data.frame(t(exp_tnf))
exp$samples <- rownames(exp)

exp_ml <- merge(samples,exp, by = "samples")
# rownames(exp_ml) = exp_ml$samples
# exp_ml1 = exp_ml[,!names(exp_ml) %in% c('samples')]

genes = colnames(exp_ml)[3:198]

library(mlr3)
library(caret)
set.seed(123)

# 创建分割索引
split_index <- createDataPartition(exp_ml$group, p = 0.6, list = FALSE)

# 创建训练集和测试集
train_data <- exp_ml[split_index, ]
test_data <- exp_ml[-split_index, ]



train_data$group = ifelse(train_data$group=="normal",0,1)
y = as.factor(train_data$group)
x = train_data[,!names(train_data) %in% c('samples','group')]
x = as.matrix(x)


###lasso
library(glmnet)
set.seed(111)
fit=glmnet(x,y,family = "binomial",maxit = 10000, nfold = 10)
plot(fit,xvar="lambda",label = TRUE)
cvfit = cv.glmnet(x,y,family="binomia",maxit = 10000, nfolds = 10)
plot(cvfit)

coef=coef(fit,s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)#查看模型的相关系数geneCoef
geneCoef

lassoGene <- lassoGene[-1]
lassoGene

##RandomForest
library(randomForest)
set.seed(100)
rf <- randomForest(y~.,  data = x , ntree = 500)
plot(rf, main = 'Random Forest', lwd = 2)

optionTrees = which.min(rf$err.rate[, 1])
#rf2 = randomForest(y~., data = x, ntree = optionTrees, importance = T)
rf2 = randomForest(y~., data = x, ntree = optionTrees)

# 计算变量重要性
importance = importance(x = rf2)
varImpPlot(rf2, main = 'Feature Importance')
# 筛选重要性大于5的变量
rfGenes = importance[order(importance[, 'MeanDecreaseGini'], decreasing = T), ]
mean_importance = mean(rfGenes)
rfGenes = names(rfGenes[rfGenes > mean_importance])
rfGenes


###Boruta
library(Boruta)
set.seed(1124)
x1 = train_data[, -1]
x1$group <- as.factor(x1$group)

Var.Selec<-Boruta(
  group~.,
  x1,
  pValue = 0.01, #confidence level. 可信水平
  mcAdj = TRUE, #是否使用Bonferroni调整
  #if set to TRUE, a multiple comparisons adjustment using the Bonferroni method will be applied.
  maxRuns = 500, #迭代最多次数
  doTrace =0,#可以选0-3，运行结果显示的详细程度，0不显示轨迹
  holdHistory = TRUE, #如果设置为TRUE，则存储完整的重要性历史记录，并将其作为结果的ImpHistory元素返回。
  getImp = getImpRfZ #用于获取属性重要性的函数。默认值是 getImpRfZ，它从 ranger 包运行随机森林并收集平均降低精度测量的 Z 分数。
)
Var.Selec

#（绿色是重要的变量，红色是不重要的变量，蓝色是影子变量，黄色是Tentative变量）。
tiff("Boruta1.tiff",width = 6*300, height = 6*300, res = 300)
plotImpHistory(Var.Selec,
               ylab="Z-Scores",
               las=2)
dev.off()

tiff("Boruta2.tiff",width = 6*300, height = 6*300, res = 300)
plot(Var.Selec,
     whichShadow=c(F,F,F),
     xlab="",
     ylab="Z-Scores",
     las=2)
dev.off()

getConfirmedFormula(Var.Selec)
getNonRejectedFormula(Var.Selec)
boruta_genes <- getSelectedAttributes(Var.Selec,withTentative=FALSE)
boruta_genes


lassoGene
rfGenes
boruta_genes

both <- Reduce(intersect, list(lassoGene, rfGenes, boruta_genes))
both


###训练集roc
roc_exp <- as.data.frame(x[,both])
roc_exp$group <- y
roc_exp$group <- as.factor(roc_exp$group)

library(pROC)
selected_genes <- both
library(ggsci)
cors = pal_lancet(alpha = 0.7)(7)



auc_values <- sapply(selected_genes, function(both) {
  roc_obj <- roc(roc_exp$group, roc_exp[[both]], levels = c("0", "1"))
  auc(roc_obj)
})

sorted_indices <- order(auc_values, decreasing = TRUE)
sorted_genes <- selected_genes[sorted_indices]

roc1 <- roc(roc_exp$group, roc_exp[[sorted_genes[1]]], levels = c("0", "1"))
plot(roc1, legacy.axes = TRUE, col = cors[1], 
     main = "ROC Curves for Selected Genes",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)",
     lwd = 2)

for (i in 2:length(sorted_genes)) {
  roc_temp <- roc(roc_exp$group, roc_exp[[sorted_genes[i]]], levels = c("0", "1"))
  plot(roc_temp, add = TRUE, col = cors[i], lwd = 2)
}

#abline(a = 0, b = 1, lty = 2, col = "darkgray")

legend_text <- sapply(seq_along(sorted_genes), function(i) {
  paste0(sorted_genes[i], " - AUC: ", sprintf("%.3f", auc_values[sorted_genes[i]]))
})

legend("bottomright", legend = legend_text, 
       col = cors[1:length(sorted_genes)], lwd = 2,cex = 0.8)


#验证集roc
test_data$group = ifelse(test_data$group=="normal",0,1)
y = as.factor(test_data$group)
x = test_data[,!names(test_data) %in% c('samples','group')]
x = as.matrix(x)


roc_exp <- as.data.frame(x[,both])
roc_exp$group <- y
roc_exp$group <- as.factor(roc_exp$group)

library(pROC)
selected_genes <- both
library(ggsci)
cors = pal_lancet(alpha = 0.7)(7)



auc_values <- sapply(selected_genes, function(both) {
  roc_obj <- roc(roc_exp$group, roc_exp[[both]], levels = c("0", "1"))
  auc(roc_obj)
})

sorted_indices <- order(auc_values, decreasing = TRUE)
sorted_genes <- selected_genes[sorted_indices]

roc1 <- roc(roc_exp$group, roc_exp[[sorted_genes[1]]], levels = c("0", "1"))
plot(roc1, legacy.axes = TRUE, col = cors[1], 
     main = "ROC Curves for Selected Genes",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)",
     lwd = 2)

for (i in 2:length(sorted_genes)) {
  roc_temp <- roc(roc_exp$group, roc_exp[[sorted_genes[i]]], levels = c("0", "1"))
  plot(roc_temp, add = TRUE, col = cors[i], lwd = 2)
}

#abline(a = 0, b = 1, lty = 2, col = "darkgray")

legend_text <- sapply(seq_along(sorted_genes), function(i) {
  paste0(sorted_genes[i], " - AUC: ", sprintf("%.3f", auc_values[sorted_genes[i]]))
})

legend("bottomright", legend = legend_text, 
       col = cors[1:length(sorted_genes)], lwd = 2,cex = 0.8)


##逻辑回归
both
library(dplyr)
train_data1 <- train_data %>%
  select(samples, group, BCL2A1, BTG2, CXCL10, DUSP2, HES1, PHLDA2, YRDC)
train_data1$group <- as.factor(train_data1$group)
str(train_data1)

test_data1 <- test_data %>%
  select(samples, group, BCL2A1, BTG2, CXCL10, DUSP2, HES1, PHLDA2, YRDC)
test_data1$group <- as.factor(test_data1$group)
str(test_data1)


##训练集
logistic_model <- glm(group ~ BCL2A1 + BTG2 + CXCL10 + DUSP2 + HES1 + PHLDA2 + YRDC,
                      data = train_data1,
                      family = binomial())
summary(logistic_model)

library(ggsci)
cors=pal_lancet(alpha = 0.5)(3)
cors

library(pROC)
predictions <- predict(logistic_model, type = "response")
roc_obj <- roc(train_data1$group, predictions)

plot(roc_obj, print.auc=TRUE, col = '#00468B7F')
auc_value <- auc(roc_obj)
print(auc_value)

##验证集
test_predictions <- predict(logistic_model, newdata = test_data1, type = "response")
roc_obj_test <- roc(test_data1$group, test_predictions)

plot(roc_obj_test, print.auc=TRUE, col = "#ED00007F")
auc_value_test <- auc(roc_obj_test)
print(auc_value_test)
