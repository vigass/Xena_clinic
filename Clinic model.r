rm(list = ls())
setwd("D:\\Project\\TCGA\\ESCA")

# Step1 -------------------------------------------------------------------
load('step2_output.Rdata')

library(tidyverse)

gene_TNF <- read.table("D:\\Project\\TCGA\\ESCA\\clinic\\TNF.txt")
genes <- gene_TNF$V1

# genes <- rownames(exp_TPM)

exp1 <- exp_TPM[rownames(exp_TPM) %in% genes,]

exp <- as.data.frame(t(exp1))
rm(exp_counts,exp_TPM)

load("D:\\Project\\TCGA\\ESCA\\clinic\\pd.rdata")

head(substr(rownames(exp), 1, 16))

rownames(exp) <- substr(rownames(exp), 1, 16)
samples$samples <- substr(samples$samples, 1,16)


samples1 <- samples[samples$samples %in% pd$sample, ]
exp2 <- exp[rownames(exp) %in% samples1$samples, ]
pd1 <- pd[pd$sample %in% samples1$samples, ]

save(pd1,exp2,samples1, file = "clinic1_output.rdata")



# Step2  ------------------------------------------------------------------
rm(list = ls())
load("clinic1_output.rdata")
exp2[] <- lapply(exp2, as.numeric)
exp2 <- log2(exp2+1)
exp2$sample <- rownames(exp2)
rownames(exp2) <- seq_along(rownames(exp2))
pd1$OS.time=pd1$OS.time/365
#机器学习
##单因素筛选
cox_matrix <- pd1[,c('sample','OS','OS.time')]

cox <- merge(cox_matrix,exp2, by="sample")

library(survival)
pFilter=0.05 #设一个p值标准，后面用
outResult=data.frame() #建一个空白数据框，后面for循环输出用
sigGenes=c("OS","OS.time") #建一个向量，后面for循环输出用，因为后面还要用到OS及OS.time，所以先放在向量里
for(i in colnames(cox[,4:ncol(cox)])){ #从第3列开始循环，因为1列2列不是gene，是OS和OS.time
  coxcox <- coxph(Surv(OS.time, OS) ~ cox[,i], data = cox)#开始逐一循环cox分析
  coxcoxSummary = summary(coxcox) #summary命令对coxcox总结，方面后面提取数据
  pvalue=coxcoxSummary$coefficients[,"Pr(>|z|)"] #提取p值，这个主要是后面提取有意义的gene用
  if(pvalue<pFilter){ # 这里我们不需要提取所有基因的数据，只需要有意义的gene的结果，所以设置如果pvalue<0.05才提取数据
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,#合并行，实际上是对循环结果的合并，前面设置的空白数据框outResult这里用，循环必须有个开始
                    cbind(id=i,#合并列，是每个基因的统计数据
                          HR=coxcoxSummary$conf.int[,"exp(coef)"],#提取单个基因的HR
                          L95CI=coxcoxSummary$conf.int[,"lower .95"],#提取单个基因的HR的95%CI低值
                          H95CI=coxcoxSummary$conf.int[,"upper .95"],#提取单个基因的HR的95%CI高值
                          pvalue=coxcoxSummary$coefficients[,"Pr(>|z|)"])#提取单个基因的p值
    )
  }
}
head(outResult)
outResult =as.data.frame(outResult)
outResult[, -1] <- lapply(outResult[, -1], as.numeric)

#可视化
library(ggplot2)
genes=outResult$id
hr <- sprintf("%.3f",outResult$"HR")
hrLow  <- sprintf("%.3f",outResult$"L95CI")
hrHigh <- sprintf("%.3f",outResult$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(outResult$pvalue<0.001, "<0.001", sprintf("%.3f", outResult$pvalue))


pdf("UniCoxSurForestPlot.pdf",width = 6, height = 6)
n <- nrow(outResult)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

layout(matrix(c(1,2),nc=2),width=c(2.5,2))

xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,genes,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex)
text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex)
text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),
       n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)

dev.off()

outResult$id


##Lasso筛选
cox_genes=outResult$id
unique_cox_genes <- unique(cox_genes)

lasso_matrix <- pd1[,c('sample','OS','OS.time')]
lasso <- merge(lasso_matrix,exp2, by="sample")

lasso_genes <- colnames(lasso)[!(colnames(lasso) %in% c("sample", "OS", "OS.time"))]
common_genes <- intersect(lasso_genes, unique_cox_genes)
lasso <- lasso[, c("sample", "OS", "OS.time", common_genes)]

library(glmnet) #install.packages("glmnet")
library(survival) #install.packages("survival")
x <- as.matrix(lasso[, !(colnames(lasso) %in% c('sample', 'OS', 'OS.time'))])
y <- Surv(lasso$OS.time, lasso$OS)


lasso_model <- glmnet(x, y, family = "cox", 
                         alpha = 1) # alpha=1表示Lasso回归

par(mfrow = c(1,1))
plot(lasso_model)
lasso_model

set.seed(100)
fitCV <- cv.glmnet(x, y,                    
                   family = "cox",                   
                   type.measure = "deviance",                   
                   nfolds = 10)
plot(fitCV)

#提取lasso genes
coef=coef(lasso_model,s = fitCV$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)#查看模型的相关系数geneCoef
geneCoef

lassoGene <- lassoGene[-1]
lassoGene

save(lassoGene,lasso, 
     pd1,samples1,exp2, file = "clinic2_output.rdata")


# Step3  ------------------------------------------------------------------
rm(list = ls())
load("clinic2_output.rdata")

coxs <- lasso[, c("sample", "OS", "OS.time", lassoGene)]

library(survival)
library(survminer)

surv_obj <- Surv(coxs$OS.time, coxs$OS)

# 拟合多因素Cox回归模型
cox_model <- coxph(surv_obj ~ ., data = coxs[, c(lassoGene)])
# 查看模型结果
cox_summary <- summary(cox_model)
cox_summary
# 获取系数值
coef_values <- cox_summary$coefficients
# 打印系数值
print(coef_values)
#ggforest(cox_model, data = coxs)

coef_values=as.data.frame(coef_values)

gene_names <- rownames(coef_values)
coef_values <- coef_values[, "coef"]
coef_df <- data.frame(gene = gene_names, coef = coef_values)
# 计算风险评分
# 对 coxs 数据框中基因表达量与系数值进行矩阵乘法
risk_scores <- as.matrix(coxs[, lassoGene]) %*% coef_values
# 将风险评分添加到 coxs 数据框中
coxs$score <- risk_scores
# 打印风险评分数据框的前几行
fivenum(coxs$score)

# 计算均值和中位数
mean_score <- mean(coxs$score)
median_score <- median(coxs$score)
cat("Mean Score:", mean_score, "\n")
cat("Median Score:", median_score, "\n")

# 绘制风险评分分布的直方图
hist(coxs$score, breaks=30, main="Risk Score Distribution", xlab="Risk Score", col="lightblue")

#使用中位数作为评分标准
coxs$Risk_level <- NA
coxs$Risk_level <- ifelse(coxs$score>mean_score,'High','Low')
table(coxs$Risk_level)


merge_pd <- pd1[, !(names(pd1) %in% 'X_PATIENT')]
merge_coxs <- coxs[,c('sample',"score", 'Risk_level')]

final_coxs <- merge(merge_pd,merge_coxs,by = 'sample')

final_coxs$age <- ifelse(final_coxs$age <= 65, "<=65",">65")
#write.csv(final_coxs,"final_coxs.csv")

cols_to_convert <- c("grade", "T_stage", "M_stage", "N_stage", "gender", "Risk_level")
final_coxs[cols_to_convert] <- lapply(final_coxs[cols_to_convert], factor)

formula <- as.formula(paste("Surv(OS.time, OS) ~",
                            paste(c("T_stage", "M_stage", "N_stage", "age", "gender", "Risk_level"),
                                  collapse = " + ")))
cox_model <- coxph(formula, data = final_coxs)
summary(cox_model)
ggforest(cox_model, data = final_coxs)


#KA
library(ggsci)
cors <- pal_lancet(alpha = 0.6)(3)
library(survminer)
km_fit <- survfit(Surv(OS.time, OS) ~ Risk_level, data = final_coxs)
km_fit
km <- ggsurvplot(km_fit, data = final_coxs, 
                 pval = TRUE,                      #是否显示p值
                 risk.table = TRUE,                #是否风险表添加
                 conf.int = TRUE,                  #是否显示置信区间
                 title = "Kaplan-Meier Survival Curve for Risk score", 
                 xlab = "Time (Years)", ylab = "Survival Probability",
                 legend.title = "Risk score", legend.labs = c("Low", "High"),
                 palette = cors)
km

pdf("KM.pdf", height = 9, width = 9)  
print(km, newpage = FALSE)
dev.off()


#time ROC
library(timeROC)
roc_time <- final_coxs


ROC <- timeROC(T = roc_time$OS.time, 
               delta = roc_time$OS, 
               marker = roc_time$score,   #预测的生死（默认较大的预测值和高事件风险相关，如果较大预测值和低风险事件相关，需要添加一个“-”反转关联）
               cause = 1,                 #阳性事件结局（这里阳性事件为死亡，对应1）
               weighting = "marginal",    #权重计算方法，默认
               times = c(1,3,5),          #时间点划分：1年、3年、5年
               ROC = T,                   #是否保存TP和FP值
               iid = T)                   #是否计算ROC曲线下面积
ROC
pdf("roc_plot.pdf", width = 6, height = 6)
plot(ROC,title = F,
     time=1, col="#00468BFF")
plot(ROC,
     time=3, col="#ED0000FF",add = T)
plot(ROC,
     time=5, col="#42B540FF",add = T)
title(main="Time Dependent ROC")

legend("bottomright",paste0("AUC of ",c(1,3,5)," Years survival: ",
                            format(round(ROC$AUC,3),nsmall=2)), 
       col=c("#00468BFF", "#ED0000FF", "#42B540FF"),lty=1,lwd=2,bty = "n")
dev.off()



###风险评分散点图
dot_cox <- final_coxs
dot_cox <- dot_cox[order(dot_cox$score),]
dot_cox$id <- 1:122
library(ggplot2)
p1 <- ggplot(dot_cox, aes(x = id, y = score)) +
  geom_point(aes(col = Risk_level)) + 
  scale_color_manual(values = cors) + 
  geom_hline(yintercept = median(dot_cox$score), colour = "grey", linetype="dashed", linewidth =0.8) +
  geom_vline(xintercept = sum(dot_cox$Risk_level == "Low"), colour="grey", linetype = "dashed", linewidth  = 0.8) +
  theme_bw() + 
  xlab("Samples")
p1

##生存时间散点图
suv_cox <- final_coxs
suv_cox <- suv_cox[order(suv_cox$Risk_level),]
suv_cox$OS <- factor(suv_cox$OS, levels = c(0, 1), labels = c("Alive", "Dead"))

suv_cox$id <- 1:122
p2 <- ggplot(suv_cox, aes(x = id, y = OS.time)) +
  geom_point(aes(col = OS)) +
  scale_colour_manual(values = cors) +
  geom_vline(xintercept = sum(suv_cox$Risk == "Low"), colour = "grey", linetype = "dashed", linewidth  = 0.8) +
  labs(col = "Status")+ 
  xlab("Samples") + 
  theme_bw()

p2

##热图
heat_cox <- final_coxs
heat_cox <- heat_cox[order(heat_cox$Risk_level),]


heatexp_cox <- lasso[, !(colnames(lasso) %in% c("OS", "OS.time"))]
rownames(heatexp_cox) <- heat_cox$sample
heatexp_cox <- heatexp_cox[,-1]


match_positions <- match(heat_cox$sample,rownames(heatexp_cox))
heatexp_cox <- heatexp_cox[match_positions,]
heatexp_cox <- t(heatexp_cox)

hot <- apply(heatexp_cox, 1, scale)
rownames(hot) <- colnames(heatexp_cox)

table(heat_cox$Risk)

library(circlize)
col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c("#00468BFF", "white", "#ED0000FF")
)
library(ComplexHeatmap)
heatmap_plot <- Heatmap(t(hot), name = "Expression", cluster_rows = F, row_km = 2,cluster_columns = F,
                        col = col_fun ,show_column_names = F , show_row_names = T,
                        row_gap = unit(0, "mm"),  # 关闭行之间的间隙
                        column_gap = unit(0, "mm"),  # 关闭列之间的间隙
                        clustering_distance_columns = 'pearson', row_title  = "Genes",
                        row_title_side = "left" , column_km = 2, column_title = NULL ,
                        top_annotation = HeatmapAnnotation(
                          foo = anno_block(# 设置填充色
                            gp = gpar(fill = c("#00468BFF","#cD0000FF")),
                            # 设置文本标签
                            labels = c("Low Risk", "High Risk"),
                            # 文本标签样式
                            labels_gp = gpar(col = "white", fontsize = 10)
                          )))
heatmap_plot

# 
# table(heat_cox$Risk)
# sample_group_info <- factor(heat_cox$Risk, levels = c("Low", "High"))
# 
# heatmap_plot <- Heatmap(t(hot), 
#                         name = "Expression", 
#                         cluster_rows = TRUE, 
#                         cluster_columns = F,
#                         row_km = 2, 
#                         col = col_fun,
#                         show_column_names = FALSE, 
#                         show_row_names = TRUE, 
#                         row_gap = unit(0, "mm"),  # No gap between rows
#                         column_gap = unit(0, "mm"),  # No gap between columns
#                         clustering_distance_columns = 'pearson', 
#                         row_title = "Genes", 
#                         row_title_side = "left",
#                         column_title = NULL,
#                         column_split = sample_group_info,  # Use your sample group info here
#                         top_annotation = HeatmapAnnotation(
#                           foo = anno_block(
#                             gp = gpar(fill = c("#00468BFF", "#cD0000FF")),
#                             labels = c("Low Risk", "High Risk"),
#                             labels_gp = gpar(col = "white", fontsize = 10)
#                           )))
# heatmap_plot


p1
ggsave("dot1_risk.pdf",p1, height = 6,width = 6,dpi = 300)
p2
ggsave("dot2_risk.pdf",p2, height = 6,width = 6,dpi = 300)
draw(heatmap_plot)
pdf("heatmap_risk.pdf", height = 6, width = 6)
draw(heatmap_plot)
dev.off()


# Step4  ------------------------------------------------------------------
rm(list = ls())
library(survival)
library(rms)
final_coxs <- read.csv("final_coxs.csv",row.names = 1)
cols_to_convert <- c("grade", "T_stage", "M_stage", "N_stage","age",
                     "gender", "Risk_level")
final_coxs[cols_to_convert] <- lapply(final_coxs[cols_to_convert], factor)

dd <- datadist(final_coxs)
options(datadist = dd)


# formula <- as.formula(paste("Surv(OS.time, OS) ~",
#                             paste(c("grade", "T_stage", "M_stage", "N_stage","age",
#                                     "gender", "Risk_level"),
#                                   collapse = " + ")))
# cox_model <- cph(formula, data = final_coxs, surv = TRUE)

cox_model1 <- coxph(Surv(OS.time,OS) ~ grade+T_stage+M_stage+N_stage+age+gender+Risk_level, 
                   data = final_coxs)

cox_model1

library(rms)
dd <- datadist(final_coxs)
options(datadist = dd)
surv_object = with(final_coxs,Surv(OS.time,OS == 1))
cox_model2 <- cph(surv_object ~ grade+T_stage+M_stage+N_stage+age+gender+Risk_level, 
                  data = final_coxs,x=TRUE,y=TRUE,surv=TRUE)
print(cox_model2)

#创建一个与生存分析相关的Survival对象
surv<-Survival(cox_model2)
#定义用于不同时间点的生存函数
surv1<- function(x) surv(1,x) #一年生存概率预测
surv3<- function(x) surv(3,x) #三年生存概率预测
surv5<- function(x) surv(5,x) #五年生存概率预测

#绘制nomogram，可视化cox比例风险模型的结果
nom <- nomogram(cox_model2,
                fun = list(surv1, surv3, surv5),
                lp = T, 
                funlabel = c("1-year survival", "3-year survival", "5-year survival"),
                fun.at = list(seq(0.05,0.95,by = 0.05), 
                              seq(0.05,0.95,by = 0.05), 
                              seq(0.05,0.95,by = 0.05))
                )# Set survival probability axis
plot(nom)



# library(regplot)
# regplot(cox_model1, 
#         observation = final_coxs[121, ],
#         points = T, 
#         plots = c("violin", "boxes"),
#         failtime = c(1,3,5), 
#         odds = F, 
#         leftlabel=T, 
#         prfail=TRUE,
#         showP=T,
#         droplines=T,
#         #colors="red",
#         rank="range",
#         interval="confidence",
#         title="Cox regression"
#         )

####Hosmer-Lemeshow拟合优度检验
#使用calibrate函数创建一个校准对象cal1
# cal1<-calibrate(cox_model2,
#                   cmethod='KM',#表示使用Kaplan-Meier（KM）估计方法进行校准
#                   method="boot",#表示使用自助法（Bootstrap）进行校准，Bootstrap是一种统计方法，它通过从原始数据中有放回地进行重采样来估计参数的不确定性和分布。在这里，Bootstrap 用于生成多个随机样本来估计校准曲线的分布，以便获得更可靠的校准结果。
#                   u=1,#设置时间间隔，需要与之前模型中定义的time.inc一致
#                   m=10,#每次抽样的样本量，根据样本量来确定，标准曲线一般将所有样本分为3组（在图中显示3个点）
#                   B=1000)#抽样次数
# 
# cal3 <- calibrate(cox_model2,
#                   cmethod = 'KM',
#                   method = "boot",
#                   u = 3,
#                   m = 10,
#                   B = 1000)
# 
# cal5 <- calibrate(cox_model2,
#                   cmethod = 'KM',
#                   method = "boot",
#                   u = 5,
#                   m = 10,
#                   B = 1000)
# 
# 
# #绘制校准曲线
# par(mar=c(6,6,3,3))
# plot(cal1,#绘制校准曲线的数据
#      lwd=1,#线条宽度为1
#      lty=1,#线条类型为1（实线）
#      conf.int=T,#是否显示置信区间
#      errbar.col="blue3",#直线和曲线的误差线颜色设置为蓝色
#      col="red3",#校准曲线的颜色设置为红色
#      xlim=c(0,1),#x轴的限制范围，从0到1
#      ylim=c(0,1),#y轴的限制范围，从0到1
#      xlab="Nomogram-Predicted Probability of 1-YearDFS",#x轴标签
#      ylab="Actual 1-Year DFS (proportion)", #y轴标签
#      subtitles = F)#不显示副标题
# 
# # 添加3年的校准曲线
# plot(cal3,
#      add = TRUE,
#      lwd = 1,
#      lty = 2,  # 使用不同的线型
#      col = "blue3")  # 使用不同的颜色
# 
# # 添加5年的校准曲线
# plot(cal5,
#      add = TRUE,
#      lwd = 1,
#      lty = 3,  # 使用不同的线型
#      col = "green3")  # 使用不同的颜色
# 
# # 添加图例
# legend("bottomright",
#        legend = c("1-Year DFS", "3-Year DFS", "5-Year DFS"),
#        col = c("red3", "blue3", "green3"),
#        lty = c(1, 2, 3),
#        lwd = 1)
