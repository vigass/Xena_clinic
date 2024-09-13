rm(list = ls())
setwd("D:/Project/TCGA/ESCA/clinic")

library(data.table)
library(tidyverse)
clinical = fread("TCGA-ESCA.GDC_phenotype.tsv.gz",header = T, sep = '\t',data.table = F)
surv = read.delim("TCGA-ESCA.survival.tsv",header = T)
# # rm(list = ls())


tmp = data.frame(colnames(clinical))

# stage
table(clinical$tumor_stage.diagnoses)
# M分期
table(clinical$pathologic_M)
# T分期
table(clinical$pathologic_T)
# N分期
table(clinical$pathologic_N)
# grade
table(clinical$tumor_grade.diagnoses)
table(clinical$neoplasm_histologic_grade)

meta <- clinical
tmp <- as.data.frame(colnames(meta))

meta = meta[,c(
  'submitter_id.samples',
  'neoplasm_histologic_grade',
  'pathologic_T' ,
  'pathologic_M',
  'pathologic_N',
  'age_at_index.demographic',
  'gender.demographic' ,
  'tumor_stage.diagnoses')]
colnames(meta)=c('ID',
                 'grade',
                 'T_stage','M_stage','N_stage','age','gender','stage')

meta_org = meta
#检查各列的信息是否规范 没有冗余信息 缺失的信息用NA代替
#——grade——————————————————————————————————
table(meta$grade)
meta <- meta %>%
  filter(!(grade %in% c("GX")))
table(meta$grade)

#——Stage——————————————————————————————————
table(meta$stage)

meta <- meta %>%
  filter(!(stage %in% c("not reported", "stage x")))
meta$stage <- gsub("stage ", "", meta$stage)
meta$stage <- gsub("[abc]", "", meta$stage)

stage_mapping <- c("i" = 1, "ii" = 2, "iii" = 3, "iv" = 4)
meta$stage <- stage_mapping[meta$stage]
table(meta$stage)

#——T_Stage——————————————————————————————————
table(meta$T_stage)
meta$T_stage <- gsub("[abcd]", "", meta$T_stage)
table(meta$T_stage)

#——M_Stage——————————————————————————————————
table(meta$M_stage)
meta$M_stage <- gsub("\\s*\\(.*?\\)", "", meta$M_stage)
table(meta$M_stage)
meta$M_stage <- ifelse(meta$M_stage == "M1a", "M1", meta$M_stage)
table(meta$M_stage)
meta <- meta %>%
  filter(M_stage %in% c("M0", "M1"))
table(meta$M_stage)

#——N_Stage——————————————————————————————————
table(meta$N_stage)
meta$N_stage <- gsub("\\s*\\(.*?\\)", "", meta$N_stage)
meta$N_stage <- gsub("[abcmi]", "", meta$N_stage)
meta <- meta %>%
  filter(N_stage != "NX")
table(meta$N_stage)


#合并信息
surv_filtered <- surv[surv$sample %in% meta$ID,]

pd <- merge(surv_filtered,meta, by.x = "sample", by.y = "ID")


save(pd, file = "pd.rdata")
