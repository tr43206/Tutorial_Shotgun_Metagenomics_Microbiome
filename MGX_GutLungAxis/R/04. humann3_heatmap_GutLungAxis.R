#install.packages(c("pheatmap", "RColorBrewer"))

rm(list = ls())

library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

setwd('C:/Users/tr432/Desktop/Roh/_Study_Design/25.05.02~_Gut-Lung_Axis/Exp1+1-S/MGX/04.humann.out/pathabundance_relab_split')

#######################################################################################
# 1. 파일 불러오기

path_df <- read_excel('pathabundance_relab_unstratified_heatmap용.xlsx',
                      sheet = 'unstrat',
#                      header = TRUE,
#                      check.names = FALSE,
#                      sep = '\t'
                      )

res_df <- read.csv('maaslin3_Group_abund_BKY.csv',
                   header = TRUE,
                   check.names = FALSE
                   )

#######################################################################################
# 1. 컬럼 타입 정리
path_df$`# Pathway` <- trimws(as.character(path_df$`# Pathway`))  #Edit: `# Pathway` or `# Gene Family`
res_df$feature <- trimws(as.character(res_df$feature))
res_df$q_abund_ind_bky <- as.numeric(res_df$q_abund_ind_bky)

res_df <- res_df %>%
  arrange(q_abund_ind_bky)# %>%
#  slice_head(n = 15)


# 2. significant pathway만 추출
sig_pathways <- unique(res_df$feature[res_df$q_abund_ind_bky < 0.05# &
#                                        res_df$coef < 0
                                      ])


# 3. abundance 테이블에서 해당 pathway만 필터링
path_sig <- path_df[path_df$`# Pathway` %in% sig_pathways, ]


dim(path_sig)
head(path_sig[, 1:3])

#######################################################################################
# 4. matrix 만들기

mat <- as.matrix(path_sig[, -1])
rownames(mat) <- path_sig$`# Pathway`  #Edit: `# Pathway` or `# Gene Family`
mode(mat) <- "numeric"

# 5. 샘플 annotation
group <- data.frame(
  Group = c(rep('Control', 6), rep('VNAM', 9))
)
rownames(group) <- colnames(mat)

ann_colors <- list(
  Group = c(Control = "#FFC0CB", VNAM = "#87CEEB")
)

#######################################################################################
# 6. heatmap

p <- pheatmap(mat,
              scale = 'row',
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              annotation_col = group,
              annotation_colors = ann_colors,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              border_color = 'grey80',  #Edit: 'grey80' or FALSE
              fontsize_row = 9.5,
              fontsize_col = 8,
              filename = "humann_heatmap_sig_v3.png",
              width = 9,
              height = 5
              )
p

dev.off()

#######################################################################################
