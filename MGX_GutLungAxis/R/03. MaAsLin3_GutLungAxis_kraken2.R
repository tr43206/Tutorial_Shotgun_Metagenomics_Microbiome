rm(list = ls())

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biobakery/maaslin3")
library(maaslin3)
library(mutoss)
library(dplyr)
library(ggplot2)

#####################################################################################
### Microbiome association detection with MaAsLin 3

setwd("C:/Users/tr432/Desktop/Roh/_Study_Design/25.05.02~_Gut-Lung_Axis/Exp1+1-S/MGX/04.kraken2.out/bracken/family")

### 상대적 풍부도 데이터
taxa_table <- read.csv("family_merged_frac.txt", sep = '\t', row.names = 1,
#                       skip = 1,
                       check.names = FALSE, comment.char = "")
taxa_table <- t(taxa_table)


metadata <- read.csv("C:/Users/tr432/Downloads/filter_GutLungAxis_Control/maaslin3_RUN45_Metadata_GutLungAxis.txt", sep = '\t', row.names = 1)


##**## (예시)
# meta$Group <- factor(meta$Group, levels = c("Control", "VNAM"))
# meta$Sex   <- factor(meta$Sex)        # 범주형 공변량 예시
# meta$Age   <- as.numeric(meta$Age)    # 연속형 공변량 예시
##**##

metadata$Group <- factor(metadata$Group, levels = c('Control', 'VNAM'))
#metadata$BM_wbc   <- as.numeric(metadata$BM_wbc)
#metadata$b.w   <- as.numeric(metadata$b.w)
metadata$ExpNo   <- factor(metadata$ExpNo)
metadata$CageNo   <- factor(metadata$CageNo)
metadata$Initial_b.w   <- as.numeric(metadata$Initial_b.w)
metadata$SequencingDepth   <- as.numeric(metadata$SequencingDepth)
metadata$MGX_Reads   <- as.numeric(metadata$MGX_Reads)

colnames(taxa_table) <- gsub("\\.", "-", colnames(taxa_table))

taxa_table[1:5, 1:5]
metadata[1:5, 1:5]

#####################################################################################

set.seed(1)
fit_out <- maaslin3(
  input_data    = taxa_table,
  input_metadata = metadata,
  output        = "maaslin3_Group+Initial_b.w+MGX_Reads_BH",  #Edit
  formula       = "~ Group + Initial_b.w + MGX_Reads",  #Edit
  normalization = "TSS",
  transform     = "LOG",  ## "LOG"는 base 2 log / "PLOG"는 pseudo-log
#  zero_threshold = 0,      ## 26.01.02 추가 - transform = "PLOG" 기준
#  evaluate_only  = "abundance",   ## 26.01.02 추가 - transform = "PLOG" 기준
#  warn_prevalence = FALSE,  ## 26.01.02 추가 - transform = "PLOG" 기준
  correction    = "BH",
  augment       = TRUE,
  standardize   = TRUE,
  max_significance = 0.1,
  median_comparison_abundance  = TRUE,  ## 상대적 풍부도엔 써도되지만 절대 풍부도엔 쓰면 안됨
  median_comparison_prevalence = TRUE,
  cores        = 1
)

## * 샘플 수가 적은데 (약 20개 내외) CageNo 종류가 많아서 오류가 났을 것으로 예상_25.12.15

#####################################################################################

## Abundance - BKY FDR
res_abund <- fit_out$fit_data_abundance$results %>%
  dplyr::select(-qval_individual, -qval_joint)

## individual p-values에 대한 BKY
keep_abund_ind <- !is.na(res_abund$pval_individual)
bky_abund_ind  <- mutoss::two.stage(res_abund$pval_individual[keep_abund_ind], alpha = 0.05)
res_abund$q_abund_ind_bky <- NA_real_
res_abund$q_abund_ind_bky[keep_abund_ind] <- bky_abund_ind$adjPValues


##**##
## (선택) 만약 유의한게 없어서 오류가 난다면 아래 부분으로 대체 (BKY 보정 직접 보정)
pv <- suppressWarnings(as.numeric(res_abund$pval_individual))
keep <- is.finite(pv)

bky_qvalues <- function(p, alpha = 0.05) {
  p <- as.numeric(p)
  m <- length(p)
  o <- order(p)
  ro <- order(o)
  ps <- p[o]
  
  # stage 1: BH at alpha' = alpha/(1+alpha)
  alpha1 <- alpha / (1 + alpha)
  crit1  <- (1:m) * alpha1 / m
  idx1   <- which(ps <= crit1)
  r1     <- if (length(idx1) == 0) 0 else max(idx1)
  
  # m0-hat = m - r1  (BKY two-step estimate)
  m0 <- m - r1
  if (m0 < 1) m0 <- 1  # 극단 케이스(거의 없음) 안전장치
  
  # stage 2: BH-style adjusted p-values with m0 in denominator
  qs <- (m0 * ps) / (1:m)
  qs <- rev(cummin(rev(qs)))
  qs[qs > 1] <- 1
  
  qs[ro]
}

res_abund$q_abund_ind_bky <- NA_real_
if (any(keep)) {
  # 1) mutoss 시도
  out <- try(mutoss::two.stage(pv[keep], alpha = 0.05), silent = TRUE)
  
  adj <- NULL
  if (!inherits(out, "try-error")) {
    adj <- out[["adjPValues"]]
    if (!is.null(adj)) adj <- as.numeric(adj)
  }
  
  # 2) mutoss가 이상하면 직접 계산으로 폴백
  if (is.null(adj) || length(adj) != sum(keep)) {
    adj <- bky_qvalues(pv[keep], alpha = 0.05)
  }
  
  res_abund$q_abund_ind_bky[keep] <- adj
}
##**##

write.csv(res_abund, "maaslin3_Group+Initial_b.w+MGX_Reads_abund_BKY.csv", row.names = FALSE)  #Edit

#####################################################################################

## Prevalence - BKY FDR
res_prev <- fit_out$fit_data_prevalence$results %>%
  dplyr::select(-qval_individual, -qval_joint)

## individual p-values에 대한 BKY
keep_prev_ind <- !is.na(res_prev$pval_individual)
bky_prev_ind  <- mutoss::two.stage(res_prev$pval_individual[keep_prev_ind], alpha = 0.05)
res_prev$q_prev_ind_bky <- NA_real_
res_prev$q_prev_ind_bky[keep_prev_ind] <- bky_prev_ind$adjPValues


##**##
## (선택) 만약 유의한게 없어서 오류가 난다면 아래 부분으로 대체 (BKY 보정 직접 보정)
pv <- suppressWarnings(as.numeric(res_prev$pval_individual))
keep <- is.finite(pv)

bky_qvalues <- function(p, alpha = 0.05) {
  p <- as.numeric(p)
  m <- length(p)
  o <- order(p)
  ro <- order(o)
  ps <- p[o]

  # stage 1: BH at alpha' = alpha/(1+alpha)
  alpha1 <- alpha / (1 + alpha)
  crit1  <- (1:m) * alpha1 / m
  idx1   <- which(ps <= crit1)
  r1     <- if (length(idx1) == 0) 0 else max(idx1)

  # m0-hat = m - r1  (BKY two-step estimate)
  m0 <- m - r1
  if (m0 < 1) m0 <- 1  # 극단 케이스(거의 없음) 안전장치

  # stage 2: BH-style adjusted p-values with m0 in denominator
  qs <- (m0 * ps) / (1:m)
  qs <- rev(cummin(rev(qs)))
  qs[qs > 1] <- 1

  qs[ro]
}

res_prev$q_prev_ind_bky <- NA_real_
if (any(keep)) {
# 1) mutoss 시도
  out <- try(mutoss::two.stage(pv[keep], alpha = 0.05), silent = TRUE)

  adj <- NULL
  if (!inherits(out, "try-error")) {
    adj <- out[["adjPValues"]]
    if (!is.null(adj)) adj <- as.numeric(adj)
  }

# 2) mutoss가 이상하면 직접 계산으로 폴백
  if (is.null(adj) || length(adj) != sum(keep)) {
    adj <- bky_qvalues(pv[keep], alpha = 0.05)
  }

  res_prev$q_prev_ind_bky[keep] <- adj
}
##**##

write.csv(res_prev, "maaslin3_Group+Initial_b.w+MGX_Reads_prev_BKY.csv", row.names = FALSE)  #Edit

#####################################################################################
## 개별 plot 수정

plots <- maaslin3::maaslin_plot_results_from_output(
  output        = "maaslin3_Group+Initial_b.w_BH",
  metadata      = metadata,
  normalization = "TSS",
  transform     = "LOG",
  median_comparison_abundance  = TRUE,
  median_comparison_prevalence = TRUE,
  max_significance = 0.1
)

## 원하는 plot 선택. assoc는 직접 찾아야됨. 아마 균주명만 바꾸면 될듯
assoc <- plots$assocation_plots$Group$`d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola_A`
names(assoc)
p <- assoc[["linear"]]

## y축 라벨에서 수정
p2 <- p + ylab("Phocaeicola_A")

## 우측 상단 텍스트(annotate 레이어) 제거 -> 보통 마지막 layer가 annotation이라 마지막 1개만 빼면 됨
p2$layers <- p2$layers[-length(p2$layers)]

## 새 텍스트로 다시 추가 (ex. BKY q값 + (in full model) 제거)
p2 <- p2 +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.05, vjust = 1.1, size = 3,
           label = sprintf("q = %.3f\nCoef: %.2f",
                           0.01027221,  #Edit: BKY에서 얻은 q 값
                           -6.8270896))  #Edit: coef 값


ggsave("Group_Phocaeicola_A_linear_rev.png",
       p2, width = 3, height = 4, dpi = 500)

#####################################################################################
