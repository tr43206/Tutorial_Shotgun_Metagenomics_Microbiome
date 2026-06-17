# =========================================================
# MGX alpha / beta diversity with metadata
# abundance: species_merged_frac(1).txt
# metadata : maaslin3_RUN45_Metadata_GutLungAxis(1).txt
# =========================================================

rm(list = ls())

library(tidyverse)
library(vegan)
library(ggplot2)

# -----------------------------
# 1) 파일 경로
# -----------------------------

setwd('C:/Users/tr432/Desktop/Roh/_Study_Design/25.05.02~_Gut-Lung_Axis/Exp1+1-S/MGX/04.kraken2.out/bracken/species')

abund_file <- "species_merged_frac.txt"
meta_file  <- "C:/Users/tr432/Downloads/filter_GutLungAxis_Control/maaslin3_RUN45_Metadata_GutLungAxis.txt"

# -----------------------------
# 2) abundance 불러오기
#    첫 열: species name
#    나머지: sample columns
# -----------------------------
abund_raw <- read.delim(
  abund_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# species를 rowname으로
abund <- abund_raw %>%
  column_to_rownames(var = "name")

# numeric 변환
abund[] <- lapply(abund, as.numeric)

# species x sample  -> sample x species
abund_t <- as.data.frame(t(as.matrix(abund)), check.names = FALSE)
abund_t$SampleID <- rownames(abund_t)

# -----------------------------
# 3) metadata 불러오기
#    #SampleID 열 이름 정리
# -----------------------------
meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(meta)[colnames(meta) == "#SampleID"] <- "SampleID"

# fecal sample만 사용
meta <- meta %>%
  filter(SampleType == "Fecal")

# 필요한 열만 우선 보존
meta_sub <- meta %>%
  select(
    SampleID,
    Group,
    SamplingWeek,
    Subject_ID,
    ExpNo,
    MouseNo2,
    CageNo,
    Age,
    Sex,
    Initial_b.w,
    b.w,
    BM_wbc,
    Spleen
  )

# -----------------------------
# 4) 공통 sample만 추출
# -----------------------------
common_samples <- intersect(abund_t$SampleID, meta_sub$SampleID)

abund_t2 <- abund_t %>%
  filter(SampleID %in% common_samples) %>%
  arrange(SampleID)

meta_sub2 <- meta_sub %>%
  filter(SampleID %in% common_samples) %>%
  arrange(SampleID)

# 순서 일치 확인
stopifnot(all(abund_t2$SampleID == meta_sub2$SampleID))

# rowname 제거 후 SampleID를 rowname으로 지정
rownames(abund_t2) <- NULL

otu_mat <- abund_t2
rownames(otu_mat) <- otu_mat$SampleID
otu_mat$SampleID <- NULL
otu_mat <- as.data.frame(otu_mat, check.names = FALSE)

# Group / Week factor 정리
meta_sub2$Group <- factor(meta_sub2$Group, levels = c("Control", "VNAM"))
meta_sub2$SamplingWeek <- factor(meta_sub2$SamplingWeek,
                                 levels = c("Week0", "Week1", "Week2"))

# -----------------------------
# 5) Alpha diversity 계산
# -----------------------------
alpha_df <- data.frame(
  SampleID  = rownames(otu_mat),
  Observed  = specnumber(otu_mat),
  Shannon   = diversity(otu_mat, index = "shannon"),
  Simpson   = diversity(otu_mat, index = "simpson")
) %>%
  left_join(meta_sub2, by = "SampleID")

write.csv(alpha_df, "alpha_diversity_results.csv", row.names = FALSE)

######################################################################################
# -----------------------------
# 10) Beta diversity
# -----------------------------
# -----------------------------
# Bray-Curtis
# -----------------------------
bray_dist <- vegdist(otu_mat, method = "bray")
pcoa_bray <- cmdscale(bray_dist, eig = TRUE, k = 2)

pcoa_bray_df <- data.frame(
  SampleID = rownames(otu_mat),
  PC1 = pcoa_bray$points[, 1],
  PC2 = pcoa_bray$points[, 2]
) %>%
  left_join(meta_sub2, by = "SampleID")

eig_bray <- pcoa_bray$eig
pos_bray <- eig_bray[eig_bray > 0]
var_bray <- round(100 * pos_bray / sum(pos_bray), 4)

write.csv(pcoa_bray_df, "pcoa_bray_coordinates.csv", row.names = FALSE)
write.csv(
  data.frame(PC1 = var_bray[1], PC2 = var_bray[2]),
  "pcoa_bray_variance.csv",
  row.names = FALSE
)

permanova_bray <- adonis2(bray_dist ~ Group, data = meta_sub2, permutations = 999)
write.csv(
  data.frame(
    Method = "Bray-Curtis",
    F = permanova_bray$F[1],
    R2 = permanova_bray$R2[1],
    Pvalue = permanova_bray$`Pr(>F)`[1],
    N = nrow(meta_sub2)
  ),
  "permanova_bray.csv",
  row.names = FALSE
)


# -----------------------------
# Jaccard (presence/absence)
# -----------------------------
otu_pa <- ifelse(otu_mat > 0, 1, 0)
jaccard_dist <- vegdist(otu_pa, method = "jaccard", binary = TRUE)

pcoa_jaccard <- cmdscale(jaccard_dist, eig = TRUE, k = 2)

pcoa_jaccard_df <- data.frame(
  SampleID = rownames(otu_mat),
  PC1 = pcoa_jaccard$points[, 1],
  PC2 = pcoa_jaccard$points[, 2]
) %>%
  left_join(meta_sub2, by = "SampleID")

eig_jaccard <- pcoa_jaccard$eig
pos_jaccard <- eig_jaccard[eig_jaccard > 0]
var_jaccard <- round(100 * pos_jaccard / sum(pos_jaccard), 4)

write.csv(pcoa_jaccard_df, "pcoa_jaccard_coordinates.csv", row.names = FALSE)
write.csv(
  data.frame(PC1 = var_jaccard[1], PC2 = var_jaccard[2]),
  "pcoa_jaccard_variance.csv",
  row.names = FALSE
)

permanova_jaccard <- adonis2(jaccard_dist ~ Group, data = meta_sub2, permutations = 999)
write.csv(
  data.frame(
    Method = "Jaccard",
    F = permanova_jaccard$F[1],
    R2 = permanova_jaccard$R2[1],
    Pvalue = permanova_jaccard$`Pr(>F)`[1],
    N = nrow(meta_sub2)
  ),
  "permanova_jaccard.csv",
  row.names = FALSE
)
