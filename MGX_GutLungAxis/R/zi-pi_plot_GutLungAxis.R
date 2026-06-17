# install.packages("MetaNet")
# install.packages("pcutils")
# install.packages("ggrepel")
# install.packages("patchwork")

rm(list = ls())

library(MetaNet)
library(pcutils)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(igraph)

setwd('C:/Users/tr432/Desktop/Roh/_Study_Design/25.05.02~_Gut-Lung_Axis/Exp1+1-S/MGX/04.metaphlan4.out')

# =========================
# Input
# =========================
abund <- read.delim('only_species_FGB+GGB 제거.txt', 
#                    skip = 1, 
                    row.names = 1, check.names = FALSE)
meta  <- read.delim('C:/Users/tr432/Downloads/filter_GutLungAxis_Control/maaslin3_RUN45_Metadata_GutLungAxis.txt', check.names = FALSE)

# abundance를 numeric으로 강제
abund[] <- lapply(abund, function(x) as.numeric(as.character(x)))

######################################################################################
# Function 1: network + Zi/Pi 계산
######################################################################################
run_pizi <- function(group_name,
                     min_prev = 0.20,     # prevalence cutoff
                     min_mean = 1e-4,     # mean relative abundance cutoff
                     r_cut = 0.50,        # correlation cutoff
                     p_cut = 0.05,        # p-value cutoff
                     p_adjust = FALSE,    # TRUE면 BH 보정 사용
                     out_prefix = group_name) {
  
  # -------------------------
  # 1) group + Week2 + MGX 샘플만 선택
  # -------------------------
  samples <- meta %>%
    filter(Group == group_name,
           SamplingWeek == "Week2",
           !is.na(MGX_Reads)) %>%
    pull("#SampleID")
  
  samples <- intersect(samples, colnames(abund))
  
  cat("\n=============================\n")
  cat("Group:", group_name, "\n")
  cat("Samples:", paste(samples, collapse = ", "), "\n")
  cat("n =", length(samples), "\n")
  
  if (length(samples) < 4) {
    stop(paste0(group_name, ": usable sample 수가 너무 적습니다."))
  }
  
  # -------------------------
  # 2) abundance matrix 생성
  # 이미 relative abundance 파일이므로 재계산 안 함
  # -------------------------
  mat <- abund[, samples, drop = FALSE]
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  # -------------------------
  # 3) taxon filtering
  # -------------------------
  keep <- rowMeans(mat > 0, na.rm = TRUE) >= min_prev &
    rowMeans(mat, na.rm = TRUE) >= min_mean &
    rowSums(mat, na.rm = TRUE) > 0 &
    apply(mat, 1, sd, na.rm = TRUE) > 0
  
  mat <- mat[keep, , drop = FALSE]
  
  cat("Taxa after filtering:", nrow(mat), "\n")
  
  if (nrow(mat) < 5) {
    stop(paste0(group_name, ": filtering 후 taxa 수가 너무 적습니다."))
  }
  
  # -------------------------
  # 4) MetaNet input 형태로 변환
  # row = sample, column = taxon
  # -------------------------
  totu <- as.data.frame(t(mat))
  
  # -------------------------
  # 5) correlation 계산
  # -------------------------
  if (p_adjust) {
    corr <- c_net_calculate(
      totu,
      method = "spearman",
      p.adjust.method = "BH"
    )
  } else {
    corr <- c_net_calculate(
      totu,
      method = "spearman",
      p.adjust.method = NULL
    )
  }
  
  # -------------------------
  # 6) network 생성
  # -------------------------
  net <- c_net_build(
    corr,
    r_threshold = r_cut,
    p_threshold = p_cut,
    use_p_adj = p_adjust,
    delete_single = TRUE
  )
  
  cat("Nodes:", igraph::vcount(net), "\n")
  cat("Edges:", igraph::ecount(net), "\n")
  
  if (igraph::ecount(net) == 0) {
    stop(paste0(group_name, ": edge가 0개입니다. threshold를 완화해야 합니다."))
  }
  
  # -------------------------
  # 7) module detection
  # -------------------------
  net_mod <- module_detect(
    net,
    method = "cluster_fast_greedy"
  )
  
  # -------------------------
  # 8) Zi / Pi 계산
  # -------------------------
  net_zp <- zp_analyse(net_mod)
  
  # -------------------------
  # 9) 결과 table 정리
  # -------------------------
  zp_table <- get_v(net_zp)
  
  # 실제 taxon name 찾기
  if ("name" %in% colnames(zp_table)) {
    zp_table$Taxon <- zp_table$name
  } else if ("label" %in% colnames(zp_table)) {
    zp_table$Taxon <- zp_table$label
  } else if ("v_name" %in% colnames(zp_table)) {
    zp_table$Taxon <- zp_table$v_name
  } else {
    zp_table$Taxon <- igraph::V(net_zp)$name
  }
  
  # degree 추가
  deg_vec <- igraph::degree(net_mod)
  zp_table$degree <- deg_vec[zp_table$Taxon]
  
  # 역할 재분류
  zp_table$role <- dplyr::case_when(
    zp_table$Zi > 2.5  & zp_table$Pi > 0.62  ~ "Network hubs",
    zp_table$Zi > 2.5  & zp_table$Pi <= 0.62 ~ "Module hubs",
    zp_table$Zi <= 2.5 & zp_table$Pi > 0.62  ~ "Connectors",
    TRUE ~ "Peripherals"
  )
  
  # 저장
  write.csv(zp_table,
            paste0(out_prefix, "_ZiPi_table.csv"),
            row.names = FALSE)
  
  return(list(
    corr = corr,
    net = net,
    net_mod = net_mod,
    net_zp = net_zp,
    zp_table = zp_table
  ))
}

######################################################################################
# Function 2: 발표사진 스타일 Pi-Zi scatter plot
######################################################################################
plot_pizi <- function(zp_table,
                      title_text = "Control",
                      point_color = "#6FA8DC",
                      label_color = "firebrick",
                      label_size = 4,
                      base_size = 15) {
  
  df <- zp_table
  
  # 혹시 role 열이 factor면 문자형으로
  df$role <- as.character(df$role)
  
  # 라벨 붙일 taxon: hub / connector만
  lab_df <- df %>%
    filter(role != "Peripherals")
  
  # degree 처리
  if (all(is.na(df$degree)) || length(unique(na.omit(df$degree))) <= 1) {
    df$degree_plot <- rep(3, nrow(df))
    use_size_scale <- FALSE
  } else {
    df$degree_plot <- df$degree
    use_size_scale <- TRUE
  }
  
  p <- ggplot(df, aes(x = Pi, y = Zi)) +
    
    geom_vline(xintercept = 0.62, linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = 2.5,  linetype = "dashed", linewidth = 0.8) +
    
    geom_point(aes(shape = role, size = degree_plot),
               color = point_color,
               alpha = 0.85) +
    
    # lab_df가 0행이어도 문제없이 동작
    geom_text_repel(
      data = lab_df,
      aes(label = Taxon),
      color = label_color,
      fontface = "bold",
      size = label_size,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey50"
    ) +
    
    annotate("text", x = 0.16, y = 4.5, label = "Module hubs", fontface = "bold", size = 5) +
    annotate("text", x = 0.82, y = 4.5, label = "Network hubs", fontface = "bold", size = 5) +
    annotate("text", x = 0.16, y = -1.7, label = "Peripherals", fontface = "bold", size = 5) +
    annotate("text", x = 0.82, y = -1.7, label = "Connectors", fontface = "bold", size = 5) +
    
    scale_shape_manual(
      values = c(
        "Peripherals"  = 16,
        "Module hubs"  = 17,
        "Connectors"   = 15,
        "Network hubs" = 18
      ),
      breaks = c("Module hubs", "Network hubs", "Peripherals", "Connectors"),
      drop = FALSE
    ) +
    
    scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(-2, 5), breaks = seq(-2, 5, 1)) +
    
    labs(
      title = title_text,
      x = "Among-module connectivity (Pi)",
      y = "Within-module degree (Zi)",
      shape = NULL,
      size = NULL
    ) +
    
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 13),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1),
      axis.ticks.length = unit(0.2, "cm"),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.key.height = unit(0.7, "cm")
    )
  
  if (use_size_scale) {
    p <- p + scale_size_continuous(range = c(2.5, 8))
  } else {
    p <- p + scale_size_identity()
  }
  
  return(p)
}

######################################################################################
# Run
######################################################################################
res_Con <- run_pizi(
  group_name = "Control",
  min_prev   = 0.00,
  min_mean   = 0,
  r_cut      = 0.50,
  p_cut      = 0.05,
  p_adjust   = FALSE,
  out_prefix = "Control"
)

res_VNAM <- run_pizi(
  group_name = "VNAM",
  min_prev   = 0.00,
  min_mean   = 0,
  r_cut      = 0.50,
  p_cut      = 0.05,
  p_adjust   = FALSE,
  out_prefix = "VNAM"
)

######################################################################################
# Plot
######################################################################################
p_con  <- plot_pizi(res_Con$zp_table,  title_text = "Control")
p_vnam <- plot_pizi(res_VNAM$zp_table, title_text = "VNAM")

# 개별 확인
p_con
p_vnam

res_Con$zp_table[
  res_Con$zp_table$name == 'Phocaeicola_vulgatus',
]

res_Con$zp_table[
  res_Con$zp_table$name == 'Enterobacter_jensenii',
]

res_Con$zp_table[
  res_Con$zp_table$role == 'Module hubs',
]

res_Con$zp_table[
  res_Con$zp_table$role == 'Connectors',
]

res_Con$zp_table[
  res_Con$zp_table$role == 'Network hubs',
]

# 저장
ggsave("Control_PiZi_plot.png", p_con, width = 6, height = 5, dpi = 500, bg = "white")
ggsave("VNAM_PiZi_plot.png",    p_vnam, width = 6, height = 5, dpi = 500, bg = "white")

# 두 패널 나란히
p_both <- p_con + p_vnam + plot_layout(ncol = 2)

p_both

ggsave("Control_VNAM_PiZi_plot_2panel.png",
       p_both,
       width = 12,
       height = 5.5,
       dpi = 500,
       bg = "white")