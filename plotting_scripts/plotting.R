library('ggplot2')
library('data.table') 
library('patchwork')
library('ggrepel')
library('dplyr')
library(tidyr)
library(eulerr)
library(ggpubr)
library(VennDiagram)
groundTruth  <- readLines("./groundTruthMinMax.txt")


go_sizes <- fread("./../go_goSizes.tsv", col.names = c("GOid", "Size"))
go_sizes[, Type:="GO"]
go_diff_sizes <- fread("./../go_goDiff.tsv", col.names = c("Size"))
go_diff_sizes[, Type:="GO"]

ens_sizes <- fread("./../ensembl_goSizes.tsv", col.names = c("GOid", "Size"))
ens_sizes[, Type:="ENSEMBL"]
ens_diff_sizes <- fread("./../ensembl_goDiff.tsv", col.names = c("Size"))
ens_diff_sizes[, Type:="ENSEMBL"]

combined_sizes <- rbind(go_sizes, ens_sizes)
combined_diff_sizes  <- rbind(go_diff_sizes, ens_diff_sizes)


p <- ggplot(combined_sizes, aes(x = Type, y = Size, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("Size") +
  ggtitle("GO Size by Mapping Type") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  )

p
# ggsave(p, filename="./../report/plots/goSizes.png", dpi=300, height=10, width = 6)

q <- ggplot(combined_diff_sizes, aes(x = Type, y = Size, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("Size") +
  ggtitle("Size Differences of Parents to their Children \n by Mapping Type") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  )
# ggsave(p, filename="./../report/plots/goSizediffByType.png", dpi=300, , height=10, width = 6)
plot  <- p | q
plot
ggsave(plot, filename="./../report/plots/goSizes.png", dpi=300, height=10, width = 20)


# min = 50, max = 500
go_out <- fread("./../go_out.out")
go_out[, shortest_path_to_a_true := NULL]
go_out[, Type := "GO"]
ens_out <- fread("./../ens_out.out")
ens_out[, shortest_path_to_a_true := NULL]
ens_out[, Type := "ENSEMBL"]

sig_fej  <- go_out[fej_fdr <= 0.05, ]

sig_fej[, term]
gProfiler <- readLines("./gProfiler_BP.txt")
gProfiler

groundTruth  <- readLines("./groundTruthMinMax.txt")


compare   <- list(
  gProfiler = unique(gProfiler),
  JAR_ORA_GO = unique(sig_fej[, term]),
  SoT = unique(groundTruth)
)


venn.plot <- venn.diagram(
  x = compare,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.1
)

png("./../report/plots/go_mappingCompgProfiler.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()

ens_fej  <- ens_out[fej_fdr <= 0.05, ]

ens_fej[, term]
gProfiler <- readLines("./gProfiler_BP.txt")
gProfiler

a <- intersect(gProfiler, groundTruth)
# [1] "GO:0051051" "GO:0040013" "GO:0048511" "GO:0045926" "GO:1903828"
# [6] "GO:1904950" "GO:0043433" "GO:0098754" "GO:2000242"
intersect(a, unique(ens_fej[, term]))
# [1] "GO:0051051" "GO:0040013" "GO:0048511" "GO:0045926" "GO:1903828"
# [6] "GO:1904950" "GO:0043433"

# → "GO:0098754" "GO:2000242" not found by our JAR


compare   <- list(
  gProfiler = unique(gProfiler),
  JAR_ORA_ENS = unique(ens_fej[, term]),
  SoT = unique(groundTruth)
)

venn.plot <- venn.diagram(
  x = compare,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.1
)


png("./../report/plots/ens_mappingCompgProfiler.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()


combined_out  <- rbind(go_out, ens_out)
p <- ggplot(combined_out, aes(x = Type, y = ks_stat, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("Enrichment Score") +
  ggtitle("Boxplot of Enrichment Score by Type") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  )

ggsave(p, filename="./../report/plots/esByType.png", dpi=300, height=10, width = 6)

top5 <- combined_out %>%
  group_by(Type) %>%
  top_n(-5, ks_fdr) %>%
  ungroup()


p <- ggplot(combined_out, aes(x = Type, y = ks_fdr, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("Enrichment Score") +
  ggtitle("Boxplot of Benjamini Hochberg Adjusted P-Values of Enrichmen Scores (min=50, max=500)") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(x = 0.8, y = 0.05, label = "alpha = 0.05"), color = "red", vjust = -1, size = 5) +
  geom_text_repel(data = top5, aes(x = Type, y = ks_fdr, label = term),
                  color = "black", size = 6, max.overlaps = Inf, 
                  box.padding = 1.7, point.padding = 0.1, force=2)
  p
ggsave(p, filename="./../report/plots/BHBoxplotFDRKS.png", dpi=300, height=10, width = 12)

top5 <- combined_out %>%
  group_by(Type) %>%
  top_n(-5, fej_fdr) %>%
  ungroup()


p <- ggplot(combined_out, aes(x = Type, y = fej_fdr, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("Enrichment Score") +
  ggtitle("Boxplot of Benjamini Hochberg Adjusted fej P-Values(min=50, max=500)") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(x = 0.8, y = 0.05, label = "alpha = 0.05"), color = "red", vjust = -1, size = 5) +
  geom_text_repel(data = top5, aes(x = Type, y = fej_fdr, label = term),
                  color = "black", size = 6, max.overlaps = Inf, 
                  box.padding = 1.7, point.padding = 0.1, force=2)
ggsave(p, filename="./../report/plots/BHBoxplotFDRFEJ.png", dpi=300, height=10, width = 12)

sig_threshold <- 0.05

sig_gos_go <- go_out %>%
  mutate(
    KS_Significant = ifelse(ks_fdr < sig_threshold, "KS", NA),
    FEJ_Significant = ifelse(fej_fdr < sig_threshold, "FEJ", NA)
  ) %>%
  pivot_longer(cols = c(KS_Significant, FEJ_Significant), names_to = "Source", values_to = "Method") %>%
  filter(!is.na(Method)) %>%
  select(term, Method)

sig_gos_ens <- ens_out %>%
  mutate(
    KS_Significant = ifelse(ks_fdr < sig_threshold, "KS", NA),
    FEJ_Significant = ifelse(fej_fdr < sig_threshold, "FEJ", NA)
  ) %>%
  pivot_longer(cols = c(KS_Significant, FEJ_Significant), names_to = "Source", values_to = "Method") %>%
  filter(!is.na(Method)) %>%
  select(term, Method)

sig_lists_go <- list(
  JAR_ORA_GO = sig_gos_go %>% filter(Method == "FEJ") %>% pull(term),
  JAR_GSEA_GO = sig_gos_go %>% filter(Method == "KS") %>% pull(term),
  SoT  = groundTruth
)

sig_lists_ens <- list(
  JAR_ORA_ENS = sig_gos_ens %>% filter(Method == "FEJ") %>% pull(term),
  JAR_GSEA_ENS = sig_gos_ens %>% filter(Method == "KS") %>% pull(term),
  SoT  = groundTruth
)

a <- intersect(sig_gos_ens %>% filter(Method == "KS") %>% pull(term), groundTruth)
#  [1] "GO:0002683" "GO:0032845" "GO:0040013" "GO:0043433" "GO:0043901"
#  [6] "GO:0045926" "GO:0048511" "GO:0051051" "GO:0099531" "GO:1903828"
# [11] "GO:1904950"
groundTruth
#  [1] "GO:0099531" "GO:0098754" "GO:0043433" "GO:1903828" "GO:2000242"
#  [6] "GO:1904950" "GO:0043901" "GO:0048511" "GO:0040013" "GO:0051051"
# [11] "GO:0045926" "GO:0002683" "GO:0032845"
setdiff(groundTruth, a)
# → [1] "GO:0098754" "GO:2000242" missing

venn.plot <- venn.diagram(
  x = sig_lists_go,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.11
)

png("./../report/plots/goOverlap.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()

venn.plot <- venn.diagram(
  x = sig_lists_ens,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.11
)

png("./../report/plots/ensOverlap.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()

# comparing GSEA and ORA

sig_lists_ora <- list(
  JAR_ORA_GO = sig_gos_go %>% filter(Method == "FEJ") %>% pull(term),
  JAR_ORA_ENS = sig_gos_ens %>% filter(Method == "FEJ") %>% pull(term),
  SoT  = groundTruth
)

sig_lists_gsea <- list(
  JAR_GSEA_GO = sig_gos_go %>% filter(Method == "KS") %>% pull(term),
  JAR_GSEA_ENS = sig_gos_ens %>% filter(Method == "KS") %>% pull(term),
  SoT  = groundTruth
)

venn.plot <- venn.diagram(
  x = sig_lists_gsea,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.11
)

png("./../report/plots/gseaOverlap.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()

venn.plot <- venn.diagram(
  x = sig_lists_ora,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.11
)

png("./../report/plots/oraOverlap.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()


go_scatt_log_fej  <- ggplot(go_out, aes(x = size, y = -log(fej_fdr))) +
  geom_point(alpha = 0.7, color = '#56949f') +
  labs(
    title = "Scatter Plot of Gene Set Size \n vs FEJ-padj (GO Mapping)",
    x = "Gene Set Size",
    y = "-log(FEJ-padj)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) 

go_scatt_log_ks  <- ggplot(go_out, aes(x = size, y = -log(ks_fdr))) +
  geom_point(alpha = 0.7, color = '#56949f') +
  labs(
    title = "Scatter Plot of Gene Set Size \n vs KS-padj (GO Mapping)",
    x = "Gene Set Size",
    y = "-log(KS-padj)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) 

p <- go_scatt_log_fej|go_scatt_log_ks
p
ggsave(p, filename="./../report/plots/goScatt.png", dpi=300, height=10, width=15)

ens_scatt_log_fej <- ggplot(ens_out, aes(x = size, y = -log10(fej_fdr))) +
  geom_point(alpha = 0.7, color = '#eb6f92') +
  labs(
    title = "Scatter Plot of Gene Set Size \nvs FEJ-padj (ENSEMBL Mapping)",
    x = "Gene Set Size",
    y = "-log(FEJ-padj)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) 

ens_scatt_log_ks <- ggplot(ens_out, aes(x = size, y = -log10(ks_fdr))) +
  geom_point(alpha = 0.7, color = '#eb6f92') +
  labs(
    title = "Scatter Plot of Gene Set Size \n vs KS-padj (ENSEMBL Mapping)",
    x = "Gene Set Size",
    y = "-log(KS-padj)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) 
p <- ens_scatt_log_fej| ens_scatt_log_ks
ggsave(p, filename="./../report/plots/ensScatt.png", dpi=300, height=10, width=15)




combined_out  <- rbind(go_out, ens_out)
p <- ggplot(combined_out, aes(x = Type, y = ks_stat, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("Enrichment Score") +
  ggtitle("Boxplot of Enrichment Score by Type") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  )

ggsave(p, filename="./../report/plots/esByType.png", dpi=300, height=10, width = 6)



# top5 <- combined_out %>%
#   group_by(Type) %>%
#   top_n(-5, ks_fdr) %>%
#   ungroup() %>%
#   mutate(text_color = ifelse(term %in% groundTruth, "blue", "black"))
#
# p <- ggplot(combined_out, aes(x = Type, y = ks_fdr, color = Type)) +
#   geom_boxplot(alpha = 0.8) +
#   scale_y_log10() +
#   xlab("Type") +
#   ylab("Enrichment Score") +
#   ggtitle("Benjamini Hochberg Adjusted P-Values of Enrichment Scores") +
#   theme_bw(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 14), 
#     legend.text = element_text(size = 12), 
#     panel.grid.major = element_line(color = "gray", size = 0.5), 
#     panel.grid.minor = element_line(color = "lightgray", size = 0.25)
#   ) +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
#   geom_text(aes(x = 0.8, y = 0.05, label = "alpha = 0.05"), color = "red", vjust = -1, size = 5) +
#   geom_text_repel(data = top5, aes(x = Type, y = ks_fdr, label = term),
#                   color = top5$text_color, size = 6, max.overlaps = Inf, 
#                   box.padding = 1.7, point.padding = 0.1, force = 2)
# ggsave(p, filename = "./../report/plots/BHBoxplotFDRKS.png", dpi = 300, height = 10, width = 12)
#
# top5 <- combined_out %>%
#   group_by(Type) %>%
#   top_n(-5, fej_fdr) %>%
#   ungroup() %>%
#   mutate(text_color = ifelse(term %in% groundTruth, "blue", "black"))
#
# p <- ggplot(combined_out, aes(x = Type, y = fej_fdr, color = Type)) +
#   geom_boxplot(alpha = 0.8) +
#   scale_y_log10() +
#   xlab("Type") +
#   ylab("Enrichment Score") +
#   ggtitle("Benjamini Hochberg Adjusted fej P-Values") +
#   theme_bw(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 14), 
#     legend.text = element_text(size = 12), 
#     panel.grid.major = element_line(color = "gray", size = 0.5), 
#     panel.grid.minor = element_line(color = "lightgray", size = 0.25)
#   ) +
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
#   geom_text(aes(x = 0.8, y = 0.05, label = "alpha = 0.05"), color = "red", vjust = -1, size = 5) +
#   geom_text_repel(data = top5, aes(x = Type, y = fej_fdr, label = term),
#                   color = top5$text_color, size = 6, max.overlaps = Inf, 
#                   box.padding = 1.7, point.padding = 0.1, force = 2)
#
# ggsave(p, filename = "./../report/plots/BHBoxplotFDRFEJ.png", dpi = 300, height = 10, width = 12)

top5 <- combined_out %>%
  group_by(Type) %>%
  top_n(-5, ks_fdr) %>%
  ungroup() %>%
  mutate(text_color = ifelse(term %in% groundTruth, "blue", "black"))

p <- ggplot(combined_out, aes(x = Type, y = ks_fdr, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("padj") +
  ggtitle("Benjamini Hochberg Adjusted P-Values of Enrichment Scores") +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(x = 0.8, y = 0.05, label = "alpha = 0.05"), color = "red", vjust = -1, size = 6) +
  geom_text_repel(data = top5, aes(x = Type, y = ks_fdr, label = term),
                  color = top5$text_color, size = 7, max.overlaps = Inf, 
                  box.padding = 2, point.padding = 0.6, force = 22)

ggsave(p, filename = "./../report/plots/BHBoxplotFDRKS.png", dpi = 300, height = 10, width = 12)

top5 <- combined_out %>%
  group_by(Type) %>%
  top_n(-5, fej_fdr) %>%
  ungroup() %>%
  mutate(text_color = ifelse(term %in% groundTruth, "blue", "black"))

p <- ggplot(combined_out, aes(x = Type, y = fej_fdr, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("padj") +
  ggtitle("Benjamini Hochberg Adjusted fej P-Values") +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18), 
    panel.grid.major = element_line(color = "gray", size = 0.5), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25)
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(x = 0.8, y = 0.05, label = "alpha = 0.05"), color = "red", vjust = -1, size = 6) +
  geom_text_repel(data = top5, aes(x = Type, y = fej_fdr, label = term),
                  color = top5$text_color, size = 7, max.overlaps = Inf, 
                  box.padding = 1.7, point.padding = 0.1, force = 2)

ggsave(p, filename = "./../report/plots/BHBoxplotFDRFEJ.png", dpi = 300, height = 10, width = 12)


sig_ks  <- go_out[ks_fdr <= 0.05, ]

sig_ks[, term]
fgsea  <- fread("../gsea/gsea_results.tsv")
groundTruth  <- readLines("./groundTruthMinMax.txt")


compare   <- list(
  fgsea  = fgsea[, gs_exact_source],
  JAR_GSEA_GO = sig_ks[, term],
  SoT = groundTruth
)

venn.plot <- venn.diagram(
  x = compare,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.1
)

png("./../report/plots/go_mappingCompfgsea.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()

ens_ks  <- ens_out[ks_fdr <= 0.05, ]

ens_ks[, term]

compare   <- list(
  fgsea  = fgsea[, gs_exact_source],
  JAR_GSEA_ENS = ens_ks[, term],
  SoT = groundTruth
)

a <- intersect(fgsea[, gs_exact_source], groundTruth)
#  [1] "GO:0051051" "GO:0002683" "GO:0040013" "GO:1903828" "GO:0048511"
#  [6] "GO:0043433" "GO:0098754" "GO:0045926" "GO:1904950" "GO:2000242"
intersect(a, unique(ens_ks[, term]))
# [1] "GO:0051051" "GO:0002683" "GO:0040013" "GO:1903828" "GO:0048511"
# [6] "GO:0043433" "GO:0045926" "GO:1904950"

# → "GO:0098754" "GO:2000242" missing

venn.plot <- venn.diagram(
  x = compare,
  filename = NULL,
  col = "black",
  fill = c("#ea9d34", "#eb6f92", "#56949f"),
  cat.col = "black",
  cat.cex = 3,
  cex = 2,
  margin = 0.11
)

png("./../report/plots/ens_mappingCompgfgsea.png", width = 10, height = 10, units = "in", res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()
