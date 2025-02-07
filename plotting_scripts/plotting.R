library('ggplot2')
library('data.table')
library('patchwork')
library('ggrepel')
library('dplyr')
library(tidyr)
library(eulerr)
library(ggpubr)


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
  ggtitle("GO Size by Type") +
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

p
ggsave(p, filename="./../report/plots/goSizesByType.png", dpi=300, height=10, width = 6)

p <- ggplot(combined_diff_sizes, aes(x = Type, y = Size, color = Type)) +
  geom_boxplot(alpha = 0.8) +
  scale_y_log10() +
  xlab("Type") +
  ylab("Size") +
  ggtitle("Boxplot of Size Differences of Parents to their Children by Type") +
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
ggsave(p, filename="./../report/plots/goSizediffByType.png", dpi=300, , height=10, width = 6)


# min = 50, max = 500
go_out <- fread("./../go_out.out")
go_out[, shortest_path_to_a_true := NULL]
go_out[, Type := "GO"]
ens_out <- fread("./../ens_out.out")
ens_out[, shortest_path_to_a_true := NULL]
ens_out[, Type := "ENSEMBL"]

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
    axis.text = element_text(size = 14),
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
  ggtitle("Boxplot of Benjamini Hochberg Adjusted feh P-Values(min=50, max=500)") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
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
  ORA = sig_gos_go %>% filter(Method == "FEJ") %>% pull(term),
  GSEA = sig_gos_go %>% filter(Method == "KS") %>% pull(term)
)

sig_lists_ens <- list(
  ORA = sig_gos_ens %>% filter(Method == "FEJ") %>% pull(term),
  GSEA = sig_gos_ens %>% filter(Method == "KS") %>% pull(term)
)

p <- plot(euler(sig_lists_go), quantities = list(cex=2), labels=list(cex=2),fills = list(fill = c("#ea9d34","#eb6f92","#56949f")))
p
ggsave(p, filename="./../report/plots/goOverlap.png", dpi=300,  height = 10, width = 12)

p <- plot(euler(sig_lists_ens), quantities = list(cex=2), labels=list(cex=2),fills = list(fill = c("#ea9d34","#eb6f92","#56949f")))
p
ggsave(p, filename="./../report/plots/ensOverlap.png", dpi=300,  height = 10, width = 12)



go_scatt_log_fej  <- ggplot(go_out, aes(x = size, y = -log(fej_fdr))) +
  geom_point(alpha = 0.7, color = 'blue') +
  labs(
    title = "Scatter Plot of Gene Set Size vs FEJ-FDR (GO Mapping)",
    x = "Gene Set Size",
    y = "-log(FEJ-FDR)"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme for better visuals
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

go_scatt_log_ks  <- ggplot(go_out, aes(x = size, y = -log(ks_fdr))) +
  geom_point(alpha = 0.7, color = 'blue') +
  labs(
    title = "Scatter Plot of Gene Set Size vs KS-FDR (GO Mapping)",
    x = "Gene Set Size",
    y = "-log(KS-FDR)"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme for better visuals
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p <- go_scatt_log_fej|go_scatt_log_ks
ggsave(p, filename="./../report/plots/goScatt.png", dpi=300, height=10, width=15)

ens_scatt_log_fej <- ggplot(ens_out, aes(x = size, y = -log10(fej_fdr))) +
  geom_point(alpha = 0.7, color = 'blue') +
  labs(
    title = "Scatter Plot of Gene Set Size vs FEJ-FDR (ENSEMBL Mapping)",
    x = "Gene Set Size",
    y = "-log(FEJ-FDR)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ens_scatt_log_ks <- ggplot(ens_out, aes(x = size, y = -log10(ks_fdr))) +
  geom_point(alpha = 0.7, color = 'blue') +
  labs(
    title = "Scatter Plot of Gene Set Size vs KS-FDR (ENSEMBL Mapping)",
    x = "Gene Set Size",
    y = "-log(KS-FDR)"
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme for better visuals
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
p <- ens_scatt_log_fej| ens_scatt_log_ks
ggsave(p, filename="./../report/plots/ensScatt.png", dpi=300, height=10, width=15)
