library(dplyr)
library(ggplot2)
library(readr)

df <- read_csv("path/to/your/project/20251204/result_self-FHL1_average_expression_20251204.csv",
               locale = locale(encoding = "UTF-8"))   # or read_csv("your_file.csv")

# ========== Modify here as needed! ==========
# Which cells to plot? Change the vector below (order = plotting order)
#plot_these_cells <- c("B Cells", "Myeloid Cells", "NK Cells", "NK-T Cells",
                      #"Neutrophil", "Plasma Cells", "Platelets", "T Cells")
plot_these_cells <- c("B Cells", "Myeloid Cells", "NK Cells",
                      "Plasma Cells", "Platelets", "T Cells")
# Colors — must correspond to data groups
my_colors <- c(
  #"01 Health Control" = "#3B9AB2",
  "01 Mild"           = "#3B9AB2",
  "02 Survival"       = "#E69F00",
  "03 None-Survival"  = "#CC0000")

# Auto-filter data + set x-axis factor order
df_plot <- df %>%
  filter(cell_type %in% plot_these_cells) %>%
  mutate(cell_type = factor(cell_type, levels = plot_these_cells))

# Generate box-plot (same logic as before)
ggplot(df_plot, aes(x = cell_type, ymin = q05, lower = q25,
                    middle = median_expression, upper = q75, ymax = q95,
                    fill = group, color = group)) +
  geom_boxplot(stat = "identity", width = 0.7, alpha = 0.85, size = 0.9) +
  scale_fill_manual(values = my_colors, name = "Group") +
  scale_color_manual(values = my_colors, name = "Group") +
  labs(x = "Cell Type", y = "Expression Level", title = unique(df$gene)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave("selected_cells_boxplot.pdf", width = 10, height = 7)
