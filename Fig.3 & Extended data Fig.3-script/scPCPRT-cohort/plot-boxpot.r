# 加载必要包
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# 读取数据
data <- read.table("C:/Users/wupin/Desktop/R-data/data.txt", header = TRUE)

# 定义需要比较的组对
comparison_pairs <- list(
  c("Base_NR", "PD1_NR"),
  c("PD1_NR", "RTPD1_NR"),
  c("Base_R1", "PD1_R1"),
  c("PD1_R1", "RTPD1_R1"),
  c("Base_R2", "PD1_R2"),
  c("PD1_R2", "RTPD1_R2")
)

# 将数据从宽格式转为长格式
long_data <- data %>%
  pivot_longer(
    cols = 3:6,  # 第3-6列是打分数据
    names_to = "score",
    values_to = "value"
  ) %>%
  # 确保sample_type是因子
  mutate(sample_type = factor(sample_type))

# 筛选只包含需要比较的组别的数据
filtered_data <- long_data %>%
  filter(sample_type %in% unique(unlist(comparison_pairs)))

# 为每个打分创建单独的图
for (score_name in unique(filtered_data$score)) {
  
  # 筛选当前打分数据
  plot_data <- filtered_data %>% filter(score == score_name)
  
  # 创建基础图形
  p <- ggplot(plot_data, aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(width = 0.6, alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, aes(color = sample_type)) +
    labs(
      title = paste("Comparison of", score_name),
      x = "Sample Type",
      y = score_name
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"  # 移除图例，因为x轴已标明
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2")
  
  # 添加统计检验
  p <- p + stat_compare_means(
    comparisons = comparison_pairs,
    method = "wilcox.test",
    label = "p.format",  # 显示具体p值
    hide.ns = TRUE,      # 隐藏不显著的结果
    label.size = 3,
    step.increase = 0.1,
    tip.length = 0.01
  )
  
  # 保存图形
  ggsave(
    paste0('C:/Users/wupin/Desktop/',score_name, "_comparison.pdf"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}