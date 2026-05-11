library(xgboost)
library(SHAPforxgboost)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)



output_dir <- "F:/karst" 

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("已创建输出文件夹: %s\n", output_dir))
}


# 1. 数据加载与预处理
data_file <- file.path(output_dir, "Sample Data.csv")
df <- read_csv(data_file)
df <- na.omit(df)

if("LF" %in% names(df) && is.character(df$LF)) {
  df$LF <- as.numeric(as.factor(df$LF)) - 1
}

y <- df$ES

urban_gradient <- df$Urbanization

X <- df %>% select(-ES, -Urbanization)
X_matrix <- as.matrix(X)


# 2. 构建与训练 XGBoost 模型
dtrain <- xgb.DMatrix(data = X_matrix, label = y)

# 设置模型参数
params <- list(
  objective = "reg:squarederror", # 回归任务
  max_depth = 6,                  # 树的最大深度
  eta = 0.05,                     # 学习率
  subsample = 0.8,
  colsample_bytree = 0.8
)


set.seed(42)
cat("正在训练 XGBoost 模型...\n")
xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 150)

# 计算训练集 R 平方
preds <- predict(xgb_model, X_matrix)
r2 <- 1 - sum((y - preds)^2) / sum((y - mean(y))^2)
cat(sprintf("模型全局 R2 得分: %.4f\n", r2))



# 3. 计算 SHAP 值
cat("正在计算 SHAP 值...\n")

shap_results <- shap.values(xgb_model = xgb_model, X_train = X_matrix)
shap_score <- shap_results$shap_score 
shap_long <- shap.prep(xgb_model = xgb_model, X_train = X_matrix)


# 4. 核心可视化图表输出
p1 <- shap.plot.summary(shap_long) +
  theme_bw() +
  ggtitle("Global Feature Importance (SHAP Summary Plot)")

save_path_1 <- file.path(output_dir, "1_Global_SHAP_Summary.png")
ggsave(save_path_1, p1, width = 8, height = 6, dpi = 600)
cat(sprintf("已生成: %s\n", save_path_1))


shap_abs <- as.data.frame(abs(shap_score))

shap_abs$Urban_Ref <- urban_gradient


shap_abs$Urban_Bin <- ntile(shap_abs$Urban_Ref, 10)


gradient_summary <- shap_abs %>%
  select(-Urban_Ref) %>%
  group_by(Urban_Bin) %>%
  summarise(across(everything(), mean), .groups = 'drop')


gradient_long <- gradient_summary %>%
  pivot_longer(cols = -Urban_Bin, names_to = "Feature", values_to = "Mean_Abs_SHAP")


top_features <- gradient_long %>%
  group_by(Feature) %>%
  summarise(overall_mean = mean(Mean_Abs_SHAP)) %>%
  top_n(6, overall_mean) %>%
  pull(Feature)

gradient_filtered <- gradient_long %>% 
  filter(Feature %in% top_features) %>%
  mutate(Feature = factor(Feature, levels = top_features))

p2 <- ggplot(gradient_filtered, aes(x = Urban_Bin, y = Mean_Abs_SHAP, fill = Feature)) +
  geom_area(alpha = 0.85, color = "white", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2") + 

  labs(title = "Evolution of Dominant Drivers along Urbanization Level Gradient",
       x = "Urbanization Level Decile (1=Low -> 10=High)",
       y = "Mean |SHAP Value| (Impact Magnitude)") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(color = "black"), 
        legend.position = "right", 
        panel.grid.minor = element_blank())

save_path_2 <- file.path(output_dir, "2_Gradient_Drivers_Evolution.png")
ggsave(save_path_2, p2, width = 10, height = 6, dpi = 600)
cat(sprintf("已生成: %s\n", save_path_2))


p3 <- shap.plot.dependence(data_long = shap_long, x = "UI", y = "UI", color_feature = "auto") +
  theme_bw() +
  ggtitle("SHAP Dependence Plot: Urbanization Intensity (UI)") +
  theme(axis.text = element_text(color = "black")) 

save_path_3 <- file.path(output_dir, "3_UI_Dependence.png")
ggsave(save_path_3, p3, width = 7, height = 5, dpi = 600)
cat(sprintf("已生成: %s\n", save_path_3))


cat("正在生成所有因子以 Urbanization 为横轴的 SHAP 依赖图...\n")


shap_vals_df <- as.data.frame(shap_score)
shap_vals_df$Urbanization_Level <- urban_gradient
shap_long_all <- shap_vals_df %>%
  pivot_longer(
    cols = -Urbanization_Level, # 排除参照列进行宽转长
    names_to = "Feature",
    values_to = "SHAP_Value"
  )


target_order <- c("PRE", "NDVI", "UI", "GDP", "PD", "DEM", "TEM", "Slope", "LF")
shap_long_all$Feature <- factor(shap_long_all$Feature, levels = target_order)


p4 <- ggplot(shap_long_all, aes(x = Urbanization_Level, y = SHAP_Value)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.7) +
  
  geom_point(alpha = 0.4, color = "#4682B4", size = 1.2) +
  
  geom_smooth(method = "loess", color = "#e74c3c", se = FALSE, linewidth = 1.2) +
  
  facet_wrap(~ Feature, scales = "free_y", ncol = 3) +
  
  labs(x = "Urbanization Level",
       y = "SHAP Value (Directional Impact on ES)") +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.background = element_rect(fill = "#e9ecef", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    panel.grid.minor = element_blank()
  )

save_path_4 <- file.path(output_dir, "4_All_Features_Urbanization_Dependence_Grid.png")
ggsave(save_path_4, p4, width = 12, height = 9, dpi = 600)
cat(sprintf("已生成: %s\n", save_path_4))

cat("\n所有分析及绘图全部完成！\n")