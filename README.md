# MK_Sen
Beyond the dichotomy: Non-linear water-soil-carbon coordination across the mountain urban-rural continuum

library(terra)
library(parallel)


# 1. 路径
input_dir  <- "F:/karst/ess_output/penalty"
output_dir <- "F:/karst/shiyan_urban_rural/Sen_slope"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
temp_dir <- file.path(output_dir, "temp")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)


terraOptions(tempdir = temp_dir, progress = 1, memfrac = 0.7)


# 2. 参数
years <- 2000:2020
alpha <- 0.05
resample_method <- "bilinear"
ncores <- max(1L, detectCores(logical = FALSE) - 1L)
wopt_float <- list(gdal = c("COMPRESS=LZW"), datatype = "FLT4S")
wopt_int   <- list(gdal = c("COMPRESS=LZW"), datatype = "INT2S")


# 3. 输入文件
input_files <- file.path(
  input_dir,
  sprintf("CCD_D_effective_%d.tif", years)
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "以下输入栅格不存在：\n",
      paste(missing_files, collapse = "\n")
    )
  )
}


# 4. 2000年模板
template_file <- file.path(input_dir, "CCD_D_effective_2000.tif")
template <- rast(template_file)

if (is.na(crs(template)) || crs(template) == "") {
  warning("2000 年模板缺少坐标系，请手工给 template 指定 WGS_1984_Albers。")
}

same_grid <- function(r1, r2, tol = 1e-10) {
  same_dim <- nrow(r1) == nrow(r2) && ncol(r1) == ncol(r2)
  same_res <- isTRUE(all.equal(res(r1), res(r2), tolerance = tol))
  same_ext <- isTRUE(all.equal(as.vector(ext(r1)), as.vector(ext(r2)), tolerance = tol))
  same_dim && same_res && same_ext
}


# 5. 对齐并缓存
aligned_dir <- file.path(output_dir, "aligned")
dir.create(aligned_dir, recursive = TRUE, showWarnings = FALSE)

aligned_files <- character(length(years))

for (i in seq_along(years)) {
  yr <- years[i]
  in_file  <- input_files[i]
  out_file <- file.path(aligned_dir, sprintf("aligned_%d.tif", yr))
  aligned_files[i] <- out_file
  
  
  if (file.exists(out_file)) {
    cat("年份", yr, "的对齐文件已存在，跳过。\n")
    next
  }
  
  cat("正在对齐年份：", yr, "\n")
  r <- rast(in_file)
  
  if (is.na(crs(r)) || crs(r) == "") {
    warning(paste0("年份 ", yr, " 缺少坐标系，已赋值为模板坐标系。"))
    crs(r) <- crs(template)
  }
  
  if (!same.crs(r, template)) {
    cat("  -> CRS 不一致，执行 project()\n")
    r_aligned <- project(
      r, template,
      method   = resample_method,
      threads  = TRUE,
      use_gdal = TRUE,
      by_util  = TRUE
    )
  } else if (!same_grid(r, template)) {
    cat("  -> 网格未对齐，执行 resample()\n")
    r_aligned <- resample(
      r, template,
      method   = resample_method,
      threads  = TRUE,
      by_util  = TRUE
    )
  } else {
    cat("  -> 已对齐，直接使用\n")
    r_aligned <- r
  }
  
  names(r_aligned) <- paste0("Y", yr)
  writeRaster(r_aligned, out_file, overwrite = TRUE, wopt = wopt_float)
  
  rm(r, r_aligned)
  gc()
}


# 6. 加载时间序列栅格
s <- rast(aligned_files)
names(s) <- paste0("Y", years)


# 7. 多年平均值：直接用 terra::mean（更快）
cat("开始计算多年平均值...\n")

mean_file <- file.path(output_dir, "eci_mean_2016_2020.tif")
if (!file.exists(mean_file)) {
  mean_raster <- mean(
    s,
    na.rm = TRUE,
    filename = mean_file,
    overwrite = TRUE,
    wopt = wopt_float
  )
} else {
  mean_raster <- rast(mean_file)
  cat("多年平均值文件已存在，跳过重算。\n")
}


# 8. 预计算配对索引（关键提速点）
pair_idx <- combn(seq_along(years), 2)
p1 <- pair_idx[1, ]
p2 <- pair_idx[2, ]
dt_full <- years[p2] - years[p1]


# 9.MK + Sen 函数
mk_sen_fast <- function(v, p1, p2, dt_full, alpha) {
  ok <- !is.na(v)
  n <- sum(ok)
  
  if (n < 2) {
    return(c(NA, NA, NA, NA, NA))
  }
  
  ok_pairs <- ok[p1] & ok[p2]
  if (!any(ok_pairs)) {
    return(c(NA, NA, NA, NA, NA))
  }
  
  dy <- v[p2[ok_pairs]] - v[p1[ok_pairs]]
  slopes <- dy / dt_full[ok_pairs]
  
  sen_slope <- median(slopes)
  
  S <- sum(sign(dy))
  
  y <- v[ok]
  ties_tab <- table(y)
  ties_tab <- ties_tab[ties_tab > 1]
  
  tie_term <- 0
  if (length(ties_tab) > 0) {
    tie_term <- sum(ties_tab * (ties_tab - 1) * (2 * ties_tab + 5))
  }
  
  varS <- (n * (n - 1) * (2 * n + 5) - tie_term) / 18
  
  if (varS == 0) {
    Z <- 0
    p <- 1
  } else {
    if (S > 0) {
      Z <- (S - 1) / sqrt(varS)
    } else if (S < 0) {
      Z <- (S + 1) / sqrt(varS)
    } else {
      Z <- 0
    }
    p <- 2 * (1 - pnorm(abs(Z)))
  }
  
  tau <- S / (n * (n - 1) / 2)
  
  trend_class <- 0L
  if (is.na(sen_slope) || is.na(p)) {
    trend_class <- NA_integer_
  } else if (sen_slope < 0 && p < alpha) {
    trend_class <- -2L
  } else if (sen_slope < 0 && p >= alpha) {
    trend_class <- -1L
  } else if (sen_slope > 0 && p >= alpha) {
    trend_class <- 1L
  } else if (sen_slope > 0 && p < alpha) {
    trend_class <- 2L
  } else {
    trend_class <- 0L
  }
  
  c(sen_slope, Z, p, tau, trend_class)
}


# 10. 趋势计算
trend_file <- file.path(output_dir, "eci_trend_2000_2020.tif")

if (!file.exists(trend_file)) {
  cat("开始计算 MK + Sen 趋势，使用", ncores, "个核...\n")
  
  trend_stack <- app(
    s,
    fun = mk_sen_fast,
    p1 = p1,
    p2 = p2,
    dt_full = dt_full,
    alpha = alpha,
    cores = ncores,
    filename = trend_file,
    overwrite = TRUE,
    wopt = wopt_float
  )
  
  names(trend_stack) <- c("sen_slope", "mk_z", "mk_p", "mk_tau", "trend_class")
} else {
  trend_stack <- rast(trend_file)
  names(trend_stack) <- c("sen_slope", "mk_z", "mk_p", "mk_tau", "trend_class")
  cat("趋势总文件已存在，跳过重算。\n")
}


# 11. 拆分输出
out_slope <- file.path(output_dir, "eci_sen_slope_2000_2020.tif")
out_z     <- file.path(output_dir, "eci_mk_z_2000_2020.tif")
out_p     <- file.path(output_dir, "eci_mk_p_2000_2020.tif")
out_tau   <- file.path(output_dir, "eci_mk_tau_2000_2020.tif")
out_class <- file.path(output_dir, "eci_trend_class_2000_2020.tif")

if (!file.exists(out_slope)) {
  writeRaster(trend_stack[[1]], out_slope, overwrite = TRUE, wopt = wopt_float)
}
if (!file.exists(out_z)) {
  writeRaster(trend_stack[[2]], out_z, overwrite = TRUE, wopt = wopt_float)
}
if (!file.exists(out_p)) {
  writeRaster(trend_stack[[3]], out_p, overwrite = TRUE, wopt = wopt_float)
}
if (!file.exists(out_tau)) {
  writeRaster(trend_stack[[4]], out_tau, overwrite = TRUE, wopt = wopt_float)
}
if (!file.exists(out_class)) {
  writeRaster(round(trend_stack[[5]]), out_class, overwrite = TRUE, wopt = wopt_int)
}

cat("\n处理完成。\n")
cat("输出目录：", output_dir, "\n")
