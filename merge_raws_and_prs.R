
#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript merge_raws_and_prs.R <effects.tsv> <out_prs.tsv> <out_merged.tsv|-> <raw_or_dir ...>\n",
      "  effects.tsv : 2 cols (no header) Variant\\tEffect\n",
      "  out_prs.tsv : output PRS (FID IID PRS N_nonNA)\n",
      "  out_merged  : merged genotype matrix (or '-' to skip)\n",
      "  raw_or_dir  : one or more .raw/.raw.gz files or directories\n", sep="")
  quit(status=1)
}
eff_file <- args[1]
out_prs  <- args[2]
out_merged <- args[3]
inputs   <- args[-(1:3)]

# 收集所有 .raw 文件
collect_raws <- function(paths) {
  out <- character(0)
  for (p in paths) {
    if (!file.exists(p)) next
    info <- file.info(p)
    if (isTRUE(info$isdir)) {
      fs <- list.files(p, pattern="\\.raw(\\.gz)?$", full.names = TRUE, recursive = FALSE)
      out <- c(out, fs)
    } else {
      out <- c(out, p)
    }
  }
  unique(out)
}
raw_files <- collect_raws(inputs)
if (length(raw_files) == 0) stop("No .raw files found")

message(sprintf("[INFO] Found %d raw files", length(raw_files)))

# 读取效应表
eff <- fread(eff_file, header = FALSE, sep = "\t",
             col.names = c("variant", "effect"), na.strings=c("NA","", "."))
if (nrow(eff) == 0) stop("effects.tsv is empty")


# 合并逻辑：按 (FID, IID) 做全外连接，列并集；同名位点列冲突时：
# 取先前值与新值：若先前为 NA 则用新值；若新值为 NA 保留旧值；若均非 NA 且不同，记一次冲突并保留旧值。


merge_two <- function(A, B) {
  stopifnot(all(c("FID","IID") %in% names(A)), all(c("FID","IID") %in% names(B)))
  setkey(A, FID, IID); setkey(B, FID, IID)
  M <- merge(A, B, by=c("FID","IID"), all=TRUE, suffixes=c("", ".y"), sort=FALSE)
  ycols <- grep("\\.y$", names(M), value=TRUE)
  conflict_cnt <- 0L
  if (length(ycols) > 0) {
    base <- sub("\\.y$", "", ycols)
    for (i in seq_along(ycols)) {
      x <- base[i]; y <- ycols[i]
      # 若原列不存在（理论上应存在），直接重命名新列
      if (!x %in% names(M)) {
        setnames(M, y, x)
        next
      }
      # 合并：优先已有值；已有为 NA 用新值；两者皆非 NA 且不同 -> 冲突（保留旧值）
      both_nonNA <- !is.na(M[[x]]) & !is.na(M[[y]]) & (M[[x]] != M[[y]])
      conflict_cnt <- conflict_cnt + sum(both_nonNA, na.rm=TRUE)
      idx <- is.na(M[[x]]) & !is.na(M[[y]])
      if (any(idx)) M[[x]][idx] <- M[[y]][idx]
      M[[y]] <- NULL
    }
  }
  attr(M, "conflicts") <- as.integer(conflict_cnt)
  M
}

# 逐文件合并
read_raw <- function(f) {
  dt <- tryCatch(
    read.table(f, sep="\t", check.names = FALSE,h=F),
    error = function(e) { message("[skip read] ", f, " : ", e$message); NULL }
  )
  colnames(dt)<-dt[1,]
  dt<-dt[-1,]
  if (is.null(dt)) return(NULL)
  if (!all(c("FID","IID") %in% names(dt))) {
    message("[skip format] missing FID/IID: ", f); return(NULL)
  }
  dt
}

merged <- NULL
total_conflicts <- 0L
for (i in seq_along(raw_files)) {
  f <- raw_files[i]
  dt <- read_raw(f)
  if (is.null(dt)) next
  if (is.null(merged)) {
    merged <- dt
  } else {
    merged <- merge_two(merged, dt)
    total_conflicts <- total_conflicts + attr(merged, "conflicts")
    attr(merged, "conflicts") <- NULL
  }
  if (i %% 10 == 0 || i == length(raw_files)) {
    message(sprintf("[merge] %d/%d files processed", i, length(raw_files)))
  }
}
if (is.null(merged)) stop("No valid .raw loaded")

message(sprintf("[INFO] merged dims: %d samples x %d columns (incl. FID,IID)",
                nrow(merged), ncol(merged)))
if (total_conflicts > 0)
  message(sprintf("[WARN] genotype conflicts encountered: %d (kept earlier values)", total_conflicts))

# 计算 PRS
all_vars <- setdiff(names(merged), c("FID","IID"))
common_vars <- intersect(all_vars, eff$variant)
if (length(common_vars) == 0) stop("No overlapping variants between merged raw and effects")

# 取矩阵并对齐效应
geno <- as.matrix(merged[, ..common_vars])
mode(geno) <- "numeric"   # 强制为数值
effects_vec <- eff$effect[match(common_vars, eff$variant)]

valid_cols <- which(!is.na(effects_vec))
if (length(valid_cols) == 0) stop("All effects are NA after matching")
geno <- geno[, valid_cols, drop = FALSE]
effects_vec <- effects_vec[valid_cols]

# 分母：每样本有效位点数（非 NA）
n_nonNA <- rowSums(!is.na(geno))

# 分子：把 NA 当 0，再做矩阵乘
geno[is.na(geno)] <- 0
num <- as.numeric(geno %*% effects_vec)
PRS <- ifelse(n_nonNA > 0, num / n_nonNA, NA_real_)

# 输出 PRS
out <- data.table(FID = merged$FID, IID = merged$IID,
                  PRS = PRS, N_nonNA = as.integer(n_nonNA))
setorder(out, IID)
fwrite(out, out_prs, sep = "\t", quote = FALSE)
message(sprintf("[OK] PRS -> %s (samples=%d, variants_used=%d)",
                out_prs, nrow(out), length(valid_cols)))

# 可选输出合并矩阵
if (!identical(out_merged, "-")) {
  fwrite(merged, out_merged, sep = "\t", quote = FALSE)
  message(sprintf("[OK] merged genotype -> %s", out_merged))
}
