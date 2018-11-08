#检验x是否为空。返回警告信息
mean_ci <- function(x, level = 0.95) {
  if (length(x) == 0) {
    warning("`x` was empty", call. = FALSE)
    return(c(-Inf, Inf))
  }
  se <- sd(x) / sqrt(length(x))
  alpha <- 1 - level
  mean(x) + se * qnorm(c(alpha / 2, 1 - alpha / 2))
}

#将向量x中的NA替换成replacement
#is.na重复两次。因此用变量ismiss替换
replace_missings <- function(x, replacement) {
  # Define is_miss
  is_miss = is.na(x)
  
  # Rewrite rest of function to refer to is_miss
  x[is_miss] <- replacement
  cat(sum(is_miss), replacement, "\n")
  x
}