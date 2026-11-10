# Run after installing the package: Rscript tests/cv-regression.R
library(CALF)

set.seed(42)
binary <- data.frame(
  target = rep(c(0, 1), each = 24),
  marker1 = c(rnorm(24), rnorm(24, mean = 2)),
  marker2 = rnorm(48),
  marker3 = rnorm(48)
)

for (fn in list(cv.calf, perm_target_cv.calf)) {
  for (criterion in c("pval", "auc")) {
    output_dir <- tempfile("calf-cv-")
    dir.create(output_dir)
    result <- fn(binary, limit = 2, proportion = 0.75, times = 2,
                 targetVector = "binary", optimize = criterion,
                 outputPath = output_dir)
    stopifnot(is.matrix(result), nrow(result) == 2,
              ncol(result) == ncol(binary) + if (criterion == "pval") 2 else 1,
              all(is.finite(result)),
              file.exists(file.path(output_dir, "binary.csv")),
              file.exists(file.path(output_dir, "binaryKept.csv")),
              file.exists(file.path(output_dir, "binaryUnkept.csv")),
              length(readLines(file.path(output_dir, "AUCUnkeptList.txt"))) == 2)
  }
}

nonbinary <- data.frame(
  target = seq_len(48) / 48 + rnorm(48, sd = 0.01),
  marker1 = seq_len(48) / 48,
  marker2 = rnorm(48)
)
for (fn in list(cv.calf, perm_target_cv.calf)) {
  result <- fn(nonbinary, limit = 1, proportion = 0.75, times = 2,
               targetVector = "nonbinary", optimize = "corr")
  stopifnot(is.matrix(result), nrow(result) == 2,
            ncol(result) == ncol(nonbinary) + 1)
}

stopifnot(inherits(try(cv.calf(binary, limit = 2, proportion = 1,
                              times = 1, targetVector = "binary"), silent = TRUE),
                   "try-error"))

# Independent numerical checks: AUC and ROC, including tied scores.
stopifnot(identical(CALF:::compute.auc(c(4, 5, 6), c(1, 2, 3)), 1),
          identical(CALF:::compute.auc(c(1, 2, 3), c(4, 5, 6)), 1),
          identical(CALF:::compute.auc(rep(1, 3), rep(1, 3)), 0.5))

one_marker <- data.frame(
  target = rep(c(0, 1), each = 4),
  marker = c(1, 2, 3, 4, 5, 6, 7, 8)
)
one_fit <- calf(one_marker, nMarkers = 1, targetVector = "binary", optimize = "auc")
stopifnot(one_fit$auc == 1, nrow(one_fit$selection) == 1L,
          one_fit$selection$Marker == "marker")
roc_coordinates <- one_fit$rocPlot$data
stopifnot(identical(unname(unlist(roc_coordinates[1, c("x", "y")])), c(0, 0)),
          identical(unname(unlist(tail(roc_coordinates[, c("x", "y")], 1))), c(1, 1)))

# The negative weight should be selected for an inversely related target.
set.seed(123)
negative_marker <- data.frame(
  target = seq_len(20),
  marker = -seq_len(20)
)
negative_fit <- calf(negative_marker, 1, "nonbinary")
stopifnot(negative_fit$selection$Weight == -1,
          isTRUE(all.equal(negative_fit$finalBest, 1)))

# Undefined scores must not produce invalid marker selections.
constant_data <- data.frame(target = rep(c(0, 1), each = 4), marker = rep(4, 8))
stopifnot(inherits(try(calf(constant_data, 1, "binary", "pval"), silent = TRUE),
                   "try-error"),
          inherits(try(calf(transform(one_marker, marker = NA_real_), 1,
                                "binary", "auc"), silent = TRUE), "try-error"),
          inherits(try(calf_exact_binary_subset(one_marker, 1, nCase = 5,
                                               nControl = 3), silent = TRUE), "try-error"))

# Output helpers must restore the caller's warning preference.
set.seed(99)
subset_fit <- calf_subset(binary, 1, targetVector = "binary", times = 2)
original_warn <- getOption("warn")
options(warn = 2)
write.calf_subset(subset_fit, tempfile(fileext = ".csv"))
stopifnot(isTRUE(getOption("warn") == 2))
options(warn = original_warn)

message("CALF CV regression tests passed")
