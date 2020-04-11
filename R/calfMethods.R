#'@method print calf
#'@export

print.calf <- function(x, ...) {
  print.data.frame(data.frame(Variable = na.omit(x$xName),
                              Weight = na.omit(x$weight)),
                   row.names = FALSE, 
                   check.names = FALSE)
  cat("\n Final score (", x$score, "): ", x$finalS, "\n", sep = "")
}

#'@method predict calf
#'@export

predict.calf <- function(object, newx) {
  bvec <- rep(0, ncol(newx))
  bvec[object$xInd] <- object$weight
  as.matrix(x) %*% bvec
}
