#'@method print calf
#'@export

print.calf <- function(x, ...) {
  print.data.frame(data.frame(Variable = na.omit(x$predName),
                              Weight = na.omit(x$weight)),
                   row.names = FALSE, 
                   check.names = FALSE)
  finScore <- tail(na.omit(x$trainVec), 1)
  cat("\n Final score (", x$score, "): ", finScore, "\n", sep = "")
}