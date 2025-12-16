# alternative models:

# adapt random forest model to be trained over different parameters:
adapt_ranger <- getModelInfo("ranger", regex = FALSE)[[1]]
adapt_ranger$method <- "adapted_ranger"
adapt_ranger$parameters <- data.frame(parameter = c("mtry", "splitrule", "min.node.size", "num.trees", "replace"),
                                     class = c("numeric", "character", "numeric", "numeric", "logical"),
                                     label = c("#Randomly Selected Predictors",
                                               "Splitting Rule",
                                               "Minimal Node Size",
                                               "number of trees to grow",
                                               "replacement"))
adapt_ranger$grid <- function(x, y, len = NULL, search = "grid") {
  if(search == "grid") {
    srule <-
      if (is.factor(y))
        "gini"
    else
      "variance"
    out <- expand.grid(mtry =
                         caret::var_seq(p = ncol(x),
                                        classification = is.factor(y),
                                        len = len),
                       min.node.size = ifelse( is.factor(y), 1, 5),
                       splitrule = c(srule, "extratrees")[1:min(2,len)],
                       num.trees = 50,
                       replace = TRUE)
  } else {
    srules <- if (is.factor(y))
      c("gini", "extratrees")
    else
      c("variance", "extratrees", "maxstat")
    out <-
      data.frame(
        min.node.size= sample(1:(min(20,nrow(x))), size = len, replace = TRUE),
        mtry = sample(1:ncol(x), size = len, replace = TRUE),
        splitrule = sample(srules, size = len, replace = TRUE,),
        num.trees = c(50, 100, 500)
      )
  }
}
