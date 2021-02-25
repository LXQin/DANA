#' Convert gene-wise clusters into nested list
listClusters <- function(geneClusters) {
  geneClusters <- as.factor(geneClusters)
  clusters <- list()
  clusters <- lapply(
    levels(geneClusters),
    function(x) names(geneClusters)[geneClusters==x]
  )
  names(clusters) <- levels(geneClusters)
  return(clusters)
}
