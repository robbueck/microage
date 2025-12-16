clr <- function(x){
  log(x)-mean(log(x))
}

# adapted clr transformation from maya

adapt_clr <- function(counts){
  pseudo <- rep(0,nrow(counts))
  n_otu <- ncol(counts)
  max_depth_index <- which(rowSums(counts)==max(rowSums(counts)))
  max_depth_sample <- counts[max_depth_index,]
  gnz_max <- sum(log(max_depth_sample[max_depth_sample!=0])) # sum of log of counts from sample with highest counts
  pseudo_max <- min(min(counts[counts!=0])/10,1) # smallest count / 10
  n_zero_max <- sum(counts[max_depth_index,]==0) # how many zeros in highest count sampe
  for(i in 1:nrow(counts)){
    temp_sample <- counts[i,]
    n_zero_temp <- sum(temp_sample==0)
    gnz_temp <- sum(log(temp_sample[temp_sample!=0]))
    pseudo[i] <- exp((1/(n_zero_temp-n_otu))*((n_zero_max-n_otu)*log(pseudo_max)+gnz_max-gnz_temp))
  }
  if(any(pseudo==0)){
    pseudo[pseudo==0]<-min(pseudo)
  }
  mzero<-(counts==0)
  counts <- counts + pseudo[row(counts)]*mzero
  counts_clr <- t(apply(counts, 1, clr))
  return(list(counts_clr, pseudo))
}



# fix merge_taxa2 naming
merge_taxa2_fixed <- function (x, taxa = NULL, pattern = NULL, name = "Merged")
{
  if (is.null(taxa) && is.null(pattern)) {
    return(x)
  }
  if (!is.null(pattern)) {
    if (!is.null(taxa)) {
      mytaxa <- taxa
    }
    else {
      mytaxa <- taxa(x)
    }
    if (length(grep(pattern, mytaxa)) == 0) {
      return(x)
    }
    mytaxa <- mytaxa[grep(pattern, mytaxa)]
  }
  else if (is.null(taxa)) {
    mytaxa <- taxa(x)
  }
  else {
    mytaxa <- taxa
  }
  x2 <- phyloseq::merge_taxa(x, mytaxa, 1)
  mytaxa <- gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mytaxa))
  taxa_names(x2) <- gsub(mytaxa[[1]], name, taxa_names(x2))
  tax_table(x2)[1, ] <- rep(name, ncol(tax_table(x2)))
  x2
}



# plot diversity indices:
# set colour to "read_counts_final" to color the points according to the number of reads used to calculate the indices
diversity_plots <- function(ps_object, x = NULL, colour = "read_counts_final", shape = NULL,
                            measures = c("observed", "chao1", "diversity_shannon", "evenness_simpson", "dominance_simpson", "PD"),
                            angle = 0, vjust = 1, hjust=1, alpha = 1, size = 1){
  require(tidyr)
  metadata <- sample_data(ps_object)
  alpha.div <- microbiome::alpha(ps_object, index = c("Observed","Chao1", "Shannon", "Simpson")) 
  alpha.div<- as.data.frame(alpha.div)
  if("PD" %in% measures) {
    set.seed(711)
    phy_tree(ps_object) <- ape::root(phy_tree(ps_object), sample(taxa_names(ps_object), 1), resolve.root = TRUE)
    pd <- btools::estimate_pd(ps_object)
    alpha.div<- cbind(alpha.div, pd)
  }
  metadata <- cbind(metadata, alpha.div)
  metadata <- cbind2(metadata, as.data.frame(sample_sums(ps_object)))
  colnames(metadata)[ncol(metadata)] <- "read_counts_final"
  metadata_long <- metadata %>% 
    pivot_longer(
      cols = all_of(measures), 
      names_to = "index",
      values_to = "value"
    ) %>%
    arrange(lifestyle)
  # mutate(lifestyle = factor(lifestyle, levels = c("non_industrialized", "industrialized")))
  richness_map <- aes_string(x = x, y = "value", colour = colour, shape = shape)
  p <- ggplot(metadata_long, richness_map) + # aes(x = age, y = value, color = age)) +
    geom_point(alpha = alpha, size = size) +
    {if(length(unique(measures)) > 1)facet_wrap(~index, scales = "free_y")} +
    {if(class(metadata[,c(colour)]) == "numeric")scale_color_gradientn(name = colour, colors = topo.colors(7))} +
    ylab("Alpha Diversity Measure") +
    theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust=hjust))
  return(p)
}


# plot a PCoA ordination
plot_pcoa <- function(ps_object, ordination,
                      color = NULL,
                      shape = NULL,
                      ellipses = F,
                      label = NULL,
                      title = NULL,
                      alpha = 1,
                      size = 1,
                      axes = 1:2,
                      plot_dens = F){
  DF <- plot_ordination(ps_object, ordination, justDF = T, axes = axes)
  if (!is.null(color)) {
    if (!color %in% names(DF)) {
      warning("Color variable was not found in the available data you provided.", 
              "No color mapped.")
      color <- NULL
    }
  }
  if (!is.null(shape)) {
    if (!shape %in% names(DF)) {
      warning("Shape variable was not found in the available data you provided.", 
              "No shape mapped.")
      shape <- NULL
    }
  }
  if (!is.null(label)) {
    if (!label %in% names(DF)) {
      warning("Label variable was not found in the available data you provided.", 
              "No label mapped.")
      label <- NULL
    }
  }
  x = colnames(DF)[1]
  y = colnames(DF)[2]
  if (ncol(DF) <= 2) {
    message("No available covariate data to map on the points for this plot `type`")
    ord_map = aes_string(x = x, y = y)
  } else {ord_map = aes_string(x = x, y = y, color = color, shape = shape, 
                               na.rm = TRUE)
  }
  p <- ggplot(DF, ord_map) + geom_point(na.rm = TRUE, alpha = alpha, size = size)
  if (!is.null(label)) {
    label_map <- aes_string(x = x, y = y, label = label)
    p = p + geom_text(label_map, data = rm.na.phyloseq(DF, 
                                                       label), size = 2, vjust = 1.5, na.rm = TRUE)
  }
  if (!is.null(title)) {
    p = p + ggtitle(title)
  }
  if (length(ordination$values$Eigenvalues[axes]) > 0) {
    eigvec = ordination$values$Eigenvalues
    fracvar = eigvec[axes]/sum(eigvec)
    percvar = round(100 * fracvar, 1)
    strivar = as(c(p$label$x, p$label$y), "character")
    strivar = paste0(strivar, "   [", percvar, "%]")
    p = p + xlab(strivar[1]) + ylab(strivar[2])
  }
  if (ellipses) {
    p <- p + stat_ellipse(ord_map)
  }
  if (plot_dens) {
    xdens <- cowplot::axis_canvas(p, axis = "x")+
      geom_density(data = DF, aes(x = Axis.1, fill = !!sym(color)),
                   alpha = 0.7, size = 0.2)
    
    ydens <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE)+
      geom_density(data = DF, aes(x = Axis.2, fill = !!sym(color)),
                   alpha = 0.7, size = 0.2)+
      coord_flip()
    
    p <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p <- insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")
    p <- ggdraw(p)
  }
  return(p)
}

