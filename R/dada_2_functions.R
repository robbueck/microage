# dada2 additional functions

library(ape)
library(btools)
library(rBLAST)
# loss_ratio <- function(reduced, total) {
#   if(total == 0 & reduced == 0){
#     return(1)
#   } else if(reduced == 0){
#     warning("Total is empty")
#   } else {
#     return(  (total - reduced)/total)
#   }
# }
get_n_cores <- function(nc){
  if(is.logical(nc)){
    if(nc) {
      possibleError <- try(
        expr = {
          nc <- parallelly::availableCores(which = "all", methods = "SGE", omit = 1) + 1
        }
      )
      if(inherits(possibleError, "try-error")){
        warning("No multithreading available, this task will run on one core. If needed, specify number of cores explicitly")
        nc <- 1
      }
    } else {
      nc <- 1
    }
  } 
  return(nc)
}

loss_ratio <- function(reduced, total) {
  res <- (total - reduced)/total
  res[is.na(res)] <- 1
  return(res)
}

# plot average quality scores for each cluster
plot_cluster_qual <- function(cluster_matrix){
  quality_df <- data.frame(position = 1:ncol(cluster_matrix),
                           mean_qual = colMeans(cluster_matrix),
                           sd = matrixStats::colSds(cluster_matrix))
  ggplot(quality_df, aes(x=position, y=mean_qual)) + 
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=mean_qual-sd, ymax=mean_qual+sd), width=.2,
                  position=position_dodge(0.05))
}

# perform the first filtering step and write the files, takes in a list of file prefixes
# paired end
# dafault number of cores = available cores for SGE-job
filtering_fastq <- function(sample_names, plotting = F, libname, trimLeft, truncLen,
                        maxN, maxEE, truncQ, cores = T, verbose = F) {
  # cores = get_n_cores(cores)
  # get filenames
  Ffs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(sample_names, "_dec_1.fastq.gz")))
  Rfs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(sample_names, "_dec_2.fastq.gz")))
  if(!all(file.exists(Ffs, Rfs))){
    print("missing files:")
    print(Ffs[!file.exists(Ffs)])
    print(Rfs[!file.exists(Rfs)])
    stop("Unfiltered files not complete")
  }
  if(plotting){
    plotQualityProfile(Ffs[c(1, 2)])
    plotQualityProfile(Rfs[c(1, 2)])
  }
  # filter reads
  fltFfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt_1.fastq.gz")))
  fltRfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt_2.fastq.gz")))
  names(fltFfs) <- sort(sample_names)
  names(fltRfs) <- sort(sample_names)
  lib_out <- filterAndTrim(Ffs, fltFfs, Rfs, fltRfs,
                           trimLeft = trimLeft, 
                           truncLen=truncLen,
                           maxN=maxN,
                           maxEE=maxEE,
                           truncQ=truncQ,
                           rm.phix=TRUE,
                           compress=TRUE, multithread=cores,
                           verbose = verbose)
  if(!all(file.exists(fltFfs, fltRfs))){
    print("Missing files:")
    fltFfs[!file.exists(fltFfs)]
    fltRfs[!file.exists(fltRfs)]
    warning("Filtered files not complete")
  }
  if(plotting){
    plotQualityProfile(fltFfs[c(1, 2)])
    plotQualityProfile(fltRfs[c(1, 2)])
  }
  return(lib_out)
}



# perform the first filtering step and write the files, takes in a list of file prefixes
# sinle end
# default cores = availabe cores of SGE job
filtering_fastq_single <- function(sample_names, plotting = F, libname, trimLeft, truncLen,
                            maxN, maxEE, truncQ, cores = T) {
  cores = get_n_cores(cores)
  # get filenames
  fs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(sample_names, "_dec.fastq.gz")))
  if(!all(file.exists(fs, fs))){
    print("missing files:")
    print(fs[!file.exists(fs)])
    stop("Unfiltered files not complete")
  }
  # filter reads
  fltfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt.fastq.gz")))
  names(fltfs) <- sort(sample_names)
  lib_out <- filterAndTrim(fs, fltfs,
                           trimLeft = trimLeft, 
                           truncLen=truncLen,
                           maxN=maxN,
                           maxEE=maxEE,
                           truncQ=truncQ,
                           rm.phix=TRUE,
                           compress=TRUE, multithread=cores)
  if(!all(file.exists(fltfs))){
    print("Missing files:")
    fltfs[!file.exists(fltfs)]
    warning("Filtered files not complete")
  }
  return(lib_out)
}



# learn error rates and infer sequence variants
# paired end
# default cores = available cores for SGE job
error_and_asv <- function(sample_names, libname = "all", plotting = T, debug_mode = F,
                          pool = "pseudo", nbases = 1e8, remove_empty = F, cores = T){
  cores = get_n_cores(cores)
  # learn error rates
  fltFfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt_1.fastq.gz")))
  fltRfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt_2.fastq.gz")))
  names(fltFfs) <- sort(sample_names)
  names(fltRfs) <- sort(sample_names)
  if(!remove_empty) {
    if(!all(file.exists(fltFfs, fltRfs))){
      print("Missing files:")
      print(fltFfs[!file.exists(fltFfs)])
      print(fltRfs[!file.exists(fltRfs)])
      stop("Filtered files not complete")
    }
  } else {
    fltFfs <- fltFfs[file.exists(fltFfs)]
    fltRfs <- fltRfs[file.exists(fltRfs)]
  }
  set.seed(100)
  print("Lern Error-rates")
  if(length(nbases) == 1) {
    nbases <- c(nbases, nbases)
  }
  errF <- learnErrors(fltFfs, nbases=nbases[1], multithread=cores)
  errR <- learnErrors(fltRfs, nbases=nbases[2], multithread=cores)
  if(plotting){
    plotErrors(errF, nominalQ=TRUE)
    ggsave(paste0(maindir, "/dada2/", libname, "_errors.pdf"), device = "pdf")
  }
  # sample inference
  print("sample inference")
  dadaFs <- dada(fltFfs, err=errF, multithread=cores, pool = pool)
  dadaRs <- dada(fltRfs, err=errR, multithread=cores, pool = pool)
  dadaFs[[1]]
  if(plotting){
    # quality along the sequence
    plot_cluster_qual(dadaRs[[1]]$quality)
    ggsave(paste0(maindir, "/dada2/lib_", libname, "_seq_qual.pdf"), device = "pdf")
  }
  # merge reads
  mergers <- mergePairs(dadaFs, fltFfs, dadaRs, fltRfs, verbose=TRUE)
  if(debug_mode) {
    getN <- function(x) sum(getUniques(x))
    print(sapply(mergers, getN) /sapply(dadaFs, getN))
  }
  # construct sequence table
  lib_seqtab <- makeSequenceTable(mergers)
  return(lib_seqtab)
}


# learn error rates and infer sequence variants
# single end
# default cores = available cores for SGE job
error_and_asv_single <- function(sample_names, libname = "all", plotting = T, 
                                 pool = "pseudo", nbases = 1e8, remove_empty = F,
                                 cores = T){
  cores = get_n_cores(cores)
  # learn error rates
  fltfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt.fastq.gz")))
  if(!remove_empty){
    if(!all(file.exists(fltfs))){
      print("Missing files:")
      print(fltfs[!file.exists(fltfs)])
      stop("Filtered files not complete")
    }
  } else {
    fltfs <- fltfs[file.exists(fltfs)]
  }

  names(fltfs) <- sort(sample_names)
  set.seed(100)
  print("Lern Error-rates")
  err <- learnErrors(fltfs, nbases = nbases, multithread=cores)
  if(plotting){
    plotErrors(err, nominalQ=TRUE)
    ggsave(paste0(maindir, "/dada2/", libname, "_errors.pdf"), device = "pdf")
  }
  # sample inference
  print("sample inference")
  dadas <- dada(fltfs, err=err, multithread=cores, pool = pool)
  if(plotting){
    # quality along the sequence
    plot_cluster_qual(dadas[[1]]$quality)
    ggsave(paste0(maindir, "/dada2/lib_", libname, "_seq_qual.pdf"), device = "pdf")
  }
  # construct sequence table
  lib_seqtab <- makeSequenceTable(dadas)
  return(lib_seqtab)
}



# build tree
build_tree <- function(x, cores = NULL) {
  sequences <- getSequences(x)
  names(sequences)<-sequences
  alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA, processors = cores)
  print("finished alignemnt")
  phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
  dm <- phangorn::dist.ml(phang.align)
  treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
  print("finished NJ")
  return(treeNJ)
}


# build tree for phyloseqobject
build_tree_ps <- function(x, cores = NULL) {
  alignment <- AlignSeqs(x@refseq, anchor=NA, processors = cores)
  print("finished alignemnt")
  phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
  dm <- phangorn::dist.ml(phang.align)
  treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
  print("finished NJ")
  return(treeNJ)
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


# load a blastdb for a file:
load_blastdbs <- function(s) {
  message("searching for a blastdb for: ", s)
  possibleError <- try(
    expr = {
      blast_db <- blast(db=s, type = "blastn")
      message("Found a blastdb: ", s)
    }
  )
  if(inherits(possibleError, "try-error")){
    message("No blastdb found, constructing one for: ", s)
    makeblastdb(file = s, dbtype = "nucl")
    blast_db <- blast(db=s,  type = "blastn")
  }
  return(blast_db)
}


# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxa_to_matrix <- function(ids){
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
  taxa <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    tx <- x$taxon[m]
    lowest_rank <- tail(na.omit(tx), 1)
    tx[is.na(tx) & seq_along(tx) > max(which(!is.na(tx)))] <- lowest_rank
    # taxa[startsWith(taxa, "unclassified_")] <- NA
    tx
  }))
  colnames(taxa) <- ranks
  rownames(taxa) <- getSequences(st.all.nochim)
  return(taxa)
}

# for a set of sequence seq_set returns the longest common subsequence
get_consensus_seq <- function(seq_set, cores = NULL) {
    if (length(seq_set) < 2) {
    warning("Encountered a sequence set with less than two sequences")
  } else if(length(seq_set) == 2){ # speed up things
    alignment <- AlignProfiles(seq_set[1], seq_set[2], restrict = c(-500, 2, 10))
  } else {
    alignment <- AlignSeqs(seq_set, anchor=NA, processors = cores, verbose = F, iterations = 1)
  }
  m <- consensusMatrix(alignment, baseOnly = T)
  # get all positions where all seqs are equal
  good_cols <- colMaxs(m) == length(alignment)
  
  # get start and end of longest seq of equal columns
  tt <- rle(good_cols)
  if(length(tt$lengths) > 1) {
    t1 <- which.max(tt$lengths[tt$values])
    tt$values[tt$values][-t1] <- FALSE
    if(which.max(tt$values) == 1){ # sequence is at the beginning
      s_start <- 1
    } else {
      s_start <- sum(tt$lengths[1:(which.max(tt$values)-1)]) + 1
    }
    s_end <- s_start - 1 + tt$lengths[which.max(tt$values)]
    return(subseq(alignment[1], start = s_start, end = s_end))
  } else {
    return(alignment[1])
  }
}


get_asv_cluster <- function(blast_out, plot = F, whole_graph = F) {
  # pick the best (longest) of mutiple alignments for a query-subject pair
  blast_out <- blast_out %>% 
    arrange(desc(length)) %>%
    filter( !duplicated(paste0(qseqid, sseqid)) ) %>%
    arrange(qseqid)
    # mutate(sseqid = paste0("subject_", sseqid))
  two_mode_g <- graph.data.frame(blast_out[,c("qseqid","sseqid")], directed = F)
  if(plot){
    V(two_mode_g)$type <- bipartite.mapping(two_mode_g)$type
    V(two_mode_g)$color <- ifelse(V(two_mode_g)$type, "lightblue", "salmon")
    V(two_mode_g)$shape <- ifelse(V(two_mode_g)$type, "circle", "square")
    E(two_mode_g)$color <- "black"
    pdf(paste(str_trunc(paste(merged_set, collapse = "_"), width = 150, side = "center", ellipsis = "_"), "igraph.pdf", sep = "_"))
    plot(two_mode_g, vertex.label = NA, vertex.size=1, layout = layout_with_graphopt)
    dev.off()
  }
  if(whole_graph) {
    return(two_mode_g)
  } else {
    cs <- components(two_mode_g, mode = "weak") %>% communities()
    return(cs)
  }
}


# merge phyloseq objects and merge ASVs according to the list of clusters given in clusters
# if consensus_seq = T, for each cluster the consensus sequence is calculated. For speedup, turn this off
# then a random sequence is chosen to represent the new cluster. not suited for phylogenetic diversity calculation afterwards
# simple_merge ignores ASV information during merging and merges on taxonomic information only
merge_multiple_clusters <- function(ps_list, clusters = NULL, consensus_seq = T, cores = 1, simple_merge = F) {
  plan(multisession, workers = min(4, cores))
  if(consensus_seq & simple_merge){warning("consensus sequence will not be generated during simple_merge")}
  if(!is.null(clusters) & simple_merge){warning("cluster information will be ignored during simple_merge")}
  if(!is.null(clusters)){
    # check if ASV-names match in phyloseq objects and clusters
    if(!all(unlist(clusters) %in% unlist(lapply(ps_list, taxa_names)))){
      stop("ASV names in clusters don't match ASV names in phyloseq objects")
    }
  }
  print("merge OTU-tables")
  merged_otu_table <- purrr::reduce(lapply(ps_list, otu_table), function(x, y) {merge_phyloseq(x, y)})
  rowsums_b  <- rowSums(merged_otu_table)
  if(!simple_merge){
    merged_otu_table %<-% {purrr::reduce(clusters, function(x, y) {merge_taxa(x, y, 1)},
                                         .init = merged_otu_table)}
  }
  gc()
  print("merge sample data")
  merge_sample_data %<-% {purrr::reduce(lapply(ps_list, sample_data),
                                     function(x, y) {merge_phyloseq(x, y)})}
  print(sort(unique(merge_sample_data$study)))
  
  print("merge tax tables")
  merged_tax_table <- purrr::reduce(lapply(ps_list, tax_table),
                                    function(x, y) {merge_phyloseq(x, y)})
  if(!simple_merge){
    merged_tax_table %<-% {purrr::reduce(clusters, function(x, y) {merge_taxa(x, y, 1)},
                                         .init = merged_tax_table)}
  }
  if(!simple_merge){
    # merge taxa within merged objects
    print("merge refseq objects")
    merged_refseq <- purrr::reduce(lapply(ps_list, refseq), function(x, y) {merge_phyloseq(x, y)})
    nclusters <- length(clusters)
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nclusters,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)
    merge_consensus_refseqs <- function(x, y) {
      pb$tick()
      x[y] <- get_consensus_seq(x[y], cores = cores)
      x <- merge_taxa(x, y, 1)
      return(x)
    }
    if(consensus_seq){
      print("obtain consensus sequences")
      merged_refseq <- {purrr::reduce(clusters, merge_consensus_refseqs,
                                      .init = merged_refseq)}
    } else {
      merged_refseq %<-% {purrr::reduce(clusters,
                                        function(x, y) {pb$tick(); merge_taxa(x, y, 1)},
                                        .init = merged_refseq)}
    }
    combined_ps <- merge_phyloseq(merged_otu_table, merged_tax_table, merged_refseq, merge_sample_data)
  } else {
    combined_ps <- merge_phyloseq(merged_otu_table, merged_tax_table, merge_sample_data)
    
  }
  print("last step")
  if(!all(rowSums(merged_otu_table) == rowsums_b)){
    stop("Error in merge_taxa(merged_otu_table, cluster_1, 1), rowsums are not equal before and after merging cluster", c)
  }
  print("finished merging clusters")
  if(!simple_merge) {
    taxa_names(combined_ps) <- paste0("ASV", 1:ntaxa(combined_ps))
  }
  return(combined_ps)
}



# performs blast search and filtering of results for two input ps-objects
blast_ps_function <- function(query_name, subject_name, ps_list, cores = 2){
  print(paste0("Perform blast search ", query_name, " vs ", subject_name))
  query_ps <- ps_list[[query_name]]
  subject_db <- load_blastdbs(list.files(paste("../study_data", subject_name, sep = "/"),
                                                     pattern = "^ASVs.fasta$",
                                                     full.names = T))
  blast_cmd <- paste0("-max_hsps 10 -num_threads ", cores, " -word_size 120 -perc_identity 99")
  print(paste0("BLAST options: ", blast_cmd))
  query_ps@refseq <- RemoveGaps(query_ps@refseq, removeGaps = "all", processors = cores)
  start_time <- Sys.time()
  blast_results_small <- predict(object = subject_db,
                                 newdata = query_ps@refseq,
                                 BLAST_args = blast_cmd,
                                 custom_format = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen slen qframe sframe") %>%
    mutate(minlen = pmin(qlen, slen),
           sseqid = paste0(subject_name, "_", sseqid),
           qseqid = paste0(query_name, "_", qseqid)) %>%  # only use alignments in which the smaller sequence is covered at least 50 %
    filter(mismatch == 0, gapopen == 0, qframe == 1, sframe == 1, length >= 0.5 * minlen, length >= 200)#,
           # ((qstart <= 10) | (qend >= qlen - 10)),
           # ((sstart <= 10) | (send >= slen - 10)))
  end_time <- Sys.time()
  print("finished BLAST in:")
  print(end_time - start_time)
  return(blast_results_small)
}



# Plot the ratio of unassigned reads at different taxonomic levels
# adds a horizontal line at the median ratio of annotated reads over all ranks
annotation_ratio_plot <- function(ps,
                                  levels = c("phylum", "class", "order", "family", "genus"),
                                  x = "rank",
                                  color = "rank"){
  get_unassigned_ratio <- function(x) {
    y <- ps %>%
      microbiome::aggregate_taxa(level = x) %>%
      otu_table() %>%
      data.frame
    return(1-(y[c("Unknown"),] / colSums(y)))
  }
  ratios <- map(levels, get_unassigned_ratio)
  ratios_df <- suppressMessages(purrr::reduce(ratios, full_join))
  rownames(ratios_df) <- levels
  ratios_df_long <- ratios_df %>% tibble::rownames_to_column(var = "rank") %>%
    pivot_longer(-rank, names_to = "Sample", values_to = "values")
  ratios_df_long$rank <- factor(ratios_df_long$rank, levels = levels)
  ratios_df_long <- left_join(ratios_df_long, data.frame(sample_data(ps)), by = c("Sample" = "run_accession"))
  median_line <- median(ratios_df_long$values)
  pl <- ggplot(ratios_df_long, aes(x = !!sym(x), y = values, color = !!sym(color))) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom() +
    ylab("Annotated reads") +
    geom_hline(yintercept=median_line)
    theme(legend.position="none")
  return(pl)
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
