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


# downsample the larger lifestyle to the size of the smaller lifestyle
# optional, attempt to have equal study sizes
get_list <- function(vl, dfrm, max = NA, study_prob = F) {
  lst <- lapply(c(industrialized = "industrialized", non_industrialized = "non_industrialized"),
                function(x) {dfrm %>%
                    filter(lifestyle == x) %>%
                    pull(vl) %>% unique}) # get all ids for that lifestyle
  len <- lapply(lst, length) %>% unlist %>% min(c(., max), na.rm = T)
  print(vl)
  if(study_prob) {  # downsample the category to more equal study sizes
    probs <- dfrm %>% 
      group_by(study) %>%
      mutate(prob = 1/length(unique(!!sym(vl)))) %>%
      ungroup() %>%
      select(prob, !!sym(vl)) %>%
      distinct() %>%
      column_to_rownames(vl)
    # lapply(lst, function(x) print(probs[x,]))
    final_list <- lapply(lst, function(x) sample(x, len, prob = probs[x,]))
  } else {
    final_list <- lapply(lst, function(x) sample(x, len))
  }
  return(final_list)
}



# create a dataframe with pairwise distances in long format combined with metadata about the pairs
get_dist_df <- function(dist_obj, mdata){
  dist_df <- dist_obj %>% 
    as.matrix %>%
    as.table %>%
    as.data.frame() %>%
    mutate(sample_1 = as.character(Var1),
           sample_2 = as.character(Var2),
           distance = Freq,
           .keep = "unused") %>%
    filter(sample_1 > sample_2)
  dist_df <- left_join(dist_df, mdata[,c("run_accession", "subject_ID", "age", "age_cat", "study", "lifestyle")],
                       by = join_by(sample_1 == run_accession)) %>%
    left_join(mdata[,c("run_accession", "subject_ID", "age", "age_cat", "study", "lifestyle")],
              by = join_by(sample_2 == run_accession), suffix = c("_1", "_2")) %>%
    filter(subject_ID_1 != subject_ID_2) %>%
    mutate(same_age_cat = age_cat_1 == age_cat_2,
           same_ls = lifestyle_1 == lifestyle_2,
           same_study = study_1 == study_2)
  return(dist_df)
}


# plot a PCoA ordination for repeated_subsampling
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


# perform the repeated subsampling, either plots pcoas or
# outputs statistical results (stats_output = T)
repeated_subsampling <- function(df, dist_obj, ps_obj, max_subj = NA,
                                 max_samples = NA, max_studies = NA,
                                 stats_output = T, plots = T) {
  stds <- get_list("study", df, max = max_studies)
  print("0")
  df_red <- df %>%
    filter(lifestyle == "non_industrialized" | study %in% stds$industrialized) %>%
    filter(lifestyle == "industrialized" | study %in% stds$non_industrialized) # subsample studies
  sbjcts <- get_list("subject_ID", df, max = max_subj, study_prob = T)
  print("1")
  df_red_sb <- df_red %>%
    filter(lifestyle == "non_industrialized" | subject_ID %in% sbjcts$industrialized) %>%
    filter(lifestyle == "industrialized" | subject_ID %in% sbjcts$non_industrialized) # subsample individuals
  df_red_y <- filter(df_red_sb, age <= 365) # same distribuition for above and below one year respectively
  df_red_o <- filter(df_red_sb, age > 365)
  print("2")
  smpls_y <- get_list("sample_ID", df_red_y, max = max_samples / 2, study_prob = T)
  df_red_y_smp <- df_red_y %>%
    filter(lifestyle == "non_industrialized" | sample_ID %in% smpls_y$industrialized) %>%
    filter(lifestyle == "industrialized" | sample_ID %in% smpls_y$non_industrialized)
  print("a")
  if(length(table(df_red_o$lifestyle)) > 1) {
    print("3")
    smpls_o <- get_list("sample_ID", df_red_o, max = max_samples / 2, study_prob = T)
    df_red_o_smp <- df_red_o %>%
      filter(lifestyle == "non_industrialized" | sample_ID %in% smpls_o$industrialized) %>%
      filter(lifestyle == "industrialized" | sample_ID %in% smpls_o$non_industrialized) # subsample samples
    df_red_smp <- rbind(df_red_y_smp, df_red_o_smp)
  } else {
    print("No old samples found for one lifestyle")
    df_red_smp <- df_red_y_smp
  }
  # filter distances
  dist_obj_matrix <- dist_obj %>% as.matrix()
  dist_obj_matrix_filter <- dist_obj_matrix[df_red_smp$run_accession, df_red_smp$run_accession]
  dist_obj_filter <- as.dist(dist_obj_matrix_filter)
  print("b")
  # print("get dist density")
  if(stats_output) {
    dist_df_long <- get_dist_df(dist_obj = dist_obj_filter, mdata = df) %>%
      filter(same_ls | same_age_cat)
    # between lifestyle distances
    inter_ls_dist <- dist_df_long %>%
      filter(!same_study) %>%
      group_by(age_cat_1, same_ls) %>%
      summarise(mean_inter_ls = mean(distance), median_inter_ls = median(distance))
    # differences within lifestyles
    df_p_val_dist <- dist_df_long %>% # p-vals
      filter(same_ls) %>%
      rstatix::wilcox_test(distance ~ lifestyle_1, detailed = T)
    mean_median_df <- dist_df_long %>% # means/medians
      filter(same_ls) %>%
      group_by(lifestyle_1) %>% 
      summarise(mean = mean(distance), median = median(distance)) %>% 
      pivot_longer(cols = c(mean, median), names_to = "stat", values_to = "value") %>%
      unite("lifestyle_stat", lifestyle_1, stat, sep = "_")%>%
      pivot_wider(names_from = lifestyle_stat, values_from = value) %>%
      cbind(df_p_val_dist,.)
  }
  if(plots) {
    if(stats_output) { # if already there, just filter it
      dist_df_long <- dist_df_long %>%
        filter(same_ls)
    } else {
      dist_df_long <- get_dist_df(dist_obj = dist_obj_filter, mdata = df) %>%
        filter(same_ls)
    }
    dist_dens <- dist_df_long %>%
      ggplot(aes(x = distance, color = lifestyle_1)) +
      geom_boxplot(aes(x = distance, y = lifestyle_1, color = lifestyle_1), alpha = 0.5) +
      geom_density() +
      theme_classic() +
      theme(legend.position= c(0.2, 0.8),
            legend.title = element_text(size=10),
            axis.text.y = element_text(angle = 90),
            legend.text = element_text(size = 8))
    # print("run pcoa")
    pcoa_filter <- pcoa(dist_obj_filter)
    pcoa_ls <- plot_pcoa(ps_obj, pcoa_filter, color="lifestyle", #title="Bray NMDS",
                         alpha = 0.5, size = 0.5, ellipses = F) +
      theme_classic() +
      theme(text = element_text(size = 18),
            legend.position= "none", # c(0.2, 0.2),
            legend.title = element_text(size=10),
            legend.text = element_text(size = 8))
    pcoa_age <- plot_pcoa(ps_obj, pcoa_filter, color="age", #title="Bray NMDS",
                          alpha = 0.5, size = 0.5, ellipses = F) +
      theme_classic() +
      scale_color_gradientn(name = "age", colors = topo.colors(7)) +
      theme(text = element_text(size = 18),
            legend.position= "none", # c(0.2, 0.2),
            legend.title = element_text(size=10),
            legend.text = element_text(size = 8)) 
    age_hist <- df_red_smp %>%
      # filter(lifestyle == "non_industrialized" | study %in% sample(industrialized_subjects, 853)) %>%
      ggplot(aes(x = age, color = lifestyle)) +
      geom_histogram() + 
      facet_wrap(~lifestyle) +
      theme_classic() +
      theme(legend.position= c(0.2, 0.9),
            legend.title = element_text(size=10),
            legend.text = element_text(size = 8))
    grid.arrange(pcoa_ls, pcoa_age, age_hist, dist_dens, layout_matrix = rbind(c(1, 2), c(3, 4)))
  }
  if(stats_output) {
    return(list(within_ls = mean_median_df, between_ls_res = inter_ls_dist))
  }
}

