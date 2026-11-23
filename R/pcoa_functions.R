##Project: Microbiome age index
##Aim: General functions for the scripts
##Author: Robert Bücking


#' Count features detected in a random subset of studies
#'
#' Randomly selects a number of studies (or other categories), and counts how
#' many features (e.g. taxa) have a value above 10 in more than one sample of
#' these studies. Calling it repeatedly with increasing \code{n_cat} gives an
#' accumulation curve of detected features against the number of studies.
#'
#' @param n_cat Integer. Number of categories (studies) to sample randomly.
#' @param df A data frame with one row per sample, one numeric column per
#'   feature, and the grouping column named in \code{cat}. All columns other
#'   than \code{cat} must be numeric.
#' @param cat Character. Name of the grouping column. Default \code{"study"}.
#' @param min_count Integer. Minimal number of reads to consier a taxon detected. Default \code{10}.
#'
#' @details
#' A feature counts as detected if it has a value greater than 10 in at least
#' two samples of the selected studies. The selection is random; call
#' \code{set.seed()} beforehand for reproducible results.
#'
#' @return An integer: the number of detected features.
#'
#' @importFrom dplyr filter select %>%
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n_detected <- sapply(rep(1:10, each = 20), subsample_fn, df = count_df)
#' }
#'
#' @export
subsample_fn <- function(n_cat, df, cat = "study", min_count = 10) {
  stds <- sample(unique(df[[cat]]), n_cat)
  df %>% filter(.data[[cat]] %in% stds) %>%
    select(-cat) %>% 
    colSums(. > min_count) %>% 
    `>`(1) %>%
    sum()
}


# fix merge_taxa2 naming, obsolete with update of microbiome package
merge_taxa2_fixed <- function (x, taxa = NULL, pattern = NULL, name = "Merged") {
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


#' Downsample IDs to equal numbers per lifestyle
#'
#' Collects the unique IDs (e.g. subjects) of each lifestyle group and
#' randomly downsamples the larger group to the size of the smaller one.
#' Optionally, sampling is weighted so that studies contribute more equally.
#'
#' @param vl Character. Name of the ID column in \code{dfrm}.
#' @param dfrm A data frame with the columns named in \code{vl},
#'   \code{Lifestyle} (values \code{"industrialized"} and
#'   \code{"non_industrialized"}) and \code{study}.
#' @param study_prob Logical. If \code{TRUE}, each ID is sampled with a
#'   probability inversely proportional to the number of IDs in its study, so
#'   that large studies don't dominate the downsampled set. Default
#'   \code{FALSE}.
#' @param max Integer. Maximum number of IDs per lifestyle. No effect, if NA or 
#'  larger than each lifestyle group. Default \code{NA}.
#'
#' @details
#' The smaller group keeps all of its IDs (in random order). The weighting with
#' \code{study_prob = TRUE} makes study sizes more equal, but not exactly
#' equal. Sampling is random; call \code{set.seed()} beforehand for
#' reproducible results.
#'
#' @return A named list with the elements \code{industrialized} and
#'   \code{non_industrialized}, each a vector of the same number of IDs.
#'
#' @importFrom dplyr filter pull group_by mutate ungroup select distinct %>%
#' @importFrom tibble column_to_rownames
#' @importFrom rlang sym !!
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' ids <- get_list("subject_ID", metadata, study_prob = TRUE)
#' lengths(ids)
#' }
#'
#' @export
get_list <- function(vl, dfrm, max = NA, study_prob = F) {
  lst <- lapply(c(industrialized = "industrialized", non_industrialized = "non_industrialized"),
                function(x) {dfrm %>%
                    filter(lifestyle == x) %>%
                    pull(vl) %>% unique}) # get all ids for that lifestyle
  len <- lapply(lst, length) %>% unlist %>% min(c(., max), na.rm = T)
  if(study_prob) {  #downsample the category to more equal study sizes
    probs <- dfrm %>% group_by(study) %>%
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


#' Downsample a dataset to balance lifestyles at study, subject and sample level
#'
#' Downsamples the data in three steps so that industrialized and
#' non-industrialized lifestyles are equally represented: first studies, then
#' subjects, then samples. Samples are balanced separately for the first year
#' of life and for older ages, to keep the age distributions comparable.
#'
#' @param otu_tbl_df A data frame with one row per sample, containing the
#'   columns \code{study}, \code{subject_ID}, \code{sample_ID}, \code{age}
#'   (in days) and \code{Lifestyle}, plus one numeric column per feature.
#'
#' @details
#' Each step uses \code{\link{get_list}}. Subjects and samples are drawn with
#' study-weighted probabilities (\code{study_prob = TRUE}). Samples are
#' balanced separately for ages up to 365 days and above 365 days. If only
#' one lifestyle has samples above 365 days, only the younger samples are
#' returned. Finally, numeric columns that sum to zero (absent features) are
#' removed.
#'
#' Sampling is random; call \code{set.seed()} beforehand for reproducible
#' results.
#'
#' @return The downsampled data frame.
#'
#' @seealso \code{\link{get_list}}
#'
#' @importFrom dplyr filter select where %>%
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' df_balanced <- downsampling_age_dist(otu_tbl_df)
#' table(df_balanced$Lifestyle, df_balanced$age > 365)
#' }
#'
#' @export
downsampling_age_dist <- function(otu_tbl_df) {
  stds <- get_list("study", otu_tbl_df)
  df_red <- otu_tbl_df %>%
    filter(study %in% stds$non_industrialized | study %in% stds$industrialized) # subsample studies
  sbjcts <- get_list("subject_ID", otu_tbl_df, study_prob = T)
  df_red_sb <- df_red %>%
    filter(subject_ID %in% sbjcts$non_industrialized | subject_ID %in% sbjcts$industrialized) # subsample individuals
  df_red_y <- filter(df_red_sb, age <= 365) # same distribuition for above and below one year respectively
  df_red_o <- filter(df_red_sb, age > 365)
  smpls_y <- get_list("sample_ID", df_red_y, study_prob = T)
  df_red_y_smp <- df_red_y %>%
    filter(sample_ID %in% smpls_y$non_industrialized | sample_ID %in% smpls_y$industrialized)
  if(length(table(df_red_o$Lifestyle)) > 1) {
    smpls_o <- get_list("sample_ID", df_red_o, study_prob = T)
    df_red_o_smp <- df_red_o %>%
      filter(Lifestyle == "non_industrialized" | sample_ID %in% smpls_o$industrialized) %>%
      filter(Lifestyle == "industrialized" | sample_ID %in% smpls_o$non_industrialized) # subsample samples
    df_red_smp <- rbind(df_red_y_smp, df_red_o_smp)
  } else {
    print("No old samples found for one Lifestyle")
    df_red_smp <- df_red_y_smp
  }
  df_red_smp <- df_red_smp %>%
    select(where(~ !is.numeric(.x) || sum(.x, na.rm = TRUE) > 0))
  return(df_red_smp)
}


#' Count taxa absent from each lifestyle
#'
#' Sums the counts of each taxon per lifestyle and returns, for each
#' lifestyle, the number of taxa with a total count of zero. If taxa absent
#' from the whole dataset were removed beforehand (as in
#' \code{\link{downsampling_age_dist}}), a taxon absent from one lifestyle is
#' unique to the other.
#'
#' @param df A data frame with one row per sample, the columns
#'   \code{Lifestyle}, \code{study}, \code{subject_ID}, \code{sample_ID},
#'   \code{age} and \code{sample_sum}, and one numeric column per taxon.
#'
#' @return A named numeric vector with one value per lifestyle: the number
#'   of taxa not observed in that lifestyle.
#'
#' @seealso \code{\link{downsampling_age_dist}}
#'
#' @importFrom dplyr select group_by summarise across everything %>%
#' @importFrom tibble column_to_rownames
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' unique_taxa_per_ls(downsampling_age_dist(otu_tbl_df))
#' }
#'
#' @export
unique_taxa_per_ls <- function(df) {
  df %>% select(-c("study", "subject_ID", "sample_ID", "age", "sample_sum")) %>%
    group_by(Lifestyle) %>%
    summarise(across(everything(), ~ sum(.), na.rm = T)) %>%
    column_to_rownames("Lifestyle") %>%
    `==` (0) %>%
    rowSums()
}



#' Convert a distance object to a long data frame of sample pairs
#'
#' Turns pairwise distances into a data frame with one row per sample pair,
#' adds metadata for both samples, and flags whether the two samples share
#' age category, lifestyle and study. Pairs from the same subject are
#' removed.
#'
#' @param dist_obj A \code{dist} object (or square distance matrix) with
#'   sample IDs as labels, e.g. from \code{vegan::vegdist()} or
#'   \code{phyloseq::distance()}.
#' @param mdata A data frame with the columns \code{run_accession} (matching
#'   the labels of \code{dist_obj}), \code{subject_ID}, \code{age},
#'   \code{age_cat}, \code{study} and \code{lifestyle}.
#'
#' @details
#' Each pair is kept once. \code{sample_1} is the ID that sorts later
#' alphabetically. Metadata columns get the suffixes \code{_1} and \code{_2}
#' for the two samples.
#'
#' @return A data frame with the columns \code{sample_1}, \code{sample_2},
#'   \code{distance}, the metadata of both samples, and the logical columns
#'   \code{same_age_cat}, \code{same_ls} and \code{same_study}.
#'
#' @importFrom dplyr mutate filter left_join join_by %>%
#'
#' @examples
#' \dontrun{
#' bc <- phyloseq::distance(ps, method = "bray")
#' dist_df <- get_dist_df(bc, data.frame(sample_data(ps)))
#' }
#'
#' @export
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


#' Plot a PCoA ordination with optional ellipses and marginal densities
#'
#' A customised version of \code{phyloseq::plot_ordination()} for PCoA
#' results. Plots samples on two ordination axes with the percentage of
#' variance explained in the axis titles, and can add group ellipses and
#' density curves along both axes.
#'
#' @param ps_object A \code{phyloseq} object whose sample data contains the
#'   variables used for \code{color}, \code{shape} and \code{label}.
#' @param ordination A PCoA result from \code{phyloseq::ordinate()} (e.g.
#'   \code{method = "PCoA"}), containing eigenvalues in
#'   \code{$values$Eigenvalues}.
#' @param color,shape,label Character. Names of sample data columns mapped to
#'   point colour, point shape and text labels. Columns that don't exist are
#'   ignored with a warning. Default \code{NULL} (not mapped).
#' @param ellipses Logical. If \code{TRUE}, adds a 95 \% confidence ellipse
#'   per colour (and shape) group. Default \code{FALSE}.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param alpha,size Transparency and size of the points. Defaults 1.
#' @param axes Integer vector of length 2: the ordination axes to plot.
#'   Default \code{1:2}.
#' @param plot_dens Logical. If \code{TRUE}, adds density curves of the
#'   samples along both axes, filled by \code{color}. Requires \code{color}.
#'   Default \code{FALSE}.
#'
#' @details
#' Axis titles show the percentage of variance explained, calculated from the
#' eigenvalues of the selected axes relative to the sum of all eigenvalues.
#'
#' @return A \code{ggplot} object. With \code{plot_dens = TRUE}, a cowplot
#'   drawing instead, which can no longer be modified with \code{+} like a
#'   regular ggplot.
#'
#' @seealso \code{\link[phyloseq]{plot_ordination}},
#'   \code{\link[phyloseq]{ordinate}}
#'
#' @importFrom phyloseq plot_ordination
#' @importFrom ggplot2 ggplot aes_string geom_point geom_text ggtitle xlab
#'   ylab stat_ellipse geom_density coord_flip
#' @importFrom cowplot axis_canvas insert_xaxis_grob insert_yaxis_grob ggdraw
#' @importFrom rlang sym !!
#'
#' @examples
#' \dontrun{
#' ord <- phyloseq::ordinate(ps, method = "PCoA", distance = "bray")
#' plot_pcoa(ps, ord, color = "lifestyle", ellipses = TRUE, plot_dens = TRUE)
#' }
#'
#' @export
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
    p = p + geom_text(label_map, data = phyloseq:::rm.na.phyloseq(DF, 
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


#' Balance lifestyles by subsampling and compare beta diversity
#'
#' Performs one round of subsampling so that industrialized and
#' non-industrialized samples are balanced at study, subject and sample level,
#' then compares pairwise distances within and between lifestyles. Depending
#' on the arguments, it returns statistics, draws plots, or both. Intended to
#' be called repeatedly (e.g. in a loop) to assess the variability across
#' subsamples.
#'
#' @param df A data frame of sample metadata with the columns \code{study},
#'   \code{subject_ID}, \code{sample_ID}, \code{run_accession}, \code{age}
#'   (in days), \code{age_cat} and \code{lifestyle} (values
#'   \code{"industrialized"} and \code{"non_industrialized"}).
#' @param dist_obj A \code{dist} object with all samples, labelled by
#'   \code{run_accession}.
#' @param ps_obj A \code{phyloseq} object with the sample data used to colour
#'   the PCoA plots.
#' @param max_studies,max_subj,max_samples Maximum numbers of studies,
#'   subjects and samples per lifestyle, passed to \code{\link{get_list}}.
#'   \code{max_samples} is split equally between samples up to and above 365
#'   days. Default \code{NA} (only balance the lifestyles).
#' @param stats_output Logical. If \code{TRUE} (default), calculates and
#'   returns distance statistics.
#' @param plots Logical. If \code{TRUE} (default), draws a 2 × 2 panel with
#'   PCoA plots coloured by lifestyle and by age, the age distribution per
#'   lifestyle, and the within-lifestyle distance distributions.
#'
#' @details
#' Subsampling follows the same steps as \code{\link{downsampling_age_dist}}:
#' studies, then subjects (weighted by study), then samples (weighted by
#' study), with samples balanced separately up to and above 365 days. If only
#' one lifestyle has samples above 365 days, only the younger samples are
#' used.
#'
#' The distances are then reduced to the selected samples and converted to
#' pairs with \code{\link{get_dist_df}}, which removes pairs from the same
#' subject. Only pairs from the same lifestyle or the same age category are
#' kept. Within-lifestyle distances are compared between the lifestyles with
#' a Wilcoxon test. Between-lifestyle summaries use only pairs from different
#' studies.
#'
#' Subsampling is random; call \code{set.seed()} beforehand for reproducible
#' results.
#'
#' @return If \code{stats_output = TRUE}, a list with:
#'   \describe{
#'     \item{within_ls}{A one-row data frame with the Wilcoxon test result
#'       for within-lifestyle distances, plus the mean and median distance
#'       per lifestyle.}
#'     \item{between_ls_res}{Mean and median distances per age category,
#'       for pairs from the same and from different lifestyles.}
#'   }
#'   Otherwise \code{NULL}. Plots are drawn as a side effect and not
#'   returned.
#'
#' @seealso \code{\link{get_list}}, \code{\link{get_dist_df}},
#'   \code{\link{plot_pcoa}}, \code{\link{downsampling_age_dist}}
#'
#' @importFrom dplyr filter group_by summarise %>%
#' @importFrom tidyr pivot_longer pivot_wider unite
#' @importFrom rstatix wilcox_test
#' @importFrom ape pcoa
#' @importFrom gridExtra grid.arrange
#' @importFrom stats as.dist median
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' res <- purrr::map(1:100, ~ repeated_subsampling(meta_df, bray_dist, ps,
#'                                                 max_samples = 500,
#'                                                 plots = FALSE))
#' }
#'
#' @export
repeated_subsampling <- function(df, dist_obj, ps_obj, max_subj = NA,
                                 max_samples = NA, max_studies = NA,
                                 stats_output = T, plots = T) {
  stds <- get_list("study", df, max = max_studies)
  df_red <- df %>%
    filter(lifestyle == "non_industrialized" | study %in% stds$industrialized) %>%
    filter(lifestyle == "industrialized" | study %in% stds$non_industrialized) # subsample studies
  sbjcts <- get_list("subject_ID", df_red, max = max_subj, study_prob = T)
  df_red_sb <- df_red %>%
    filter(lifestyle == "non_industrialized" | subject_ID %in% sbjcts$industrialized) %>%
    filter(lifestyle == "industrialized" | subject_ID %in% sbjcts$non_industrialized) # subsample individuals
  df_red_y <- filter(df_red_sb, age <= 365) # same distribuition for above and below one year respectively
  df_red_o <- filter(df_red_sb, age > 365)
  smpls_y <- get_list("sample_ID", df_red_y, max = max_samples / 2, study_prob = T)
  df_red_y_smp <- df_red_y %>%
    filter(lifestyle == "non_industrialized" | sample_ID %in% smpls_y$industrialized) %>%
    filter(lifestyle == "industrialized" | sample_ID %in% smpls_y$non_industrialized)
  if(length(table(df_red_o$lifestyle)) > 1) {
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
  if(stats_output) {
    dist_df_long <- get_dist_df(dist_obj = dist_obj_filter, mdata = df_red) %>%
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

