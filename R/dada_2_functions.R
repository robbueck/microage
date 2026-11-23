##Project: Microbiome age index
##Aim: Additional functions for dada2 scripts
##Author: Robert Bücking

library(ape)
library(btools)
# library(rBLAST)

#' Add rarefied ASV-level alpha diversity to a phyloseq object
#'
#' Removes mitochondrial ASVs, rarefies the samples, and
#' calculates ASV richness, Shannon diversity and the number of unclassified
#' genera per sample. The results are added as new columns to the sample data
#' of the original (non-rarefied) phyloseq object.
#'
#' @param ps_object A \code{phyloseq} object. The taxonomy table must contain
#'   a rank called \code{genus} (lower case), the sample data must contain a
#'   column \code{run_accession} with unique sample IDs, and the OTU table must
#'   have taxa as columns (\code{taxa_are_rows = FALSE}).
#' @param keep_low_counts \code{logical} Whether to keep samples below the
#'  rarefaction threshold.
#' @param rarefaction_depth \code{numeric} Read depths for rarefaction to 
#'  calculate diversity indices.
#'
#' @details
#' \itemize{
#'   \item ASVs with \code{genus == "Mitochondria"} are removed before all
#'     calculations. Note that \code{subset_taxa()} also drops ASVs whose genus
#'     is \code{NA}.
#'   \item Samples with at least 2000 reads are rarefied to 2000 reads with
#'     \code{phyloseq::rarefy_even_depth()}. Samples with fewer reads are
#'     \strong{not} rarefied; their diversity is calculated on the raw counts.
#'   \item Rarefaction is random. Call \code{set.seed()} before this function
#'     for reproducible results.
#'   \item ASVs are counted as unclassified if their genus contains the string
#'     \code{"unclassified"}.
#'   \item Existing columns \code{shannon_asv}, \code{richness_asv} and
#'     \code{n_unassigned_asvs_sample_raref} are overwritten.
#' }
#'
#' @return The input \code{phyloseq} object (with mitochondrial ASVs removed)
#'   with three additional sample data columns:
#'   \describe{
#'     \item{shannon_asv}{Shannon diversity of the rarefied ASV counts.}
#'     \item{richness_asv}{Number of observed ASVs after rarefaction.}
#'     \item{n_unassigned_asvs_sample_raref}{Number of ASVs with an
#'       unclassified genus present in the rarefied sample.}
#'   }
#'
#' @importFrom phyloseq subset_taxa sample_sums prune_samples
#'   rarefy_even_depth merge_phyloseq sample_data
#' @importFrom dplyr filter select mutate left_join %>%
#' @importFrom tibble column_to_rownames
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' ps <- add_asv_alpha_div(ps)
#' head(phyloseq::sample_data(ps)[, c("shannon_asv", "richness_asv")])
#' }
#'
#' @export
add_asv_alpha_div <- function(ps_object, keep_low_counts = T, rarefaction_depth = 2000) {
  ps_object <- subset_taxa(ps_object, genus != "Mitochondria")
  if(min(sample_sums(ps_object)) < rarefaction_depth & keep_low_counts) {
    ps_low_counts <- ps_object %>% prune_samples(samples = (sample_sums(.) < rarefaction_depth))
    ps_raref <- ps_object %>% rarefy_even_depth(., sample.size = rarefaction_depth) %>%
      merge_phyloseq(.,ps_low_counts)
  } else {
    ps_raref <- ps_object %>% rarefy_even_depth(., sample.size = rarefaction_depth)
  }
  matching_taxa <- ps_raref@tax_table %>% data.frame() %>% filter(grepl("unclassified", genus)) %>% rownames()
  raref_otu_table <- ps_raref@otu_table %>% data.frame()
  unassigned_taxa_present_per_sample_raref <- apply(raref_otu_table[,matching_taxa], 1, function(x) sum(x > 0))
  ps_raref@sam_data$n_unassigned_asvs_sample_raref <- unassigned_taxa_present_per_sample_raref
  alpha.div_gen <- microbiome::alpha(ps_raref, index = c("Observed", "Shannon"))
  ps_raref@sam_data$shannon_asv <- alpha.div_gen$diversity_shannon
  ps_raref@sam_data$richness_asv <- alpha.div_gen$observed
  mdat_raref <- ps_raref@sam_data %>% data.frame() %>% select(shannon_asv, richness_asv, 
                                                              n_unassigned_asvs_sample_raref,
                                                              run_accession)
  ps_object@sam_data <- ps_object@sam_data %>% data.frame() %>%
    mutate(shannon_asv = NULL,
           richness_asv = NULL, 
           n_unassigned_asvs_sample_raref = NULL) %>%
    left_join(., mdat_raref, by = "run_accession") %>%
    mutate(rownames_1 = run_accession) %>%
    column_to_rownames("rownames_1") %>%
    sample_data()
  return(ps_object)
}


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


#' Calculate loss ratio
#'
#' Calculates loss ratio for different filtering steps.
#'
#' @param reduced \code{integer} Number of reads after filtering step.
#' @param total \code{integer} Number of reads before filtering step.
#'
#'
#' @return Depending on the input either a single value or vector of loss ratios.
#' @export
loss_ratio <- function(reduced, total) {
  res <- (total - reduced)/total
  res[is.na(res)] <- 1
  return(res)
}

#' Plot mean quality scores per position for a read cluster
#'
#' Calculates the mean and standard deviation of the quality scores at each
#' read position and plots them as a line with error bars.
#'
#' @param cluster_matrix A numeric matrix of quality scores, with one row per
#'   read and one column per position in the read.
#'
#' @details
#' Error bars show the mean ± one standard deviation at each position. Missing
#' values are not removed, so a position with any \code{NA} (e.g. from reads
#' of different lengths) gets \code{NA} as its mean and is not plotted.
#'
#' @return A \code{ggplot} object with position on the x-axis and mean quality
#'   score on the y-axis. It can be modified further with \code{+}, e.g. to
#'   add labels or a theme.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_errorbar
#'   position_dodge
#' @importFrom matrixStats colSds
#'
#' @examples
#' qual <- matrix(sample(20:40, 50 * 100, replace = TRUE), nrow = 50)
#' plot_cluster_qual(qual)
#'
#' @export
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



#' Wrapper to filter and trim paired-end FASTQ files with dada2
#'
#' Checks that the decontaminated paired-end FASTQ files exist for all
#' samples, filters and trims them with \code{dada2::filterAndTrim()}, and
#' writes the filtered files. Quality profiles can optionally be plotted
#' before and after filtering.
#'
#' @param sample_names Character vector of sample prefixes. Input files are
#'   expected as \code{<prefix>_dec_1.fastq.gz} and
#'   \code{<prefix>_dec_2.fastq.gz}.
#' @param plotting Logical. If \code{TRUE}, quality profiles of the first two
#'   samples are plotted before and after filtering. Default \code{FALSE}.
#' @param libname Currently not used.
#' @param trimLeft,truncLen,maxN,maxEE,truncQ Filtering parameters passed to
#'   \code{\link[dada2]{filterAndTrim}}. Give two values (forward, reverse)
#'   where dada2 allows it, e.g. \code{truncLen = c(240, 200)}.
#' @param cores Passed to \code{multithread} in \code{filterAndTrim()}:
#'   \code{TRUE} uses all available cores, a number sets the core count.
#'   Default \code{TRUE}.
#' @param verbose Logical, passed to \code{filterAndTrim()}. Default
#'   \code{FALSE}.
#'
#' @details
#' The function uses a global variable \code{maindir}, which must be defined
#' before calling it. Files are read from and written to:
#' \itemize{
#'   \item Input: \code{<maindir>/fastq_files/decontaminated/}
#'   \item Output: \code{<maindir>/fastq_files/filtered/}, as
#'     \code{<prefix>_dec_filt_1.fastq.gz} and
#'     \code{<prefix>_dec_filt_2.fastq.gz}
#' }
#' PhiX reads are removed (\code{rm.phix = TRUE}) and the output is
#' compressed. The function stops with an error if any input file is missing,
#' and warns if any output file is missing afterwards. dada2 does not write a
#' file for a sample in which no reads pass the filter, so a missing output
#' usually means that sample was filtered out completely.
#'
#' @return A matrix with one row per sample and the columns \code{reads.in}
#'   and \code{reads.out}, as returned by \code{filterAndTrim()}.
#'
#' @importFrom dada2 filterAndTrim plotQualityProfile
#'
#' @examples
#' \dontrun{
#' maindir <- "/path/to/project"
#' out <- filtering_fastq(c("sample1", "sample2"),
#'                        trimLeft = c(10, 10), truncLen = c(240, 200),
#'                        maxN = 0, maxEE = c(2, 2), truncQ = 2)
#' }
#'
#' @export
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


#' Wrapper to learn error rates and infer ASVs from filtered paired-end reads
#'
#' Runs the core dada2 workflow on filtered paired-end FASTQ files: learns
#' error rates, infers sequence variants for forward and reverse reads,
#' merges the read pairs and builds a sequence table.
#'
#' @param sample_names Character vector of sample prefixes. Input files are
#'   expected as \code{<prefix>_dec_filt_1.fastq.gz} and
#'   \code{<prefix>_dec_filt_2.fastq.gz}.
#' @param libname Character. Library name used in the file names of saved
#'   plots. Default \code{"all"}.
#' @param plotting Logical. If \code{TRUE}, saves the forward error-rate plot
#'   and the per-position quality of the first reverse sample as PDFs.
#'   Default \code{TRUE}.
#' @param debug_mode Logical. If \code{TRUE}, prints per sample the fraction
#'   of denoised forward reads that were merged. Default \code{FALSE}.
#' @param pool Pooling mode passed to \code{\link[dada2]{dada}}: \code{FALSE},
#'   \code{TRUE} or \code{"pseudo"}. Default \code{"pseudo"}.
#' @param nbases Minimum number of bases used to learn error rates, passed to
#'   \code{\link[dada2]{learnErrors}}. One value is used for both directions;
#'   two values set forward and reverse separately. Default \code{1e8}.
#' @param remove_empty Logical. If \code{FALSE}, the function stops when any
#'   filtered file is missing. If \code{TRUE}, missing files (e.g. samples
#'   with no reads left after filtering) are skipped. Default \code{FALSE}.
#' @param cores \code{TRUE}, \code{FALSE} or a number of cores; converted with
#'   \code{\link{get_n_cores}}. Default \code{TRUE} (cores of the SGE job).
#'
#' @details
#' The function uses a global variable \code{maindir}, which must be defined
#' before calling it. Files are read from
#' \code{<maindir>/fastq_files/filtered/}. With \code{plotting = TRUE}, plots
#' are saved to \code{<maindir>/dada2/} (the folder must exist) as
#' \code{<libname>_errors.pdf} and \code{lib_<libname>_seq_qual.pdf}.
#'
#' The function calls \code{set.seed(100)} before learning error rates. This
#' makes results reproducible but also resets the random number generator
#' in the calling session.
#'
#' @return A sequence table: an integer matrix with one row per sample and
#'   one column per ASV (column names are the sequences), as returned by
#'   \code{\link[dada2]{makeSequenceTable}}.
#'
#' @seealso \code{\link{filtering_fastq}} for the preceding filtering step,
#'   \code{\link{plot_cluster_qual}}
#'
#' @importFrom dada2 learnErrors plotErrors dada mergePairs getUniques
#'   makeSequenceTable
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' \dontrun{
#' maindir <- "/path/to/project"
#' seqtab <- error_and_asv(c("sample1", "sample2"), libname = "run1",
#'                         remove_empty = TRUE)
#' }
#'
#' @export
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



#' Wrapper to build a neighbour-joining tree from ASV sequences
#'
#' Aligns the ASV sequences with DECIPHER, calculates maximum-likelihood
#' distances, and builds a neighbour-joining tree, e.g. for adding to a
#' phyloseq object for UniFrac distances or phylogenetic diversity.
#'
#' @param x An object containing ASV sequences that
#'   \code{\link[dada2]{getSequences}} accepts, typically a sequence table
#'   from \code{\link{error_and_asv}} or a character vector of sequences.
#' @param cores Number of processors used for the alignment, passed to
#'   \code{processors} in \code{\link[DECIPHER]{AlignSeqs}}. \code{NULL}
#'   (default) uses all available cores.
#'
#' @details
#' Distances are calculated with \code{phangorn::dist.ml()} using its default
#' substitution model (JC69). The tree is built with \code{phangorn::NJ()}
#' and is not optimised further. It is unrooted; for UniFrac, phyloseq picks
#' a random root unless you root it yourself, e.g. with
#' \code{phangorn::midpoint()}.
#'
#' Tip labels are the full ASV sequences, matching the column names of a
#' dada2 sequence table. The order of the tips differs from the order of the
#' input sequences. Alignment time grows quickly with the number of ASVs.
#'
#' @return An unrooted tree of class \code{phylo}.
#'
#' @seealso \code{\link{error_and_asv}}
#'
#' @importFrom dada2 getSequences
#' @importFrom DECIPHER AlignSeqs
#' @importFrom Biostrings DNAStringSet
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
#' tree <- build_tree(seqtab, cores = 8)
#' ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
#'                phy_tree(tree))
#' }
#'
#' @export
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


#' Convert IdTaxa output to a taxonomy matrix
#'
#' Converts the output of \code{\link[DECIPHER]{IdTaxa}} (class
#' \code{"Taxa"}) into a character matrix with one row per ASV and one column
#' per rank, analogous to the output of \code{dada2::assignTaxonomy()}.
#'
#' @param ids An object of class \code{"Taxa"} as returned by
#'   \code{DECIPHER::IdTaxa()}, with one element per ASV.
#'
#' @details
#' The ranks domain, phylum, class, order, family, genus and species are
#' extracted by name, so the rank names in the training set must match these
#' (lower case). Ranks below the lowest assigned rank are filled with the
#' name of that lowest assigned taxon. For example, an ASV classified only to
#' family gets the family name in the genus and species columns too. Missing
#' ranks between two assigned ranks stay \code{NA}.
#'
#' Labels such as \code{"unclassified_Lachnospiraceae"} produced by
#' \code{IdTaxa()} are kept as they are.
#'
#' The row names are taken from the global object \code{st.all.nochim}
#' (the chimera-free sequence table), which must exist and contain the ASVs
#' in the same order as \code{ids}.
#'
#' @return A character matrix with the columns \code{domain}, \code{phylum},
#'   \code{class}, \code{order}, \code{family}, \code{genus} and
#'   \code{species}, and the ASV sequences as row names.
#'
#' @seealso \code{\link[DECIPHER]{IdTaxa}}, \code{\link{build_tree}}
#'
#' @importFrom dada2 getSequences
#' @importFrom utils tail
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' ids <- DECIPHER::IdTaxa(DNAStringSet(getSequences(st.all.nochim)),
#'                         trainingSet, processors = NULL)
#' taxa <- taxa_to_matrix(ids)
#' }
#'
#' @export
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



#' Merge phyloseq objects and combine ASV clusters
#'
#' Merges a list of phyloseq objects into one and combines groups of ASVs
#' (clusters) into single taxa. For each cluster, counts are summed and the
#' taxonomy of the first ASV is kept. Optionally, a consensus sequence is
#' calculated as the representative sequence of each cluster.
#'
#' @param ps_list A list of \code{phyloseq} objects, each with an OTU table,
#'   taxonomy table and sample data. Unless \code{simple_merge = TRUE}, each
#'   must also contain reference sequences (\code{refseq}).
#' @param clusters A list of character vectors. Each vector contains the
#'   names of ASVs to be merged into one taxon; all names must exist in the
#'   phyloseq objects. Ignored if \code{simple_merge = TRUE}.
#' @param consensus_seq Logical. If \code{TRUE} (default), each cluster is
#'   represented by the consensus sequence of its ASVs, calculated with
#'   \code{get_consensus_seq()}. If \code{FALSE}, the sequence of the first
#'   ASV in the cluster is kept, which is much faster but not suitable for
#'   phylogenetic analyses afterwards (e.g. trees, UniFrac, phylogenetic
#'   diversity).
#' @param cores Number of cores. Up to 4 are used as parallel workers via
#'   \pkg{future}; the full number is passed to \code{get_consensus_seq()}.
#'   Default 1.
#' @param simple_merge Logical. If \code{TRUE}, the phyloseq objects are only
#'   merged with \code{phyloseq::merge_phyloseq()}, without merging clusters
#'   or reference sequences. Taxa are combined only where their names are
#'   identical across objects. Default \code{FALSE}.
#'
#' @details
#' The OTU tables, taxonomy tables and sample data are merged separately and
#' then combined. Clusters are merged with \code{phyloseq::merge_taxa()} using
#' the first ASV of each cluster as the archetype. After merging, the function
#' checks that the total counts per sample are unchanged and stops otherwise.
#'
#' Unless \code{simple_merge = TRUE}, the taxa are renamed to \code{ASV1},
#' \code{ASV2}, ... in the result; the sequences remain available through
#' \code{refseq()}.
#'
#' The function calls \code{future::plan(multisession)}, which changes the
#' parallel setup of the calling session. It prints progress messages and the
#' unique values of the \code{study} column, so the sample data must contain
#' a column \code{study}.
#'
#' @return A merged \code{phyloseq} object. Unless \code{simple_merge = TRUE},
#'   it contains one taxon per cluster plus all unclustered ASVs, and
#'   reference sequences.
#'
#' @seealso \code{\link[phyloseq]{merge_phyloseq}},
#'   \code{\link[phyloseq]{merge_taxa}}
#'
#' @importFrom phyloseq otu_table sample_data tax_table refseq merge_phyloseq
#'   merge_taxa taxa_names ntaxa
#' @importFrom future plan multisession %<-%
#' @importFrom purrr reduce
#'
#' @examples
#' \dontrun{
#' clusters <- list(c("ASV1", "ASV17"), c("ASV5", "ASV8", "ASV23"))
#' ps_all <- merge_multiple_clusters(list(ps_study1, ps_study2),
#'                                   clusters = clusters, cores = 4)
#'
#' # fast version without consensus sequences
#' ps_all <- merge_multiple_clusters(list(ps_study1, ps_study2),
#'                                   clusters = clusters,
#'                                   consensus_seq = FALSE)
#' }
#'
#' @export
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


# fix merge_taxa2 naming, obsolete with update of microbiome package
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
