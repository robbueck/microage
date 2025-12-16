rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/studies/16S/roswall_2021/"
setwd(maindir)
library(dada2)
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("phangorn")
library("picante")
library("btools")
library("microbiome")
library("umap")
library("gridExtra")
library("Rtsne")


# load additional custom functions
source("/fast/AG_Forslund/rob/mm_index/R_scripts/dada_2_functions.R")

# switch parts of the pipeline on/off
filter_step <- F
error_step <- F
asv_step <- F
merge_step <- F
remove_chimeras_step <- F
id_taxa_step <- F
tree_step <- F
phyloseq_step <- T
ordination_step <- F


# load metadata
metadata_roswall <- read.table("metadata_roswall_2021_healthy.csv", sep = ",", header = T)


# filter, learn errors and denoise
print(paste("analyzing", length(metadata_roswall$run_accession), sep = " "))
print("Filter")
if(filter_step) {
  out <- filtering_fastq(sample_names = metadata_roswall$run_accession,
                     trimLeft = 10,
                     truncLen=190,
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=2)
  saveRDS(out, "dada2/roswall_out_table.rds")
} else {
  out <- readRDS("dada2/roswall_out_table.rds")
}

total_bases <- sum(out[,2]) * 190
sample_names = metadata_roswall$run_accession
nbases = 0.25*total_bases
cores = T

fltFfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt_1.fastq.gz")))
fltRfs <- sort(file.path(maindir, "fastq_files/filtered", paste0(sample_names, "_dec_filt_2.fastq.gz")))
names(fltFfs) <- sort(sample_names)
names(fltRfs) <- sort(sample_names)
if(!all(file.exists(fltFfs, fltRfs))){
  print("Missing files:")
  print(fltFfs[!file.exists(fltFfs)])
  print(fltRfs[!file.exists(fltRfs)])
  stop("Filtered files not complete")
}

if(error_step) {
  # learn error rates
  set.seed(100)
  print("Lern Error-rates")

  errF <- learnErrors(fltFfs, nbases=nbases, multithread=cores)
  errR <- learnErrors(fltRfs, nbases=nbases, multithread=cores)
  plotErrors(errF, nominalQ=TRUE)
  ggsave(paste0(maindir, "/dada2/", "all_errors.pdf"), device = "pdf")
  saveRDS(errF, "dada2/roswall_errf.rds")
  saveRDS(errR, "dada2/roswall_errr.rds")
} else {
  errF <- readRDS("dada2/roswall_errf.rds")
  errR <- readRDS("dada2/roswall_errr.rds")
}

if(asv_step) {
  # sample inference
  print("sample inference")
  dadaFs <- dada(fltFfs, err=errF, multithread=cores, pool = "pseudo")
  dadaRs <- dada(fltRfs, err=errR, multithread=cores, pool = "pseudo")

  # quality along the sequence
  plot_cluster_qual(dadaRs[[1]]$quality)
  ggsave(paste0(maindir, "/dada2/all_seq_qual.pdf"), device = "pdf")
  saveRDS(dadaFs, "dada2/roswall_dadafs.rds")
  saveRDS(dadaRs, "dada2/roswall_dadars.rds")
} else {
  dadaFs <- readRDS("dada2/roswall_dadafs.rds")
  dadaRs <- readRDS("dada2/roswall_dadars.rds")
}

if(merge_step) {
  # merge reads
  mergers <- mergePairs(dadaFs, fltFfs, dadaRs, fltRfs, verbose=TRUE)
  # construct sequence table
  st.all <- makeSequenceTable(mergers)
  saveRDS(st.all, "dada2/roswall_seqtab.rds")
} else {
  st.all <- readRDS("dada2/roswall_seqtab.rds")
}


# how much got filtered out?
head(out)
cdf_data_filter <- melt(loss_ratio(out[,2], out[,1]),
                        value.name = "filter_ratio")
ggplot(cdf_data_filter, aes(filter_ratio)) +
  stat_ecdf() +
  ylab("samples")+
  xlim(1,0)
ggsave("dada2/cfd_filter_ratios.pdf", device = "pdf")

# assess filter criteria for samples from each library
# fnFs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_roswall$run_accession, "_dec_1.fastq.gz")))
# fnRs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_roswall$run_accession, "_dec_2.fastq.gz")))
# plotQualityProfile(c(fnFs[1], fnRs[1]))
# filtFs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_roswall$run_accession, "_dec_filt_1.fastq.gz")))
# filtRs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_roswall$run_accession, "_dec_filt_2.fastq.gz")))
# plotQualityProfile(c(filtFs[1], filtRs[1]))
# plotQualityProfile(c(fnFs[1], filtFs[1], fnRs[1], filtRs[1]))
# ggsave("dada2/quality_profile.pdf", device = "pdf")


# Inspect distribution of sequence lengths
cdf_data_seqlen <- melt(nchar(getSequences(st.all)),
                        value.name = "seq_length")
ggplot(cdf_data_seqlen, aes(seq_length)) +
  stat_ecdf() +
  ylab("samples")
ggsave("dada2/cfd_seq_lengths_per_library.pdf", device = "pdf")



# remove chimeras
print("remove chimeras")
if(remove_chimeras_step) {
  st.all.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=T, verbose=TRUE)
  saveRDS(st.all.nochim, "dada2/st_all_nochim.rds")
} else {
  st.all.nochim <- readRDS("dada2/st_all_nochim.rds")
}
dim(st.all.nochim)
sum(st.all.nochim)/sum(st.all)


# where do most reads get lost?
# load read counts for decontamination
read_counts_raw <- read.table("read_counts_raw.txt", header = F,
                              col.names = c("sample_ID", "read_counts_raw"), row.names = 1)
rownames(out) <- gsub("_dec_1.fastq.gz", "", rownames(out))
loss_table <- merge(out, read_counts_raw, by = 0)
head(loss_table)
loss_table <- merge(loss_table,
                    cbind(denoised = rowSums(st.all), removed_chimeras = rowSums(st.all.nochim)),
                    by.x = "Row.names", by.y = 0)

track <- cbind(loss_ratio(loss_table$reads.in, loss_table$read_counts_raw),  # decontamination
               loss_ratio(loss_table$reads.out, loss_table$reads.in),        # filtering
               loss_ratio(loss_table$denoised, loss_table$reads.out),        # denoising
               loss_ratio(loss_table$removed_chimeras, loss_table$denoised), # remove chimeras
               loss_ratio(loss_table$removed_chimeras, loss_table$read_counts_raw))        # total
colnames(track) <- c("decontamination", "filtering", "denoising", "remove_chimeras", "total_loss")
rownames(track) <- loss_table$Row.names
head(track)

cdf_data_step_loss <- melt(track,
                 value.name = "lost_reads")
head(cdf_data_step_loss)
colnames(cdf_data_step_loss)[1:2] <- c("sample", "step")
# violin plot:
ggplot(cdf_data_step_loss, aes(x=step, y = lost_reads, color = step)) + 
  geom_violin() +
  ylim(0,1) +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/violin_filter_ratios_per_step.pdf", device = "pdf")

cdf_data_step_loss <- ddply(cdf_data_step_loss, .(step), transform, ecd=ecdf(lost_reads)(lost_reads))
ggplot(cdf_data_step_loss, aes(x=lost_reads, col = step)) +
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
cdf_data_step_loss[is.na(cdf_data_step_loss$lost_reads),]
ggsave("dada2/cfd_filter_ratios_per_step.pdf", device = "pdf")

# save loss table to be used for further stats
loss_table <- dplyr::transmute(loss_table, Row.names = Row.names, read_counts_raw = read_counts_raw,
                               decontaminated = reads.in, filtered = reads.out,
                               denoised = denoised, removed_chimeras = removed_chimeras)
write.csv(loss_table, "dada2/read_loss_per_step.csv")


# assign taxonomy
library(DECIPHER)
print("id taxa")
dna <- DNAStringSet(getSequences(st.all.nochim)) # Create a DNAStringSet from the ASVs
if(id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/silva/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids <- IdTaxa(dna, trainingSet, processors=NULL, verbose=FALSE)
  saveRDS(ids, "dada2/taxa_ids.rds")
} else {
  ids <- readRDS("dada2/taxa_ids.rds")
}

taxa <- taxa_to_matrix(ids)

if (tree_step) {
  tree <- build_tree(st.all.nochim)
  saveRDS(tree, "dada2/tree.rds")
} else {
  tree <- readRDS("dada2/tree.rds")
}


# handoff to phyloseq
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())
if(phyloseq_step) {
  metadata_roswall <- read.table("metadata_roswall_2021_healthy.csv", sep = ",", header = T)
  #add metadata
  head(metadata_roswall)
  rownames(metadata_roswall) <- metadata_roswall$run_accession
  all(rownames(metadata_roswall) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata_roswall))
  metadata_roswall <- metadata_roswall[match(rownames(st.all.nochim), rownames(metadata_roswall)),]
  all(rownames(metadata_roswall) == rownames(st.all.nochim))
  # number of ASVs unassigned in whole study
  metadata_roswall$n_unassigned_asvs_study <-   sum(grepl("unclassified", taxa[,6]))
  metadata_roswall$n_total_asvs_study <- nrow(taxa)
  # number of ASVs not assigned to genus level per sample
  matching_taxa <- rownames(taxa)[grepl("unclassified", taxa[,6])]
  unassigned_taxa_present_per_sample <- apply(st.all.nochim[,matching_taxa], 1, function(x) sum(x > 0))
  total_taxa_present_per_sample <- apply(st.all.nochim, 1, function(x) sum(x > 0))
  metadata$n_unassigned_asvs_sample <- unassigned_taxa_present_per_sample
  metadata$n_total_asvs_sample <- total_taxa_present_per_sample

  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata_roswall), 
                 tax_table(taxa), tree)
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  # rarefaction for alpha diversity
  if(min(sample_sums(ps)) < 2000) {
    ps_low_counts <- ps %>% prune_samples(samples = (sample_sums(.) < 2000))
    ps_raref <- ps %>% rarefy_even_depth(., sample.size = 2000, trimOTUs = F) %>%
      merge_phyloseq(.,ps_low_counts)
  } else {
    ps_raref <- ps %>% rarefy_even_depth(., sample.size = 2000)
  }
  # remove low abundant taxa:
  remove_taxa_asv <- ps_raref %>%
    prevalence(., detection = 10)
  remove_taxa_asv <- names(remove_taxa_asv[remove_taxa_asv < 5e-5])
  ps_raref <- ps_raref %>%
    merge_taxa2_fixed(.,taxa = remove_taxa_asv, name = "Other")
  # add alpha diversity on ASV level
  alpha.div_gen <- microbiome::alpha(ps_raref, index = c("Observed", "Shannon"))
  alpha.div_gen <- alpha.div_gen[sample_names(ps),]
  if (all(sample_names(ps) == rownames(alpha.div_gen))) {
    ps@sam_data$shannon_asv <- alpha.div_gen$observed
    ps@sam_data$richness_asv <- alpha.div_gen$diversity_shannon
  } else {stop(" Names don't match")}
  
  # save object
  saveRDS(ps, file = "/fast/AG_Forslund/rob/studies/16S/roswall_2021/dada2/phyloseq_roswall_2021.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/studies/16S/roswall_2021/dada2/phyloseq_roswall_2021.rds")
}


# compare a-diversity between different groups
diversity_plots(ps_object = ps, x = "age", colour = "age")
ggsave("dada2/a_div_age.pdf", device = "pdf")

# check metadata in PCoA
# remove zero count samples
ps <- prune_samples(sample_sums(ps) >= 10, ps)
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
if(ordination_step) {
  dist.bray <- phyloseq::distance(ps.prop, method = "bray")
  dist.unifrac <- phyloseq::distance(ps.prop, method = "wunifrac")
  # ordinate
  pcoa_bray <- pcoa(dist.bray)
  pcoa_unifr <- pcoa(dist.unifrac)
  tsne_bray <- Rtsne(dist.bray, is_distance = TRUE)
  tsne_unifrac <- Rtsne(dist.unifrac, is_distance = TRUE)
  umap <- umap(otu_table(ps.prop))
  save(dist.bray, dist.unifrac, pcoa_bray, pcoa_unifr, tsne_bray, tsne_unifrac, umap, file = "dada2/ordinations.RData")
} else {
  load("dada2/ordinations.RData")
}

rownames(metadata_roswall) <- metadata_roswall$run_accession
cor(pcoa_bray$vectors,
    metadata_roswall[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head
cor(pcoa_unifr$vectors,
    metadata_roswall[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head

# plot reduced dimensions
plot_ordination(ps.prop, pcoa_bray, color="age", title="Bray NMDS") +
  scale_color_gradientn(name = "age", colors = topo.colors(7))
#PCoA
# bray-curtis
p1 <- plot_ordination(ps.prop, pcoa_bray, color="age", title="Bray NMDS") +
  scale_color_gradientn(name = "age", colors = topo.colors(7))
# unifrac dist:
p2 <- plot_ordination(ps.prop, pcoa_unifr, color="age", title="Unifrac NMDS") +
  scale_color_gradientn(name = "age", colors = topo.colors(7))
pf <- grid.arrange(p1, p2)
ggsave("dada2/PCoA_age.pdf", device = "pdf", plot = pf)

fit_adonis <- adonis2(dist.bray ~ age + antibiotics_before + food + birthmode + sex,
                      data = metadata_roswall, by="margin", na.action = na.omit)
saveRDS(fit_adonis,
        file = "dada2/permanova.RData")
fit_adonis <- readRDS("dada2/permanova.RData")
print(fit_adonis)



# t-SNE
tsnedata <- data.frame(tsne_bray$Y) %>%
  dplyr::rename(tSNE1 = X1, tSNE2 = X2) %>%
  bind_cols(data.frame(sample_data(ps)))
ps1 <- ggplot(tsnedata, aes(x = tSNE1, y = tSNE2, color = age)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) + 
  ggtitle("t-SNE Bray")
# unifrac
tsnedata <- data.frame(tsne_unifrac$Y) %>%
  dplyr::rename(tSNE1 = X1, tSNE2 = X2) %>%
  bind_cols(data.frame(sample_data(ps)))
ps2 <- ggplot(tsnedata, aes(x = tSNE1, y = tSNE2, color = age)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) +
  ggtitle("t-SNE Unifrac")
psf <- grid.arrange(ps1, ps2, ncol = 2)
ggsave("dada2/tSNE.pdf", device = "pdf", plot = psf)

# UMAP
umapdata <- umap$layout %>%
  data.frame() %>%
  dplyr::rename(UMAP1 = X1, UMAP2 = X2) %>%
  bind_cols(data.frame(sample_data(ps)))
ggplot(umapdata, aes(x = UMAP1, y = UMAP2, color = age)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) +
  ggtitle("UMAP")
ggsave("dada2/umap.pdf", device = "pdf")


# export to be used with picrust:
writeXStringSet(ps@refseq, "dada2/ASVs.fasta", format = "fasta")
otu_table_out <- t(as.data.frame(ps@otu_table))
write.table(data.frame(OTU = rownames(otu_table_out), otu_table_out), file = "dada2/otu_table.tsv", sep = "\t", row.names = F, quote = F)

