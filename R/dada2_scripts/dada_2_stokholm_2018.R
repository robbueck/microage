rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/studies/16S/stokholm_2018"
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
error_and_asv_step <- F
remove_chimeras_step <- F
id_taxa_step <- F
tree_step <- F
phyloseq_step <- T
ordination_step <- T


# load metadata
metadata_stokholm <- read.table("metadata_stokholm_2018_healthy.csv", sep = ",", header = T)

# metadata_stokholm <- metadata_stokholm[c(T, F, F, F, F, F, F, F, F, F, F),]
# filter, learn errors and denoise
print(paste("analyzing", length(metadata_stokholm$run_accession), sep = " "))
print("Filter")
if (filter_step) {
  out <- filtering_fastq(sample_names = metadata_stokholm$run_accession,
                     trimLeft = 25,
                     truncLen=c(220,170),
                     maxN=0,
                     maxEE=c(4,4), # to be checked, if this changes something
                     truncQ=6)  # also check this, if increasing lowers the chimera amount
  saveRDS(out, "dada2/stokholm_out_table.rds")
} else {
  out <- readRDS("dada2/stokholm_out_table.rds")
}
if (error_and_asv_step) {
  total_bases <- c(sum(out[,2]) * 210, sum(out[,2]) * 160)
  st.all <- error_and_asv(sample_names = metadata_stokholm$run_accession, nbases = total_bases*0.25)
  saveRDS(st.all, "dada2/stokholm_seqtab.rds")
} else {
  print("D")
  st.all <- readRDS("dada2/stokholm_seqtab.rds")
}
print("E")

# how much got filtered out?
cdf_data_filter <- melt(loss_ratio(out[,2], out[,1]),
                        value.name = "filter_ratio")
ggplot(cdf_data_filter, aes(filter_ratio)) +
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
ggsave("dada2/cfd_filter_ratios.pdf", device = "pdf")

# assess filter criteria for samples from each library
fnFs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_stokholm$run_accession, "_dec_1.fastq.gz")))
fnRs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_stokholm$run_accession, "_dec_2.fastq.gz")))
filtFs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_stokholm$run_accession, "_dec_filt_1.fastq.gz")))
filtRs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_stokholm$run_accession, "_dec_filt_2.fastq.gz")))
plotQualityProfile(c(fnFs[1], filtFs[1], fnRs[1], filtRs[1]))
ggsave("dada2/quality_profile.pdf", device = "pdf")


# Inspect distribution of sequence lengths
cdf_data_seqlen <- melt(nchar(getSequences(st.all)),
                        value.name = "seq_length")
ggplot(cdf_data_seqlen, aes(seq_length)) +
  stat_ecdf() +
  ylab("samples")
ggsave("dada2/cfd_seq_lengths_per_library.pdf", device = "pdf")



# remove chimeras
print("remove chimeras")
if (remove_chimeras_step) {
  st.all.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=T, verbose=TRUE)
  saveRDS(st.all.nochim, "st_all_nochim.rds")
} else {
  st.all.nochim <- readRDS("st_all_nochim.rds")
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
print("assign taxonomy")
dna <- DNAStringSet(getSequences(st.all.nochim)) # Create a DNAStringSet from the ASVs
load("/fast/AG_Forslund/rob/references/silva/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
print("id taxa")
if (id_taxa_step) {
  ids <- IdTaxa(dna, trainingSet, processors=NULL, verbose=FALSE)
  saveRDS(ids, "taxa_ids.rds")
} else {
  ids <- readRDS("taxa_ids.rds")
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
if ( phyloseq_step) {
  metadata_stokholm <- read.table("metadata_stokholm_2018_healthy.csv", sep = ",", header = T)
  #add metadata
  head(metadata_stokholm)
  rownames(metadata_stokholm) <- metadata_stokholm$run_accession
  all(rownames(metadata_stokholm) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata_stokholm))
  metadata_stokholm <- metadata_stokholm[match(rownames(st.all.nochim), rownames(metadata_stokholm)),]
  all(rownames(metadata_stokholm) == rownames(st.all.nochim))
  # number of ASVs unassigned in whole study
  metadata_stokholm$n_unassigned_asvs_study <-   sum(grepl("unclassified", taxa[,6]))
  metadata_stokholm$n_total_asvs_study <- nrow(taxa)
  # number of ASVs not assigned to genus level per sample
  matching_taxa <- rownames(taxa)[grepl("unclassified", taxa[,6])]
  unassigned_taxa_present_per_sample <- apply(st.all.nochim[,matching_taxa], 1, function(x) sum(x > 0))
  total_taxa_present_per_sample <- apply(st.all.nochim, 1, function(x) sum(x > 0))
  metadata$n_unassigned_asvs_sample <- unassigned_taxa_present_per_sample
  metadata$n_total_asvs_sample <- total_taxa_present_per_sample

  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata_stokholm), 
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
  saveRDS(ps, file = "/fast/AG_Forslund/rob/studies/16S/stokholm_2018/dada2/phyloseq_stokholm_2018.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/studies/16S/stokholm_2018/dada2/phyloseq_stokholm_2018.rds")
}


# compare a-diversity between different groups
diversity_plots(ps_object = ps, x = "age", colour = "age")
ggsave("dada2/a_div_age.pdf", device = "pdf")

# check metadata in PCoA
# remove zero count samples
ps <- prune_samples(sample_sums(ps) >= 10, ps)
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
print("ordination step")
if(ordination_step) {
  dist.bray <- phyloseq::distance(ps.prop, method = "bray")
  dist.unifrac <- phyloseq::distance(ps.prop, method = "wunifrac")
  # ordinate
  pcoa_bray <- pcoa(dist.bray)
  pcoa_unifr <- pcoa(dist.unifrac)
  tsne_bray <- Rtsne(dist.bray, is_distance = TRUE, perplexity = 80, num_threads = 10, dims = 3)
  tsne_unifrac <- Rtsne(dist.unifrac, is_distance = TRUE)
  umap <- umap(otu_table(ps.prop))
  save(pcoa_bray, pcoa_unifr, tsne_bray, tsne_unifrac, umap, file = "dada2/ordinations.RData")
} else {
  load("dada2/ordinations.RData")
}

rownames(metadata_stokholm) <- metadata_stokholm$run_accession
cor(pcoa_bray$vectors,
    metadata_stokholm[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head
cor(pcoa_unifr$vectors,
    metadata_stokholm[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head


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

# t-SNE
tsnedata <- data.frame(tsne_bray$Y) %>%
  dplyr::rename(tSNE1 = X1, tSNE2 = X2) %>%
  bind_cols(data.frame(sample_data(ps)))
ps1 <- ggplot(tsnedata, aes(x = tSNE1, y = tSNE2, color = age)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) + 
  ggtitle("t-SNE Bray")

data.frame(tsne_bray$Y) %>%
  dplyr::rename(tSNE1 = X1, tSNE2 = X2, tSNE3 = X3) %>%
  bind_cols(data.frame(sample_data(ps))) %>%
  plot_ly(., x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~age)

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

# patterns to see: 


# export to be used with picrust:
writeXStringSet(ps@refseq, "dada2/ASVs.fasta", format = "fasta")
otu_table_out <- t(as.data.frame(ps@otu_table))
write.table(data.frame(OTU = rownames(otu_table_out), otu_table_out), file = "dada2/otu_table.tsv", sep = "\t", row.names = F, quote = F)

