rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/studies/16S/gehrig_2019/"
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
library("RColorBrewer")
library(tidyverse)



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
metadata_gehr <- read.table("metadata_gehrig_2019_healthy.csv", sep = ",", header = T)

# filter, learn errors and denoise
print(paste("analyzing", length(metadata_gehr$run_accession), sep = " "))
print("Filter")
if(filter_step) {
  out <- filtering_fastq(sample_names = metadata_gehr$run_accession,
                     trimLeft = 10,
                     truncLen=c(220, 150),
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=2,
                     cores = T)
  saveRDS(out, "dada2/gehrig_out_table.rds")
} else {
  out <- readRDS("dada2/gehrig_out_table.rds")
}
if(error_and_asv_step){
  total_bases <- c(sum(out[,2]) * 220, sum(out[,2]) * 150) # because error estimation is done for forward and reverse
  st.all <- error_and_asv(sample_names = metadata_gehr$run_accession, nbases = 0.25* total_bases)
  saveRDS(st.all, "dada2/gehrig_seqtab.rds")
} else {
  st.all <- readRDS("dada2/gehrig_seqtab.rds")
}
  


# how much got filtered out?
head(out)
cdf_data_filter <- melt(loss_ratio(out[,2], out[,1]),
                        value.name = "filter_ratio")
ggplot(cdf_data_filter, aes(filter_ratio)) +
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
ggsave("dada2/cfd_filter_ratios.pdf", device = "pdf")

# assess filter criteria for samples from each library
fnFs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_gehr$run_accession, "_dec_1.fastq.gz")))
fnRs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_gehr$run_accession, "_dec_2.fastq.gz")))
filtFs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_gehr$run_accession, "_dec_filt_1.fastq.gz")))
filtRs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_gehr$run_accession, "_dec_filt_2.fastq.gz")))
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
if(remove_chimeras_step){
  st.all.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=10, verbose=TRUE)
  saveRDS(st.all.nochim, "dada2/st_all_nochim.rds")
} else{
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
print("assign taxonomy")
dna <- DNAStringSet(getSequences(st.all.nochim)) # Create a DNAStringSet from the ASVs
start.time_silva <- Sys.time()
if(id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/silva/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids <- IdTaxa(dna, trainingSet, processors=10, verbose=FALSE)
  saveRDS(ids, "dada2/taxa_ids.rds")
} else {
  ids <- readRDS("dada2/taxa_ids.rds")
}
end.time_silva <- Sys.time()
time_silva <- end.time_silva - start.time_silva

taxa <- taxa_to_matrix(ids)

# start.time_gtdb <- Sys.time()
# if(id_taxa_step) {
#   load("/fast/AG_Forslund/rob/references/GTDB_r207-mod_April2022/GTDB_r207-mod_April2022.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
#   ids_gtdb <- IdTaxa(dna, trainingSet, processors=10, verbose=FALSE)
#   saveRDS(ids_gtdb, "dada2/taxa_ids_GTDB.rds")
# } else {
#   ids_gtdb <- readRDS("dada2/taxa_ids_GTDB.rds")
# }
# end.time_gtdb <- Sys.time()
# time_gtdb <- end.time_gtdb - start.time_gtdb
# 
# taxa_gtdb <- taxa_to_matrix(ids_gtdb)
# 
# taxdb_table <- t(data.frame(silva = (!is.na(taxa)) %>% colSums() / nrow(taxa),
#              gtdb = (!is.na(taxa_gtdb)) %>% colSums() / nrow(taxa_gtdb)))


if (tree_step) {
  tree <- build_tree(st.all.nochim)
  saveRDS(tree, "dada2/tree.rds")
} else {
  tree <- readRDS("dada2/tree.rds")
}


# compare compositional plots between databases:
metadata_gehr <- read.table("metadata_gehrig_2019_healthy.csv", sep = ",", header = T)
#add metadata
rownames(metadata_gehr) <- metadata_gehr$run_accession
# all(rownames(metadata_gehr) %in% rownames(st.all.nochim))
# all(rownames(st.all.nochim) %in% rownames(metadata_gehr))
metadata_gehr <- metadata_gehr[match(rownames(st.all.nochim), rownames(metadata_gehr)),]
# all(rownames(metadata_gehr) == rownames(st.all.nochim))
level <- "genus"

# SILVA
ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata_gehr), 
               tax_table(taxa), tree)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps <- aggregate_taxa(ps, level = level)
sample_names(ps) <- paste0(sample_names(ps), "_SILVA")
sample_data(ps)$annot_DB <- "SILVA"

#GTDB
# ps_GTDB <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
#                     sample_data(metadata_gehr), 
#                     tax_table(taxa_gtdb), tree)
# dna <- Biostrings::DNAStringSet(taxa_names(ps_GTDB))
# names(dna) <- taxa_names(ps_GTDB)
# ps_GTDB <- merge_phyloseq(ps_GTDB, dna)
# taxa_names(ps_GTDB) <- paste0("ASV", seq(ntaxa(ps_GTDB)))
# sample_names(ps_GTDB) <- paste0(sample_names(ps_GTDB), "_GTDB")
# sample_data(ps_GTDB)$annot_DB <- "GTDB"
# ps_GTDB <- aggregate_taxa(ps_GTDB, level = level)
# 
# 
# ps_combined <- merge_phyloseq(ps_GTDB, ps)
# ps_combined <- ps_combined %>%
#   subset_samples(age > 100) %>%
#   subset_taxa(get(level) != "Unknown") %>%
#   microbiome::transform(transform = "compositional") %>%
#   filter_taxa(function(x) mean(x) > 1e-3, TRUE)
# n_taxa <- ps_combined@tax_table[,level] %>% unique() %>% length()
# plot_composition(ps_combined,
#                  average_by = "annot_DB", 
#                  transform = "compositional") +
#   scale_fill_manual(level, values = colorRampPalette(brewer.pal(12, "Paired"))(n_taxa)) +
#   theme(axis.text.x = element_text(angle=90, hjust=1),
#         legend.text = element_text(face = "italic"))
# 
# 
# # compare b-div between the two dbs
# dist.bray <- phyloseq::distance(ps_combined, method = "bray")
# pcoa_bray <- pcoa(dist.bray)
# 
# plot_ordination(ps_combined, pcoa_bray, color="age", title="Bray NMDS", shape = "annot_DB") +
#   scale_color_gradientn(name = "age", colors = topo.colors(7))
# 
# # intra-vs-outra-group beta diversity:
# x<- as.matrix(dist.bray)
# 
# y <- t(combn(colnames(x), 2))
# BC.Inf<- data.frame(y, dist=x[y])
# 
# BC.Inf$Sample_pair<- paste(BC.Inf$X1, BC.Inf$X2, sep = "-")
# 
# 
# BC.Inf%>%
#   separate(X1, c("Baby_A", "DB_A"))%>%
#   separate(X2, c("Baby_B", "DB_B"))%>%
#   dplyr::mutate(Same_Individual = case_when(Baby_A == Baby_B  ~ T,
#                                             Baby_A != Baby_B ~ F))%>%
#   dplyr::mutate(Same_DB = case_when(DB_A == DB_B  ~ T,
#                                              DB_A != DB_B ~ F))%>%
#   dplyr::select(c("Sample_pair", "Baby_A", "Baby_B","Same_Individual", "Same_DB", "dist"))-> BC.Inf
# 
# library(rstatix)
# BC.Inf%>%
#   wilcox_test(dist ~ Same_Individual)%>%
#   add_significance()%>%
#   add_xy_position(x = "Same_Individual")-> stats.test 
# 
# BC.Inf%>%
#   wilcox_test(dist ~ Same_DB)%>%
#   add_significance()%>%
#   add_xy_position(x = "Same_DB")-> stats.test 
# 
# BC.Inf%>%
#   ggplot(aes(x= Same_DB, y= dist, fill= Same_DB))+
#   geom_boxplot(aes(),outlier.shape=NA)+
#   geom_point(position = position_jitterdodge(), alpha= 0.1)+
#   scale_color_manual(values = c("black", "black"))+
#   scale_fill_manual(values = c("#CC6677", "#117733"))+
#   ylab("Bray-Curtis intersample distances")+
#   guides(fill = FALSE, color= FALSE)+
#   theme_classic()+
#   theme(text = element_text(size=16), axis.title.x = element_blank())+
#   scale_x_discrete(labels=c("A" = "Different_annotation_DB", 
#                             "B" = "Same_annotation_DB"))+
#   scale_y_continuous(limits=c(0, 1.2))
# 
# 
# BC.Inf%>%
#   ggplot(aes(x= Same_Individual, y= dist, fill= Same_Individual))+
#   geom_boxplot(aes(),outlier.shape=NA)+
#   geom_point(position = position_jitterdodge(), alpha= 0.1)+
#   scale_color_manual(values = c("black", "black"))+
#   scale_fill_manual(values = c("#CC6677", "#117733"))+
#   ylab("Bray-Curtis intersample distances")+
#   guides(fill = FALSE, color= FALSE)+
#   theme_classic()+
#   theme(text = element_text(size=16), axis.title.x = element_blank())+
#   scale_x_discrete(labels=c("A" = "Different_annotation_DB", 
#                             "B" = "Same_annotation_DB"))+
#   scale_y_continuous(limits=c(0, 1.2))


# handoff to phyloseq
library(phyloseq)
library(Biostrings)
theme_set(theme_bw())
if(phyloseq_step) {
  metadata_gehr <- read.table("metadata_gehrig_2019_healthy.csv", sep = ",", header = T)
  #add metadata
  head(metadata_gehr)
  rownames(metadata_gehr) <- metadata_gehr$run_accession
  all(rownames(metadata_gehr) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata_gehr))
  metadata_gehr <- metadata_gehr[match(rownames(st.all.nochim), rownames(metadata_gehr)),]
  all(rownames(metadata_gehr) == rownames(st.all.nochim))
  # number of ASVs unassigned in whole study
  metadata_gehr$n_unassigned_asvs_study <-   sum(grepl("unclassified", taxa[,6]))
  metadata_gehr$n_total_asvs_study <- nrow(taxa)
  # number of ASVs not assigned to genus level per sample
  matching_taxa <- rownames(taxa)[grepl("unclassified", taxa[,6])]
  unassigned_taxa_present_per_sample <- apply(st.all.nochim[,matching_taxa], 1, function(x) sum(x > 0))
  total_taxa_present_per_sample <- apply(st.all.nochim, 1, function(x) sum(x > 0))
  metadata$n_unassigned_asvs_sample <- unassigned_taxa_present_per_sample
  metadata$n_total_asvs_sample <- total_taxa_present_per_sample

  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata_gehr), 
                 tax_table(taxa), tree)
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  # rarefaction for alpha diversity
  if(min(sample_sums(ps)) < 2000) {
    ps_low_counts <- ps %>% prune_samples(samples = (sample_sums(.) < 2000))
    ps_raref <- ps %>% rarefy_even_depth(., sample.size = 2000) %>%
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
  if (all(sample_names(ps) == rownames(alpha.div_gen))) {
    ps@sam_data$shannon_asv <- alpha.div_gen$observed
    ps@sam_data$richness_asv <- alpha.div_gen$diversity_shannon
  } else {stop(" Names don't match")}
  
  saveRDS(ps, file = "/fast/AG_Forslund/rob/studies/16S/gehrig_2019/dada2/phyloseq_gehrig_2019.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/studies/16S/gehrig_2019/dada2/phyloseq_gehrig_2019.rds")
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
  tsne_bray <- Rtsne(dist.bray, is_distance = TRUE, perplexity = 17)
  tsne_unifrac <- Rtsne(dist.unifrac, is_distance = TRUE, perplexity = 17)
  umap <- umap(otu_table(ps.prop))
  save(pcoa_bray, pcoa_unifr, tsne_bray, tsne_unifrac, umap, file = "dada2/ordinations.RData")
} else {
  load("dada2/ordinations.RData")
}

rownames(metadata_gehr) <- metadata_gehr$run_accession
cor(pcoa_bray$vectors,
    metadata_gehr[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head
cor(pcoa_unifr$vectors,
    metadata_gehr[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head


ps.prop@sam_data$age_fact <- as.factor(metadata_gehr$age)
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

print("db-NA comparisons: SILVA vs. GTDB")
print(taxdb_table)
print("time comparisons: SILVA vs. GTDB")
print(time_silva)
print(time_gtdb)
