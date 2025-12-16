rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/studies/16S/bockulich_2016/"
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
ordination_step <- F


# load metadata
metadata_bockulich <- read.table("metadata_boc_2016_healthy.csv", sep = ",", header = T)
library_names <- unique(metadata_bockulich$run_prefix_prep)
if(any(is.na(library_names))){
  stop("Library names contain NAs")
}
out <- vector("list", length(library_names))
seqtab <- vector("list", length(library_names))
names(out) <- library_names
names(seqtab) <- library_names

# filter, learn errors and denoise seperately
for(library in library_names) {
  library_samples <- metadata_bockulich$sample_ID[metadata_bockulich$run_prefix_prep == library]
  print(paste("analyzing", length(library_samples), "samples in",
              library, which(library_names == library), "of",
              length(library_names), "libraries",
              sep = " "))
  print("Filter")
  if(filter_step) {
    out[[library]] <- filtering_fastq(sample_names = library_samples,
                                             trimLeft = 10,
                                             truncLen=135,  # for adapter trimming
                                             maxN=0,
                                             maxEE=2,
                                             truncQ=2)
    saveRDS(out[[library]], paste0(maindir, "/dada2/lib_", library, "_out_table.rds")) # save for each run separately
  } else {
    out[[library]] <- readRDS(paste0(maindir, "/dada2/lib_", library, "_out_table.rds"))
  }
  if(error_and_asv_step) {
    total_bases <- sum(out[[library]][,2]) * 125 # because error estimation is done for forward and reverse
    seqtab[[library]] <- error_and_asv(sample_names = library_samples, libname = library, nbases = 0.25 * total_bases)
    saveRDS(seqtab[[library]], paste0(maindir, "/dada2/lib_", library, "_seqtab.rds")) # save for each run separately
  } else {
    seqtab[[library]] <- readRDS(paste0(maindir, "/dada2/lib_", library, "_seqtab.rds"))
  }
}

st.all <- mergeSequenceTables(tables = seqtab)


# how much got filtered out in each library?
out[["combined"]] <-  do.call("rbind", out)
cdf_data_filter <- melt(lapply(out, function(x){loss_ratio(x[,2], x[,1])}),
                 value.name = "filter_ratio")
colnames(cdf_data_filter)[2] <- "Library"
# violin plot:
sample_size = cdf_data_filter %>% group_by(Library) %>% dplyr::summarize(num=n())
cdf_data_filter <- cdf_data_filter %>% 
  left_join(sample_size) %>%
  mutate(Library = paste0(Library, "\n", "n=", num))
ggplot(cdf_data_filter, aes(x=Library, y = filter_ratio, color = Library)) + 
  geom_violin() +
  ylim(0,1) +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/violin_filter_ratios_per_library.pdf", device = "pdf")
# cdf_data_filter <- ddply(cdf_data_filter, .(Library), transform, ecd=ecdf(filter_ratio)(filter_ratio))
ggplot(cdf_data_filter, aes(x=filter_ratio, col = Library)) + 
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
ggsave("dada2/cfd_filter_ratios_per_library.pdf", device = "pdf")
# some libraries seem to be quite bad, check qualities along seqs


# 
# 
# # assess filter criteria for samples from each library
# library_samples <- metadata_bockulich$sample_ID[metadata_bockulich$run_prefix_prep == "SeqRun1"]
# fnFs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(library_samples, "_dec_1.fastq.gz")))
# fnRs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(library_samples, "_dec_2.fastq.gz")))
# filtFs <- sort(file.path(maindir, "fastq_files/filtered", paste0(library_samples, "_dec_filt_1.fastq.gz")))
# filtRs <- sort(file.path(maindir, "fastq_files/filtered", paste0(library_samples, "_dec_filt_2.fastq.gz")))
# plotQualityProfile(c(fnFs[1], filtFs[1], fnRs[1], filtRs[1]))
# ggsave("dada2/quality_profile.pdf", device = "pdf")


# Inspect distribution of sequence lengths
cdf_data_seqlen <- melt(lapply(seqtab, function(x){nchar(getSequences(x))}),
                 value.name = "seq_length")
colnames(cdf_data_seqlen)[2] <- "Library"
# violin plot:
sample_size = cdf_data_seqlen %>% group_by(Library) %>% dplyr::summarize(num=n())
cdf_data_seqlen <- cdf_data_seqlen %>% 
  left_join(sample_size) %>%
  mutate(Library = paste0(Library, "\n", "n=", num))
ggplot(cdf_data_seqlen, aes(x=Library, y = seq_length, color = Library)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/violin_seq_lengths_per_library.pdf", device = "pdf")

# cdf_data_seqlen <- ddply(cdf_data_seqlen, .(Library), transform, ecd=ecdf(seq_length)(seq_length))
ggplot(cdf_data_seqlen, aes(x=seq_length, col = Library)) + 
  stat_ecdf() +
  ylab("samples") 
ggsave("dada2/cfd_seq_lengths_per_library.pdf", device = "pdf")



# remove chimeras
print("remove chimeras")
if (remove_chimeras_step) {
  st.all.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=10, verbose=TRUE)
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
rownames(out$combined) <- gsub("_dec_1.fastq.gz", "", rownames(out$combined))
loss_table <- merge(out$combined, read_counts_raw, by = 0)
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
colnames(cdf_data_step_loss)[1:2] <- c("sample", "step")
# violin plot:
ggplot(cdf_data_step_loss, aes(x=step, y = lost_reads, color = step)) + 
  geom_violin() +
  ylim(0,1) +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/violin_filter_ratios_per_step.pdf", device = "pdf")

# cdf plot
# cdf_data_step_loss <- ddply(cdf_data_step_loss, .(step), transform, ecd=ecdf(lost_reads)(lost_reads))
ggplot(cdf_data_step_loss, aes(x=lost_reads, col = step)) + 
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
cdf_data_step_loss[is.na(cdf_data_step_loss$lost_reads),]
loss_table[loss_table$Row.names %in% cdf_data_step_loss$sample[is.na(cdf_data_step_loss$lost_reads)],]
ggsave("dada2/cfd_filter_ratios_per_step.pdf", device = "pdf")

# save loss table to be used for further stats
loss_table <- dplyr::transmute(loss_table, Row.names = Row.names, read_counts_raw = read_counts_raw,
                               decontaminated = reads.in, filtered = reads.out,
                               denoised = denoised, removed_chimeras = removed_chimeras)
write.csv(loss_table, "dada2/read_loss_per_step.csv")


# assign taxonomy
library(DECIPHER)
dna <- DNAStringSet(getSequences(st.all.nochim)) # Create a DNAStringSet from the ASVs
print("id taxa")
start.time_silva <- Sys.time()
if (id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/silva/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids <- IdTaxa(dna, trainingSet, processors=10, verbose=FALSE)
  saveRDS(ids, "dada2/taxa_ids.rds")
} else {
  ids <- readRDS("dada2/taxa_ids.rds")
}
end.time_silva <- Sys.time()
time_silva <- end.time_silva - start.time_silva

taxa <- taxa_to_matrix(ids)
# 
# print("GTDB")
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
#                             gtdb = (!is.na(taxa_gtdb)) %>% colSums() / nrow(taxa_gtdb)))


if (tree_step) {
  tree <- build_tree(st.all.nochim)
  saveRDS(tree, "dada2/tree.rds")
} else {
  tree <- readRDS("dada2/tree.rds")
}

# 
# # compare compositional plots between databases:
# #add metadata
# rownames(metadata_bockulich) <- metadata_bockulich$sample_ID
# metadata_bockulich <- metadata_bockulich[match(rownames(st.all.nochim), rownames(metadata_bockulich)),]
# level <- "family"
# 
# # SILVA
# ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
#                sample_data(metadata_bockulich), 
#                tax_table(taxa), tree)
# dna <- Biostrings::DNAStringSet(taxa_names(ps))
# names(dna) <- taxa_names(ps)
# ps <- merge_phyloseq(ps, dna)
# taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
# sample_names(ps) <- paste0(sample_names(ps), "_SILVA")
# sample_data(ps)$annot_DB <- "SILVA"
# ps <- ps %>%
#   aggregate_taxa(level = level) %>%
#   microbiome::transform(transform = "compositional") %>%
#   filter_taxa(function(x) mean(x) > 1e-5, TRUE) %>%
#   filter_taxa(function(x){sum(x > 0) > (0.05*length(x))}, TRUE) #%>%
#   #aggregate_taxa(level = level)
# 
# 
# #GTDB
# ps_GTDB <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
#                     sample_data(metadata_bockulich), 
#                     tax_table(taxa_gtdb), tree) 
# dna <- Biostrings::DNAStringSet(taxa_names(ps_GTDB))
# names(dna) <- taxa_names(ps_GTDB)
# ps_GTDB <- merge_phyloseq(ps_GTDB, dna)
# taxa_names(ps_GTDB) <- paste0("ASV", seq(ntaxa(ps_GTDB)))
# sample_names(ps_GTDB) <- paste0(sample_names(ps_GTDB), "_GTDB")
# sample_data(ps_GTDB)$annot_DB <- "GTDB"
# ps_GTDB <- ps_GTDB %>% 
#   aggregate_taxa(level = level) %>%
#   microbiome::transform(transform = "compositional") %>%
#   filter_taxa(function(x) mean(x) > 1e-5, TRUE) %>%
#   filter_taxa(function(x){sum(x > 0) > (0.05*length(x))}, TRUE) #%>%
#   
# 
# 
# taxdb_table <- t(data.frame(silva = (tax_table(ps) != "Unknown") %>% colSums() / ntaxa(ps),
#                             gtdb = (tax_table(ps_GTDB) != "Unknown") %>% colSums() / ntaxa(ps_GTDB)))
# taxdb_table <- t(data.frame(silva = (!is.na(tax_table(ps))) %>% colSums() / ntaxa(ps),
#                             gtdb = (!is.na(tax_table(ps_GTDB))) %>% colSums() / ntaxa(ps_GTDB)))
# 
# 
# 
# # merge annotations for plotting
# ps_combined <- merge_phyloseq(ps_GTDB, ps) %>%
#   subset_samples(age > 100) %>%
#   subset_taxa(get(level) != "Unknown") %>%
#   # microbiome::transform(transform = "compositional") %>%
#   filter_taxa(function(x) mean(x) > 1e-5, TRUE) %>%
#   filter_taxa(function(x){sum(x > 0) > (0.05*length(x))}, TRUE)
# 
# n_taxa <- ps_combined@tax_table[,level] %>% unique() %>% length()
# plot_composition(ps_combined,
#                  average_by = "annot_DB", 
#                  transform = "compositional") +
#   scale_fill_manual(level, values = colorRampPalette(brewer.pal(12, "Paired"))(n_taxa)) +
#   guides(fill=guide_legend(ncol=2)) +
#   theme(axis.text.x = element_text(angle=90, hjust=1),
#         legend.text = element_text(face = "italic", size = 7),
#         legend.title = element_text(size=10),
#         legend.key.size = unit(0.2, 'cm'))
# 
# 
# # compare b-div between the two dbs
# dist.bray <- phyloseq::distance(ps_combined, method = "bray")
# pcoa_bray <- pcoa(dist.bray)
# 
# plot_ordination(ps_combined, pcoa_bray, color="age", title="Bray NMDS", shape = "annot_DB") +
#   scale_color_gradientn(name = "age", colors = topo.colors(7))
# 
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
#   separate(X1, c("Baby_A", "DB_A"), sep = "[_]")%>%
#   separate(X2, c("Baby_B", "DB_B"), sep = "[_]")%>%
#   dplyr::mutate(Same_Individual = case_when(Baby_A == Baby_B  ~ T,
#                                             Baby_A != Baby_B ~ F))%>%
#   dplyr::mutate(Same_DB = case_when(DB_A == DB_B  ~ T,
#                                     DB_A != DB_B ~ F))%>%
#   dplyr::select(c("Sample_pair", "Baby_A", "Baby_B","Same_Individual", "Same_DB", "dist"))-> BC.Inf
# 
# library(rstatix)
# BC.Inf%>%
#   wilcox_test(dist ~ Same_Individual)%>%
#   add_significance()%>%
#   add_xy_position(x = "Same_Individual")
# 
# BC.Inf%>%
#   wilcox_test(dist ~ Same_DB)%>%
#   add_significance()%>%
#   add_xy_position(x = "Same_DB")
# 
# BC.Inf%>%
#   ggplot(aes(x= Same_DB, y= dist, fill= Same_DB))+
#   geom_boxplot(aes(),outlier.shape=NA)+
#   geom_violin()+
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
#   geom_violin()+
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

if (phyloseq_step) {
  #add metadata
  metadata_bockulich <- read.table("metadata_boc_2016_healthy.csv", sep = ",", header = T)
  rownames(metadata_bockulich) <- metadata_bockulich$sample_ID
  rownames(st.all.nochim) <- paste0(rownames(st.all.nochim), "_bockulich_2016")
  all(rownames(metadata_bockulich) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata_bockulich))
  metadata_bockulich <- metadata_bockulich[match(rownames(st.all.nochim), rownames(metadata_bockulich)),]
  all(rownames(metadata_bockulich) == rownames(st.all.nochim))
  # number of ASVs unassigned in whole study
  metadata_bockulich$n_unassigned_asvs_study <-   sum(grepl("unclassified", taxa[,6]))
  metadata_bockulich$n_total_asvs_study <- nrow(taxa)
  # number of ASVs not assigned to genus level per sample
  matching_taxa <- rownames(taxa)[grepl("unclassified", taxa[,6])]
  unassigned_taxa_present_per_sample <- apply(st.all.nochim[,matching_taxa], 1, function(x) sum(x > 0))
  total_taxa_present_per_sample <- apply(st.all.nochim, 1, function(x) sum(x > 0))
  metadata$n_unassigned_asvs_sample <- unassigned_taxa_present_per_sample
  metadata$n_total_asvs_sample <- total_taxa_present_per_sample

  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata_bockulich), 
                 tax_table(taxa), tree)
  sample_names(ps) <- sample_data(ps)$run_accession
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
  
  # save object
  saveRDS(ps, file = "/fast/AG_Forslund/rob/studies/16S/bockulich_2016/dada2/phyloseq_bockulich_2014.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/studies/16S/bockulich_2016/dada2/phyloseq_bockulich_2014.rds")
}

# compare a-diversity between different groups
diversity_plots(ps_object = ps, x = "run_prefix_prep", colour = "age")
ggsave("dada2/a_div_cohort.pdf", device = "pdf")
diversity_plots(ps_object = ps, x = "age", colour = "log(read_counts_raw)")
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
  tsne_bray <- Rtsne(dist.bray, is_distance = TRUE)
  tsne_unifrac <- Rtsne(dist.unifrac, is_distance = TRUE)
  umap <- umap(otu_table(ps.prop))
  save(dist.bray, dist.unifrac, pcoa_bray, pcoa_unifr, tsne_bray, tsne_unifrac, umap, dist.bray, dist.unifrac, file = "dada2/ordinations.RData")
} else {
  load("dada2/ordinations.RData")
}

rownames(metadata_bockulich) <- metadata_bockulich$sample_ID
cor(pcoa_bray$vectors,
    metadata_bockulich[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head
cor(pcoa_unifr$vectors,
    metadata_bockulich[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head

plot_ordination(ps.prop, pcoa_bray, title="Bray NMDS", color="run_prefix_prep")+
  geom_point(color = "white") +
  geom_point(alpha = 0.3, size = 0.5)


#PCoA
# bray-curtis
p1 <- plot_ordination(ps.prop, pcoa_bray, color="age", title="Bray NMDS") +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_gradientn(name = "age", colors = topo.colors(7))
p2 <- plot_ordination(ps.prop, pcoa_bray, color="run_prefix_prep", title="Bray NMDS")+
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5)
# unifrac dist:
p3 <- plot_ordination(ps.prop, pcoa_unifr, color="age", title="Unifrac NMDS") +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_gradientn(name = "age", colors = topo.colors(7))
p4 <- plot_ordination(ps.prop, pcoa_unifr, color="run_prefix_prep", title="Unifrac NMDS") +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5)
pf <- grid.arrange(p1, p2, p3, p4, ncol = 2)
ggsave("dada2/PCoA_age.pdf", device = "pdf", plot = pf)


fit_adonis <- adonis2(dist.bray ~ age + birthmode + sex,
                      data = metadata_bockulich, by="margin", na.action = na.omit)
saveRDS(fit_adonis,
        file = paste0("dada2/", "permanova.RData"))
fit_adonis <- readRDS(paste0("dada2/", "permanova.RData"))
print(fit_adonis)



# t-SNE
tsnedata <- data.frame(tsne_bray$Y) %>%
  dplyr::rename(tSNE1 = X1, tSNE2 = X2) %>%
  bind_cols(data.frame(sample_data(ps)))
print("A")
ps1 <- ggplot(tsnedata, aes(x = tSNE1, y = tSNE2, color = age)) +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  theme(aspect.ratio = 1) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) + 
  ggtitle("t-SNE Bray")
ps2 <- ggplot(tsnedata, aes(x = tSNE1, y = tSNE2, color = run_prefix_prep)) +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  theme(aspect.ratio = 1) +
  ggtitle("t-SNE Bray")
# unifrac
tsnedata <- data.frame(tsne_unifrac$Y) %>%
  dplyr::rename(tSNE1 = X1, tSNE2 = X2) %>%
  bind_cols(data.frame(sample_data(ps)))
ps3 <- ggplot(tsnedata, aes(x = tSNE1, y = tSNE2, color = age)) +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  theme(aspect.ratio = 1) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) +
  ggtitle("t-SNE Unifrac")
ps4 <- ggplot(tsnedata, aes(x = tSNE1, y = tSNE2, color = run_prefix_prep)) +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  theme(aspect.ratio = 1) +
  ggtitle("t-SNE Unifrac")
psf <- grid.arrange(ps1, ps2, ps3, ps4, ncol = 2)
ggsave("dada2/tSNE.pdf", device = "pdf", plot = psf)

# UMAP
umapdata <- umap$layout %>%
  data.frame() %>%
  dplyr::rename(UMAP1 = X1, UMAP2 = X2) %>%
  bind_cols(data.frame(sample_data(ps)))
pu1 <- ggplot(umapdata, aes(x = UMAP1, y = UMAP2, color = age)) +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  theme(aspect.ratio = 1) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) +
  ggtitle("UMAP")
pu2 <- ggplot(umapdata, aes(x = UMAP1, y = UMAP2, color = run_prefix_prep)) +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.5) +
  theme(aspect.ratio = 1) +
  ggtitle("UMAP")
puf <- grid.arrange(pu1, pu2, ncol = 2)
ggsave("dada2/umap.pdf", device = "pdf", plot = puf)


# are the age groups equally distributed among the libraries?
libraries <- unique(metadata_bockulich$run_prefix_prep)
lapply(libraries, function(x) mean(metadata_bockulich$age[metadata_bockulich$run_prefix_prep == x]))
lapply(libraries, function(x) sum(metadata_bockulich$run_prefix_prep == x))

lib_age <- metadata_bockulich[,c("age", "run_prefix_prep")]
colnames(lib_age)[2] <- "library"
lib_age_total <- data.frame(age = lib_age$age, library = "total")
lib_age <- rbind(lib_age, lib_age_total)
# violin plot:
sample_size = lib_age %>% group_by(library) %>% dplyr::summarize(num=n())
lib_age <- lib_age %>% 
  left_join(sample_size) %>%
  mutate(library = paste0(library, "\n", "n=", num))

# cdf_data_seqlen <- ddply(cdf_data_seqlen, .(Library), transform, ecd=ecdf(seq_length)(seq_length))
ggplot(lib_age, aes(x=age, col = library)) + 
  stat_ecdf() +
  ylab("samples") 

ggplot(data=lib_age, aes(x=age, group=library, fill=library)) +
  geom_density(adjust=1.5, alpha=.6)
ggsave("dada2/library_age_density.pdf")
# test age distribution across libraries with anova:
res_anova <- aov(age ~ run_prefix_prep, data = metadata_bockulich)
summary(res_anova)

# check age and cohort influence by PERMANOVA
library(vegan)
res_per <- adonis2(dist.bray ~ age*run_prefix_prep, data = data.frame(sample_data(ps)), permutations = 999)
res_per
summary(res_per)



# export to be used with picrust:
writeXStringSet(ps@refseq, "dada2/ASVs.fasta", format = "fasta")
otu_table_out <- t(as.data.frame(ps@otu_table))
write.table(data.frame(OTU = rownames(otu_table_out), otu_table_out), file = "dada2/otu_table.tsv", sep = "\t", row.names = F, quote = F)


print("db-NA comparisons: SILVA vs. GTDB")
print(taxdb_table)
print("time comparisons: SILVA vs. GTDB")
print(time_silva)
print(time_gtdb)

