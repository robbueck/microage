rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/studies/16S/hickman_2024"
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
library("gridExtra")
library("RColorBrewer")


# load additional custom functions
source("/fast/AG_Forslund/rob/mm_index/R_scripts/dada_2_functions.R")

# switch parts of the pipeline on/off
filter_step <- F
error_and_asv_step <- T
remove_chimeras_step <- T
id_taxa_step <- T
tree_step <- T
phyloseq_step <- T
diversity_step <- T
ordination_step <- T


# load metadata
metadata_hick <- read.table("metadata_hickman_2024_healthy.csv", sep = ",", header = T)

# remove samples with too high adapter content:
fastqc_forward <- read.table("/fast/AG_Forslund/rob/studies/16S/hickman_2024/fastqc/decontam/forward/multiqc_data/multiqc_fastqc.txt",
                             header = T, sep = "\t") %>%
  mutate(sample_ID = gsub("_dec_1", "", Sample))
fastqc_reverse <- read.table("/fast/AG_Forslund/rob/studies/16S/hickman_2024/fastqc/decontam/reverse/multiqc_data/multiqc_fastqc.txt",
                             header = T, sep = "\t") %>%
  mutate(sample_ID = gsub("_dec_2", "", Sample))
fastqc_merge <- left_join(fastqc_forward, fastqc_reverse, by = "sample_ID")
good_files <- fastqc_merge %>%
  filter(adapter_content.x == "pass", adapter_content.y == "pass") %>%
  pull(sample_ID)
metadata_hick <- metadata_hick %>%
  filter(run_accession %in% good_files) %>%
  mutate(group = ifelse(Instrument == "Illumina MiSeq", 
         yes = sample(c("A", "B"), 3995, replace = T),
         no = sample(c("C", "D", "E", "F", "G", "H", "I", "J", "K", "L"), 3995, replace = T)))
library_names <- unique(metadata_hick$group)
seqtab <- vector("list", length(library_names))
names(seqtab) <- library_names


# filter, learn errors and denoise
print(paste("analyzing", length(metadata_hick$run_accession), sep = " "))
print("Filter")
if(filter_step) {
  out <- filtering_fastq(sample_names = metadata_hick$run_accession,
                         trimLeft = 40,
                         truncLen=245,
                         maxN=0,
                         maxEE=c(2,2),
                         truncQ=2,
                         cores = 20)
  saveRDS(out, "dada2/hickman_out_table.rds")
} else {
  out <- readRDS("dada2/hickman_out_table.rds")
}


for(library in library_names) {
  print("error estimation and asv step")
  library_samples <- metadata_hick$run_accession[metadata_hick$group == library]
  print(paste("analyzing", length(library_samples), "samples in",
              library, which(library_names == library), "of",
              length(library_names), "libraries",
              sep = " "))
  out_file <- paste0(maindir, "/dada2/lib_", library, "_seqtab.rds")
  if(error_and_asv_step & !file.exists(out_file)){
    total_bases <- sum(out[metadata_hick$group == library,2]) * 205 # because error estimation is done for forward and reverse
    seqtab[[library]] <- error_and_asv(sample_names = library_samples, libname = library, nbases = 0.25*total_bases)
    saveRDS(seqtab[[library]], out_file) # save for each run separately
  } else {
    seqtab[[library]] <- readRDS(out_file)
  }
}
 # checked until here, need to adapt script further from here on
st.all <- mergeSequenceTables(tables = seqtab)

# how much got filtered out in each library?
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


# assess filter criteria for samples from each library
library_samples <- metadata_hick$run_accession[metadata_hick$cohort == "karelia"]
fnFs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(library_samples, "_dec_1.fastq.gz")))
fnRs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(library_samples, "_dec_2.fastq.gz")))
filtFs <- sort(file.path(maindir, "fastq_files/filtered", paste0(library_samples, "_dec_filt_1.fastq.gz")))
filtRs <- sort(file.path(maindir, "fastq_files/filtered", paste0(library_samples, "_dec_filt_2.fastq.gz")))
plotQualityProfile(c(fnFs[1], filtFs[1], fnRs[1], filtRs[1]))
ggsave("dada2/quality_profile.pdf", device = "pdf")

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
if(remove_chimeras_step){
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

# cdf_data_step_loss <- ddply(cdf_data_step_loss, .(step), transform, ecd=ecdf(lost_reads)(lost_reads))
ggplot(cdf_data_step_loss, aes(x=lost_reads, col = step)) + 
  stat_ecdf() +
  ylab("samples")
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

start.time_gtdb <- Sys.time()
if(id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/GTDB_r207-mod_April2022/GTDB_r207-mod_April2022.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids_gtdb <- IdTaxa(dna, trainingSet, processors=10, verbose=FALSE)
  saveRDS(ids_gtdb, "dada2/taxa_ids_GTDB.rds")
} else {
  ids_gtdb <- readRDS("dada2/taxa_ids_GTDB.rds")
}
end.time_gtdb <- Sys.time()
time_gtdb <- end.time_gtdb - start.time_gtdb

taxa_gtdb <- taxa_to_matrix(ids_gtdb)

taxdb_table <- t(data.frame(silva = (!is.na(taxa)) %>% colSums() / nrow(taxa),
                            gtdb = (!is.na(taxa_gtdb)) %>% colSums() / nrow(taxa_gtdb)))



print("tree step")
if (tree_step) {
  tree <- build_tree(st.all.nochim)
  saveRDS(tree, "dada2/tree.rds")
} else {
  tree <- readRDS("dada2/tree.rds")
}

# compare compositional plots between databases:
#add metadata
# rownames(metadata_hick) <- metadata_hick$run_accession
# metadata_hick <- metadata_hick[match(rownames(st.all.nochim), rownames(metadata_hick)),]
# level <- "phylum"

# SILVA
# ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE),
#                sample_data(metadata_hick),
#                tax_table(taxa), tree)
# dna <- Biostrings::DNAStringSet(taxa_names(ps))
# names(dna) <- taxa_names(ps)
# ps <- merge_phyloseq(ps, dna)
# taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
# 
# ps <- aggregate_taxa(ps, level = level)
# sample_names(ps) <- paste0(sample_names(ps), "_SILVA")
# sample_data(ps)$annot_DB <- "SILVA"

# #GTDB
# ps_GTDB <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
#                     sample_data(metadata_hick), 
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
#   tidyr::separate(X1, c("Baby_A", "DB_A"))%>%
#   tidyr::separate(X2, c("Baby_B", "DB_B"))%>%
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
#   add_xy_position(x = "Same_Individual")-> stats.test 
# 
# BC.Inf%>%
#   wilcox_test(dist ~ Same_DB)%>%
#   add_significance()%>%
#   add_xy_position(x = "Same_DB")-> stats.test_DB
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
library(ggplot2)
theme_set(theme_bw())
if(phyloseq_step) {
  metadata_hick <- read.table("metadata_hick_2018_healthy.csv", sep = ",", header = T)
  #add metadata
  head(metadata_hick)
  rownames(metadata_hick) <- metadata_hick$run_accession
  all(rownames(metadata_hick) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata_hick))
  metadata_hick <- metadata_hick[match(rownames(st.all.nochim), rownames(metadata_hick)),]
  all(rownames(metadata_hick) == rownames(st.all.nochim))
  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata_hick), 
                 tax_table(taxa), tree)
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  # save object
  saveRDS(ps, file = "/fast/AG_Forslund/rob/studies/16S/hickman_2024/dada2/phyloseq_hickman_2024.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/studies/16S/hickman_2024/dada2/phyloseq_hickman_2024.rds")
}

if (diversity_step) {# compare a-diversity between different groups
diversity_plots(ps_object = ps, x = "cohort", colour = "age")
ggsave("dada2/a_div_cohort.pdf", device = "pdf")
diversity_plots(ps_object = ps, x = "age", colour = "age")
ggsave("dada2/a_div_age.pdf", device = "pdf")
}

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
  save(dist.bray, dist.unifrac, pcoa_bray, file = "dada2/ordinations.RData")
} else {
  load("dada2/ordinations.RData")
}

rownames(metadata_hick) <- metadata_hick$run_accession
cor(pcoa_bray$vectors,
    metadata_hick[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head
cor(pcoa_unifr$vectors,
    metadata_hick[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head



#PCoA
# bray-curtis
p1 <- plot_ordination(ps.prop, pcoa_bray, color="age", title="Bray NMDS") +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.3) +
  scale_color_gradientn(name = "age", colors = topo.colors(7)) 
p2 <- plot_ordination(ps.prop, pcoa_bray, color="cohort", title="Bray NMDS") +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.3)
# unifrac dist:
p3 <- plot_ordination(ps.prop, pcoa_unifr, color="age", title="Unifrac NMDS") +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.3) +
  scale_color_gradientn(name = "age", colors = topo.colors(7))
p4 <- plot_ordination(ps.prop, pcoa_unifr, color="cohort", title="Unifrac NMDS") +
  geom_point(color = "white") +
  geom_point(alpha = 0.5, size = 0.3)
pf <- grid.arrange(p1, p2, p3, p4, ncol = 2)
ggsave("dada2/PCoA_age.pdf", device = "pdf", plot = pf)


fit_adonis <- adonis2(dist.bray ~ age + exclusive_breastfeeding + food + birthmode + sex,
                      data = metadata_hick, by="margin", na.action = na.omit)
saveRDS(fit_adonis,
        file = paste0("dada2/", "permanova.RData"))
fit_adonis <- readRDS(paste0("dada2/", "permanova.RData"))
print(fit_adonis)


libraries <- unique(metadata_hick$cohort)
lapply(libraries, function(x) mean(metadata_hick$age[metadata_hick$cohort == x]))
lapply(libraries, function(x) sum(metadata_hick$cohort == x))

lib_age <- metadata_hick[,c("age", "cohort")]
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
res_anova <- aov(age ~ cohort, data = metadata_hick)
summary(res_anova)

# check age and cohort influence by PERMANOVA
library(vegan)
res_per <- adonis2(dist.bray ~ age*cohort, data = data.frame(sample_data(ps)), permutations = 999)
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
