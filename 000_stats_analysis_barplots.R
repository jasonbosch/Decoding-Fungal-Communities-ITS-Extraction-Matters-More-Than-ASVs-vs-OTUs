##########################################################################################################################################################################

# all the species available from ensembl fungi (available at https://fungi.ensembl.org/index.html) were downloaded. release 57 or version 110????

# for each species in the set, genomic regions corresponding to different primer pairs were extracted.

# for each primer pair, 10 different communities were randomly built, at various levels of richness (namely 50, 150, 250, 350, 450, 500, 600, 700, and 800). each level or richness comprise 10 communities, for a total of 90 samples for each primer pair.

# at this point, grinder was used to simulate the actual sequencing using three different library size that are: small (10k reads), medium (25k reads), and large (50k reads). at this point, 90 samples x 3 library size accounts for a total of 270 samples.

# generated reads length was set at 300bp, joining of forward and reverse was performed with 40bp overlap and 15% of mismatch allowed.

# true communities were built for each sample, i.e. the combination of primer pair, community, richness, and library size. to build the true community, all the species that were obtained as an output from grinder were considered as the true community. these were used, later on, as a benchmark to evaluate the performances of ASVs and OTUs.

# finally, the mock sequencing data were finally processed using either ASVs or OTUs, with or without ITS extraction.

# for each sample, alpha diversity and beta diversity were estimated. alpha diversity was estimated after performing rarefaction by means of repeated multiple iterations of rarefaction implemented by mirlyn. beta diversity was estimated using Shannon, Simpson, and Aitchison metrics. specificity and sensitivity were also computed, using the true community as a reference.

# finally, some statistics were computed to test for differences between ASVs and OTUs methodologies.

##########################################################################################################################################################################

# some info to remember. the construction of the community follows a "gradient":
# some features are set and some are a consequence. so, given:

# a fixed set of available organisms, depending on the database
# Diversity Levels: 50 < 150 < 250 < 350 < 450 < 500 < 600 < 700 < 800
# Library_Size Levels: Small < Medium < Large

# and the fact that communities were analysed using:

# Treatment Levels: FWD_ASV_Rolling FWD_ASV FWD_OTU ITS_ASV_Rolling ITS_ASV ITS_OTU

# leads to having a specific representation for each mock community.
# eack communiti is then represented by it's components at various levels:

# Level Levels: Kingdom < Phylum < Class < Order < Family < Genus

# for this reason it may be reasonable to apply some model to compare one level
# with another, to see how diversity, lib_size, and treatment affect the diversity.

# with respect to the other comparisons, it is probably sufficient to apply hypotesis
# tests because the aim would be to see if the differences between one another, e.g.
# between different library sizes, are significant.

##########################################################################################################################################################################
# Rscript 000_stats_analysis_barplots.R
# source("000_stats_analysis_barplots.R")

# load functions
source("/mnt/cinqueg/gabriele/work/microbiology/other_analyses/jason/000_stats_analysis_functions.R")

# load libraries
library("phyloseq")
library("RColorBrewer")
library("ggpubfigs")
library("MetBrewer")
library("ggplot2")
library("cowplot")
library("ggtext")
library("ggrepel")
library("dplyr")
library("multcomp")
library("FSA")
library("ggpubr")
library("multcompView")
library("viridis")
library("ggnewscale")

# set paths
root_path <- "/mnt/cinqueg/gabriele/work/microbiology/other_analyses/jason/data/"
save_means <- paste0(root_path, "tables_final/means/")
save_stats <- paste0(root_path, "tables_final/stats/")
save_abundances <- paste0(root_path, "tables_final/abundances/")
supplementary_figs <- paste0(root_path, "tables_final/supplementary_figs/")
ifelse(dir.exists(supplementary_figs), "all cool", dir.create(supplementary_figs))
ifelse(dir.exists(save_stats), "all cool", dir.create(save_stats))
ifelse(dir.exists(save_means), "all cool", dir.create(save_means))
ifelse(dir.exists(save_abundances), "all cool", dir.create(save_abundances))

##########################################################################################################################################################################

print("plotting barplots")

##########################################################################################################################################################################

print("ITS1 and 2")

# load true communities
true_its <- list()
true_its[["ITS1"]] <- readRDS(paste0(root_path, "final_analysis_results_ITS1f_ITS2/", "physeq_true.Rds"))
true_its[["ITS2"]] <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/", "physeq_true.Rds"))

# counter for the number of supplementary figs
count <- 3

# loop through levels
for (taxa_level in c("Phylum", "Class", "Order", "Family", "Genus")) {

#	taxa_level <- "Genus" ; marker <- "ITS1"

	print("Final plot")
	# barplots for supplementary
	barplot_supplementary(true_its, taxa_level, count)
	# increase counter
	count <- count + 1
}

# barplot main text, Order only
#barplot_main_phylum(true_its, "Order")

print("now, Pauvert")

# now pauvert
true_its[["Pauvert_combined"]] <- readRDS(paste0(root_path, "final_analysis_results_Pauvert/", "physeq_combined.Rds"))

# load pauvert mock community
pauvert_mock <- read.table("/mnt/cinqueg/gabriele/work/microbiology/other_analyses/jason/data/1-s2.0-S1754504818302800-mmc2.csv", header = T, sep=";")

# add in presumed abundance given equimolar pooling
pauvert_mock$Abundance <- 1/nrow(pauvert_mock)

# simplify
colnames(pauvert_mock)[1] <- "Strain"
pauvert_mock <- pauvert_mock[, c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain", "Abundance")]

# loop through levels
# update count
count <- 27
for (taxa_level in c("Phylum", "Class", "Order", "Family", "Genus")) {
	# loop through treatments
	all_abundance_treats <- as.data.frame(unique(tax_table(true_its[["Pauvert_combined"]])[, taxa_level]))
	colnames(all_abundance_treats) <- "taxa"
	for (a_treat in unique(sample_data(true_its[["Pauvert_combined"]])$Treatment)) {
#		taxa_level <- "Genus" ; marker <- "Pauvert_combined" ; a_treat <- "TRUE"

		print("Final plot")
		# only for true community. the rest we only need the abundance table
		a_plot <- ifelse(a_treat == "TRUE", TRUE, FALSE)

		# barplots for supplementary
		abundance_treat <- barplot_supplementary_pauvert(phylo_obj=true_its[["Pauvert_combined"]], taxa_level=taxa_level, treatment=a_treat, count=count, plotting=a_plot)
		colnames(abundance_treat) <- c("taxa", a_treat)
		all_abundance_treats <- merge(all_abundance_treats, abundance_treat, by="taxa", all=TRUE)
	}

	# pauvert mock
	pauv_mock_aggregated <- aggregate(as.formula(paste("Abundance ~", taxa_level, sep = " ")), pauvert_mock[, c("Abundance", taxa_level)], sum)
	colnames(pauv_mock_aggregated) <- c("taxa", "MOCK")

	# add mock to all
	all_abundance_treats_mock <- merge(pauv_mock_aggregated, all_abundance_treats, by="taxa", all=TRUE)

	# increase counter for each taxa_level
	count <- count + 1
	# finally export
	write.table(all_abundance_treats_mock, paste0(save_abundances, taxa_level, "_abundance_table_Pauvert.csv"), quote=F, sep="\t", row.names=F)
}

#graphics.off() ; print("done, ciao!")

