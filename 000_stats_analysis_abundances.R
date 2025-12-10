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

# source("000_stats_analysis.R")

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

# empty var
full_table <- vector()

for (a_dir in c("final_analysis_results_gITS7_ITS4", "final_analysis_results_ITS1f_ITS2", "final_analysis_results_Pauvert", "final_analysis_results_V6_V8")) {

#	a_dir <- "final_analysis_results_Pauvert" ; a_dir <- "final_analysis_results_V6_V8" ; a_dir <- "final_analysis_results_gITS7_ITS4"

	print(paste0("Analysing ", a_dir))

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

	# load physeq_combined
	physeq_combined <- readRDS(paste0(root_path, a_dir, "/physeq_combined.Rds"))
	# create table for plotting
	abundance <- colSums(otu_table(physeq_combined))
	if (a_dir == "final_analysis_results_Pauvert") {
		lib_size <- "large"

		# get treatment
		treatment <- unlist(lapply(sapply(strsplit(colnames(otu_table(physeq_combined)), "_"), "[", -1), paste, collapse="_"))

		# export taxa table
		write.table(tax_table(physeq_combined), paste0(save_abundances, a_dir, "_tax_table.csv"), quote=F, sep="\t", row.names=F)
	} else {
		# get lib size
		lib_size <- sapply(strsplit(colnames(otu_table(physeq_combined)), "_"), "[", 2)

		# get treatment
		treatment <- unlist(lapply(sapply(strsplit(colnames(otu_table(physeq_combined)), "_"), "[", -c(1:4)), paste, collapse="_"))
	}
	# get marker
	marker <- paste0(strsplit(a_dir, "_")[[1]][c(3,4)], collapse="_")
	# merge with the rest of the experiments
	full_table <- rbind.data.frame(full_table, cbind.data.frame(abundance, lib_size, treatment, marker, row.names=NULL))
}

##########################################################################################################################################################################

print("Library sizes from Pauvert")

##########################################################################################################################################################################

# ready an empty list
ab_plot <- list()

# loop through lib_sizes and markers
for (lib_size in unique(full_table$lib_size)) {
	# loop through markers
	for (a_marker in unique(full_table$marker)) {

		# lib_size <- "large" ; a_marker <- "results_Pauvert"

		if (a_marker == "results_Pauvert") {
			# subset table
			subset_table <- full_table[full_table$marker == a_marker & full_table$treatment!="T", ]
			# make title
			plot_title <- paste0("Library size Large ", marker)
			# make plots names
			lbl <- paste0("large_", a_marker)
		} else {
			# subset table
			subset_table <- full_table[full_table$lib_size == lib_size & full_table$marker == a_marker & full_table$treatment!="T", ]
			# make title
			plot_title <- paste0("Library size ", lib_size, " ", a_marker)
			# make plots names
			lbl <- paste0(lib_size, "_", a_marker)
		}

		# set seed for stats
		set.seed(123)
		# compute stats
		current_aov <- aov(abundance ~ treatment, data=subset_table)
		# summary to export info
		current_aov_tukey <- TukeyHSD(current_aov, ordered=TRUE)

		# get a table for tukey results to plot multiple comparisons
		current_tukey_table <- as.data.frame(current_aov_tukey$treatment)
		colnames(current_tukey_table) <- c("diff", "lwr", "uppr", "padj")

		# get adjusted pvals
		current_res <- current_tukey_table$padj

		# get rownames
		names(current_res) <- rownames(current_tukey_table)
		# compute stats letters
		current_letters <- multcompLetters(current_res)
		# put things together
		current_letters_df <- data.frame("treatment"=names(current_letters$Letters), "Letter"=as.vector(current_letters$Letters))
		# factor for plotting
		current_letters_df$treatment <- gsub("_", " ", current_letters_df$treatment)
		current_letters_df$treatment <- factor(current_letters_df$treatment, levels=unique(current_letters_df$treatment))
		if (a_marker == "results_Pauvert") {
			current_letters_df$max_y <- 49000
		} else {
			current_letters_df$max_y <- unique(ifelse(subset_table$lib_size=="small", 6000, ifelse(subset_table$lib_size=="medium", 13000, 27000)))
		}
		# remove unwanted
		current_letters_df <- current_letters_df[current_letters_df$treatment != "T", ]

		# set factors
		subset_table$treatment <- gsub("_", " ", subset_table$treatment)
		subset_table$treatment <- factor(subset_table$treatment, levels=unique(subset_table$treatment))

		# merge letters with subset_table
		merged_letters_table <- merge(subset_table, current_letters_df, by="treatment")

		# define palette
		treatment_palette <- brewer.pal(7, "Dark2")
		names(treatment_palette) <- levels(merged_letters_table$treatment)

		# plot
		ab_plot[[lbl]] <- ggplot(merged_letters_table, aes(x=treatment, y=abundance, fill=treatment)) +
					geom_boxplot(size=0.5, lwd=1, outlier.size = 10) +
					scale_fill_manual(values=treatment_palette) +
					guides(fill=guide_legend(title="Treatment", override.aes=list(shape=21, size=20), nrow=1)) +
					geom_jitter(colour="black", size=1, alpha=0.9) +
					geom_text(aes(x=treatment, y=.data$max_y, label=Letter), inherit.aes=F, size=20) +
					ggtitle(NULL) +
					ggplot_theme(leg_pos="bottom", fnt_sz=60) +
					ylab(NULL) +
					xlab(NULL) +
					theme(axis.text.x=element_blank())
	}
}

# get legend
legend <- get_legend(ab_plot[["small_results_gITS7"]])
# put plots together
ab_plot_merge <- cowplot::plot_grid(
			cowplot::plot_grid(ab_plot[["small_results_ITS1f"]] + theme(legend.position="none") + ggtitle("ITS1") + ylab("Small") + ylim(-1, 6500),
						ab_plot[["small_results_gITS7"]] + theme(axis.text.y=element_blank(), legend.position="none") + ggtitle("ITS2") + ylim(-1, 6500),
						ab_plot[["small_results_V6"]] + theme(axis.text.y=element_blank(), legend.position="none") + ggtitle("18S") + ylim(-1, 6500),
						ncol=3),
			cowplot::plot_grid(ab_plot[["medium_results_ITS1f"]] + theme(legend.position="none") + ylab("Medium") + ylim(-1, 14500),
						ab_plot[["medium_results_gITS7"]] + theme(axis.text.y=element_blank(), legend.position="none") + ylim(-1, 14500),
						ab_plot[["medium_results_V6"]] + theme(axis.text.y=element_blank(), legend.position="none") + ylim(-1, 14500),
						ncol=3),
			cowplot::plot_grid(ab_plot[["large_results_ITS1f"]] + theme(legend.position="none") + ylab("Large") + ylim(-1, 29000),
						ab_plot[["large_results_gITS7"]]  + theme(axis.text.y=element_blank(), legend.position="none") + ylim(-1, 29000),
						ab_plot[["large_results_V6"]] + theme(axis.text.y=element_blank(), legend.position="none") + ylim(-1, 29000),
						ncol=3), legend, nrow=4, rel_heights=c(1,1,1,0.3))

# export plot
export_final_svg(paste0(supplementary_figs, "fig_s1"), ab_plot_merge, base=F, width=168*10, height=168*8)

# pauvert
ab_plot_merge_pauvert <- ab_plot[["large_results_Pauvert"]] + facet_wrap(~treatment, nrow=1, scale="free_x")

# export plot
export_final_svg(paste0(supplementary_figs, "fig_s2"), ab_plot_merge_pauvert, base=F, width=168*8, height=168*4)

