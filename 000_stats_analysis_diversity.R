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
# Rscript 000_stats_analysis_diversity.R
# source("000_stats_analysis_diversity.R")

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

# library FSA is used for dunnTest
# loop through ITS1, 2, Pauvert. ORDER IS IMPORTANT FOR EXPORTING WITH CORRECT NAMES
for (a_dir in c("final_analysis_results_ITS1f_ITS2", "final_analysis_results_gITS7_ITS4", "final_analysis_results_Pauvert")) {

#	a_dir <- "final_analysis_results_gITS7_ITS4" ; # a_dir <- "final_analysis_results_Pauvert"

	print(paste0("Analysing ", a_dir))

	# load data
	similarity_scores <- readRDS(paste0(root_path, a_dir, "/similarity_scores.Rds"))
	colnames(similarity_scores) <- c("Diversity", "Library_Size", "Replicate", "Processing", "Method", "Treatment", "Level", "Sensitivity", "Specificity", "Aitchison", "Hellinger", "BrayCurtis")

	# load combined_alphadiv_df
	combined_alphadiv_df <- readRDS(paste0(root_path, a_dir, "/combined_alphadiv_df.Rds"))

	# set levels
	combined_alphadiv_df$Treatment <- gsub("_", " ", combined_alphadiv_df$Treatment)
	combined_alphadiv_df$Treatment <- factor(combined_alphadiv_df$Treatment, levels=unique(combined_alphadiv_df$Treatment))

	# compute means
	combined_alphadiv_df_averaged <- combined_alphadiv_df %>%
						group_by(Library_Size, Diversity, Treatment, Replicate) %>%
						summarise(simpson_mean = mean(Simpson), shannon_mean = mean(Shannon), observed_mean = mean(Observed)) %>%
						data.frame()

	# as numeric
	combined_alphadiv_df_averaged$Diversity <- as.numeric(as.character(combined_alphadiv_df_averaged$Diversity))
	combined_alphadiv_df$Diversity <- as.numeric(as.character(combined_alphadiv_df$Diversity))

	# compute averages to add to plot
	combined_alphadiv_df %>%
		group_by(Treatment) %>%
		summarize(mean_observed=mean(Observed), mean_diversity=mean(Diversity)) %>%
		data.frame() %>%
		write.table(paste0(save_means, a_dir, "_observed_average.csv"), quote=F, sep="\t", row.names=F)

	# set factors
	combined_alphadiv_df_averaged$Treatment <- as.character(combined_alphadiv_df_averaged$Treatment)
	combined_alphadiv_df_averaged$Treatment <- factor(combined_alphadiv_df_averaged$Treatment, levels=unique(combined_alphadiv_df_averaged$Treatment))
	combined_alphadiv_df_averaged$Diversity <- factor(combined_alphadiv_df_averaged$Diversity, levels=unique(combined_alphadiv_df_averaged$Diversity))

	##########################################################################################################################################################################

	print("plot diversities")

	##########################################################################################################################################################################

	if (a_dir == "final_analysis_results_Pauvert") {

		print(a_dir)

		# subset alphadiv large
		lib_size_large <- combined_alphadiv_df_averaged[combined_alphadiv_df_averaged$Library_Size=="Large", ]

		# define palette
		treatment_palette <- brewer.pal(7, "Dark2")
		names(treatment_palette) <- unique(lib_size_large$Treatment)

		# factor
		lib_size_large$Treatment <- gsub("_", " ", lib_size_large$Treatment)
		similarity_scores$Treatment <- gsub("_", " ", similarity_scores$Treatment)
		lib_size_large$Treatment <- factor(lib_size_large$Treatment, levels=unique(lib_size_large$Treatment))
		similarity_scores$Treatment <- factor(similarity_scores$Treatment, levels=unique(similarity_scores$Treatment))

		# compute observed differences
		lib_size_large$Observed_Diff <- as.numeric(lib_size_large$observed_mean)-as.numeric(as.character(lib_size_large$Diversity))

		# Draw richness plots
		alpha_observed <- ggplot(lib_size_large, aes(x=Treatment, y=Observed_Diff, fill=Treatment)) +
					geom_point(colour="black", size=10, shape=21) +
					scale_fill_manual(values=treatment_palette) +
					ylab("Difference in richness") +
					ggtitle(NULL) +
					guides(fill=guide_legend(override.aes=list(shape=21, size=10))) +
					facet_wrap(~Treatment, nrow=1, scale="free_x") +
					ggplot_theme(ang_le=45, leg_pos="none", fnt_sz=50) +
					theme(axis.text.x=element_blank()) +
					xlab(NULL)

		# export
		export_final_svg(paste0(supplementary_figs, "fig_S33"), alpha_observed, base=F, width=168*8, height=168*4)
	} else {

		print(a_dir)

################		# plot
################		diversity_plot <- ggplot(similarity_scores[similarity_scores$Level!="Kingdom" & similarity_scores$Diversity==500, ], aes(x=Specificity, y=Sensitivity, fill=Treatment, shape=Level)) +
################					geom_point(colour="black", size=20, stroke=1) +
################					scale_shape_manual(values=c(21:25)) +
################					guides(fill=guide_legend(override.aes=list(shape=21, size=30)), shape=guide_legend(override.aes=list(size=30))) +
################					ggtitle(NULL) +
################					facet_wrap(~Treatment) +
################					ggplot_theme(fnt_sz=50)
################		# export
################		export_final_svg(paste0(supplementary_figs, a_dir, "_diversity"), diversity_plot, base=F, width=168*6, height=168*4)

		# subset alphadiv large and diversity equal to 500
		lib_size_large <- combined_alphadiv_df_averaged[combined_alphadiv_df_averaged$Library_Size=="Large", ]

		# define palette
		diversity_palette <- brewer.pal(9, "Paired")
		names(diversity_palette) <- unique(lib_size_large$Diversity)

		# factor
		lib_size_large$Treatment <- gsub("_", " ", lib_size_large$Treatment)
		similarity_scores$Treatment <- gsub("_", " ", similarity_scores$Treatment)
		lib_size_large$Treatment <- factor(lib_size_large$Treatment, levels=unique(lib_size_large$Treatment))
		similarity_scores$Treatment <- factor(similarity_scores$Treatment, levels=unique(similarity_scores$Treatment))

		# compute observed differences
		lib_size_large$Observed_Diff <- as.numeric(lib_size_large$observed_mean)-as.numeric(as.character(lib_size_large$Diversity))
		lib_size_large$Diversity <- factor(lib_size_large$Diversity, levels=c("50", "150", "250", "350", "450", "500", "600", "700", "800"))

		# Draw richness plots
		alpha_observed <- ggplot(lib_size_large, aes(x=Diversity, y=Observed_Diff, colour=Diversity)) +
#					geom_jitter(colour="black", size=1, alpha=0.5) +
					geom_boxplot(size=1, lwd=1) +
					scale_colour_manual(values=diversity_palette) +
					ylab("Difference in richness") +
					ggtitle(NULL) +
					facet_wrap(~Treatment, nrow=1, scale="free_x") +
					ggplot_theme(ang_le=45, leg_pos="bottom", fnt_sz=15, leg_size=1, top=0.2, right=0.2, bottom=0, left=0.2) +
					guides(colour=guide_legend(override.aes=list(shape=21, size=5), nrow=1)) +
					theme(axis.text.x=element_blank(), panel.spacing = unit(0.5, "lines"), axis.ticks.x=element_blank()) +
					xlab(NULL)

		# plot
		levels_palette <- friendly_pal("tol_eight")[1:7][c(1,2,3,6,7)]
		names(levels_palette) <- unique(similarity_scores$Level)[-1]
		
		shape_palette <- c(21:25)
		names(shape_palette) <- unique(similarity_scores$Level)[-1]
		
		diversity_plot_combo <- ggplot(similarity_scores[similarity_scores$Treatment!="TRUE" & similarity_scores$Diversity=="500" & similarity_scores$Library_Size=="Large" & similarity_scores$Level!="Kingdom", ], aes(x=Specificity, y=Sensitivity, shape=Level, fill=Level)) +
					geom_point(colour="black", size=2, stroke=0.5) +
					scale_shape_manual(values=shape_palette) +
					scale_fill_manual(values=levels_palette, guide = "none") +
					guides(fill=guide_legend(override.aes=list(size=5))) +
					ggtitle(NULL) +
					facet_wrap(~Treatment, nrow=1) +
					ggplot_theme(leg_pos="bottom", fnt_sz=15, ang_le=45, leg_size=1, top=0.2, right=0.2, bottom=0, left=0.2) +
					theme(panel.spacing = unit(0.5, "lines"))
		# merge plots
		connected_diversity <- cowplot::plot_grid(cowplot::plot_grid(alpha_observed, diversity_plot_combo, labels=c("A", "B"), nrow=2, align="vh", axis="tblr", label_size=20))

		# export for ITS2 only
		if (a_dir == "final_analysis_results_gITS7_ITS4") {
			export_final_svg(paste0(root_path, "tables_final/fig_2"), connected_diversity, base=F, width=168*2, height=237)
		}

		# compute observed differences
		combined_alphadiv_df_averaged$Observed_Diff <- as.numeric(combined_alphadiv_df_averaged$observed_mean)-as.numeric(as.character(combined_alphadiv_df_averaged$Diversity))
		combined_alphadiv_df_averaged$Diversity <- factor(combined_alphadiv_df_averaged$Diversity, levels=c("50", "150", "250", "350", "450", "500", "600", "700", "800"))

		# Draw richness plots
		alpha_observed_all <- ggplot(combined_alphadiv_df_averaged, aes(x=Diversity, y=Observed_Diff, fill=Diversity)) +
					geom_boxplot(size=0.5, lwd=1) +
					geom_jitter(colour="black", size=2, alpha=0.9) +
					scale_fill_manual(values=diversity_palette) +
					ylab("Difference in richness") +
					ggtitle(NULL) +
					facet_grid(rows=vars(Library_Size), cols=vars(Treatment), scales="free", space="free_x") +
					ggplot_theme(ang_le=45, leg_pos="none", fnt_sz=50) +
					guides(fill=guide_legend(override.aes=list(shape=21, size=20)))
		# supplementary
		count <- ifelse(a_dir=="final_analysis_results_ITS1f_ITS2", 10, 11)
		export_final_svg(paste0(supplementary_figs, "fig_S", count), alpha_observed_all, base=F, width=168*12, height=168*12)
	}
}
