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
# Rscript 000_stats_analysis_distances.R
# source("000_stats_analysis_distances.R")

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

	print("plot distances")

	##########################################################################################################################################################################

	if (a_dir == "final_analysis_results_Pauvert") {

		# define palette
		treatment_palette <- brewer.pal(7, "Dark2")
		names(treatment_palette) <- unique(combined_alphadiv_df_averaged$Treatment)

		# plot shannon and simpson indexes
		shannon <- ggplot(combined_alphadiv_df_averaged, aes(x=Treatment, y=shannon_mean, fill=Treatment)) +
				geom_point(colour="black", size=10, shape=21) +
				scale_fill_manual(values=treatment_palette, labels=gsub("_", " ", levels(combined_alphadiv_df_averaged$Treatment))) +
				scale_y_continuous(limits = c(0, NA)) +
				ggplot_theme(leg_pos="none") +
				facet_wrap(~Treatment, nrow=1, scale="free_x") +
				ylab("Shannon index") +
				xlab(NULL) +
				theme(axis.text.x=element_blank())

		simpson <- ggplot(combined_alphadiv_df_averaged, aes(x=Treatment, y=simpson_mean, fill=Treatment)) +
				geom_point(colour="black", size=10, shape=21) +
				scale_fill_manual(values=treatment_palette, labels=gsub("_", " ", levels(combined_alphadiv_df_averaged$Treatment))) +
				scale_y_continuous(limits = c(0, NA)) +
				ggplot_theme(leg_pos="none") +
				facet_wrap(~Treatment, nrow=1, scale="free_x") +
				ylab("Simpson index") +
				xlab(NULL) +
				theme(axis.text.x=element_blank())

		richness <- ggplot(combined_alphadiv_df_averaged, aes(x=Treatment, y=observed_mean, fill=Treatment)) +
				geom_point(colour="black", size=10, shape=21) +
				scale_fill_manual(values=treatment_palette, labels=gsub("_", " ", levels(combined_alphadiv_df_averaged$Treatment))) +
				scale_y_continuous(limits = c(0, 200), breaks = waiver(), n.breaks=5) +
				ggplot_theme(leg_pos="none") +
				facet_wrap(~Treatment, nrow=1, scale="free_x") +
				ylab("Estimated Richness") +
				xlab(NULL) +
				theme(axis.text.x=element_blank())
		# export plot
		export_final_svg(paste0(supplementary_figs, "fig_S32"), richness, base=F, width=168*8, height=168*4)
		export_final_svg(paste0(supplementary_figs, "fig_S34"), shannon, base=F, width=168*8, height=168*4)
		export_final_svg(paste0(supplementary_figs, "fig_S35"), simpson, base=F, width=168*8, height=168*4)
	} else {
		# define palette
		diversity_palette <- brewer.pal(9, "Paired")
		names(diversity_palette) <- unique(combined_alphadiv_df_averaged$Diversity)

		# plot shannon and simpson indexes
		shannon <- ggplot(combined_alphadiv_df_averaged, aes(x=Diversity, y=shannon_mean, fill=Diversity)) +
				geom_boxplot(size=4, lwd=1) +
				geom_jitter(colour="black", size=2, alpha=0.9) +
				scale_fill_manual(values=diversity_palette, labels=gsub("_", " ", levels(combined_alphadiv_df_averaged$Diversity))) +
				scale_y_continuous(limits = c(0, NA)) +
				guides(fill=guide_legend(override.aes=list(shape=21, size=10), nrow=1)) +
				ggplot_theme(leg_pos="bottom") +
				theme(axis.text.x=element_blank()) +
				ylab("Shannon index") +
				xlab(NULL) +
				facet_wrap(~Library_Size+Treatment, ncol=7)

		simpson <- ggplot(combined_alphadiv_df_averaged, aes(x=Diversity, y=simpson_mean, fill=Diversity)) +
				geom_boxplot(size=4, lwd=1) +
				geom_jitter(colour="black", size=2, alpha=0.9) +
				scale_fill_manual(values=diversity_palette, labels=gsub("_", " ", levels(combined_alphadiv_df_averaged$Diversity))) +
				scale_y_continuous(limits = c(0, NA)) +
				guides(fill=guide_legend(override.aes=list(shape=21, size=10), nrow=1)) +
				ggplot_theme(leg_pos="bottom") +
				theme(axis.text.x=element_blank()) +
				ylab("Simpson index") +
				xlab(NULL) +
				facet_wrap(~Library_Size+Treatment, ncol=7)

		richness <- ggplot(combined_alphadiv_df_averaged, aes(x=Diversity, y=observed_mean, fill=Diversity)) +
				geom_boxplot(size=4, lwd=1) +
				geom_jitter(colour="black", size=2, alpha=0.9) +
				scale_fill_manual(values=diversity_palette, labels=gsub("_", " ", levels(combined_alphadiv_df_averaged$Diversity))) +
				scale_y_continuous(limits = c(0, NA)) +
				guides(fill=guide_legend(override.aes=list(shape=21, size=10), nrow=1)) +
				ggplot_theme(leg_pos="bottom") +
				theme(axis.text.x=element_blank()) +
				ylab("Estimated Richness") +
				xlab(NULL) +
				facet_wrap(~Library_Size+Treatment, ncol=7)

		# export plot
		count <- ifelse(a_dir=="final_analysis_results_ITS1f_ITS2", 8, 9)
		export_final_svg(paste0(supplementary_figs, "fig_S", count), richness, base=F, width=168*12, height=168*8)
		export_final_svg(paste0(supplementary_figs, "fig_S", count+4), shannon, base=F, width=168*12, height=168*8)
		export_final_svg(paste0(supplementary_figs, "fig_S", count+6), simpson, base=F, width=168*12, height=168*8)
	}

	# get table ready for all averages
	final_alphadiv <- vector()

	# compute averages for jason
	for (l_size in unique(combined_alphadiv_df_averaged$Library_Size)) {
		# subset by lib size
		lib_size_sub <- combined_alphadiv_df_averaged[combined_alphadiv_df_averaged$Library_Size==l_size, ]
		# continue looping
		for (a_treat in unique(lib_size_sub$Treatment)) {
			for (a_div in unique(lib_size_sub$Diversity)) {

				# get subsetted data
				subset_alpha <- lib_size_sub[lib_size_sub$Treatment == a_treat & lib_size_sub$Diversity == a_div, ]
				# compute mean
				avg_diff <- mean(subset_alpha$observed_mean)
				median_diff <- median(subset_alpha$observed_mean)
				std_err <- plotrix::std.error(subset_alpha$observed_mean)
				# put all together
				final_alphadiv <- rbind.data.frame(final_alphadiv, cbind.data.frame(lib_size=l_size, treatment=a_treat, diversity=a_div, avg_obs=avg_diff, median_obs=median_diff, se=std_err, se_up=avg_diff+std_err, se_down=avg_diff-std_err))
			}
		}
	}
	# export table
	write.table(final_alphadiv, paste0(save_abundances, a_dir, "_diversity_table.csv"), quote=F, sep="\t", row.names=F)
}
