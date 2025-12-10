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
# Rscript 000_stats_analysis_linear_models.R
# source("000_stats_analysis_linear_models.R")

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
library("nlme")

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

	print("linear modelling")

	##########################################################################################################################################################################

	if (a_dir == "final_analysis_results_Pauvert") {

		print(a_dir)

		# set factors
		similarity_scores$Treatment <- factor(similarity_scores$Treatment, levels=unique(similarity_scores$Treatment), ordered=F)
		similarity_scores$Level <- factor(similarity_scores$Level, levels=unique(similarity_scores$Level), ordered=F)
		colnames(similarity_scores) <- c("Diversity", "Library_Size", "Replicate", "Processing", "Method", "Treatment", "Level", "Sensitivity", "Specificity", "Aitchison", "Hellinger", "BrayCurtis")
		# some stats
		specificity_stats_full_model <- summary(lme(Specificity ~ Processing + Method, random=list(~1|Diversity, ~1|Level), data=similarity_scores))$tTable
		sensitivity_stats_full_model <- summary(lme(Sensitivity ~ Processing + Method, random=list(~1|Diversity, ~1|Level), data=similarity_scores))$tTable

		# linear modelling
		write.table(cbind.data.frame(Factor=rownames(specificity_stats_full_model), specificity_stats_full_model), paste0(save_stats, a_dir, "_specificity_stats_full_MIXED_EFFECTS.csv"), quote=F, sep="\t", row.names=F)
		write.table(cbind.data.frame(Factor=rownames(sensitivity_stats_full_model), sensitivity_stats_full_model), paste0(save_stats, a_dir, "_sensitivity_stats_full_MIXED_EFFECTS.csv"), quote=F, sep="\t", row.names=F)

		# loop through dists
		# declare empty
		pauv_all_full_models <- vector()
		for (a_dist in c("Aitchison", "BrayCurtis", "Hellinger")) {

#			a_dist <- "Aitchison"

			# define palette
			treatment_palette <- brewer.pal(7, "Dark2")
			names(treatment_palette) <- unique(similarity_scores$Treatment)

			# testing the existence of an effect due to stuff
			pauv_a_model <- broom::tidy((aov(get(a_dist)~Level+Treatment, data=similarity_scores)))
			pauv_full_model <- broom::tidy((aov(get(a_dist)~Level+Processing+Method, data=similarity_scores)))
			# store models
			pauv_all_full_models <- rbind.data.frame(pauv_all_full_models, cbind.data.frame(a_dist, rsquare=pauv_full_model$sumsq/sum(pauv_full_model$sumsq), pauv_full_model))
		}

		# plotting the stats
		pauv_plot_all_full <- pauv_all_full_models[pauv_all_full_models$term %in% c("Level", "Processing", "Method", "Library_Size"), ]
		# fix labels and stuff
		pauv_plot_all_full$term <- gsub("_", " ", pauv_plot_all_full$term)
		pauv_plot_all_full$term <- factor(pauv_plot_all_full$term, levels=unique(pauv_plot_all_full$term))
		colnames(pauv_plot_all_full) <- c("a_dist", "rsquare", "Term", "df", "sumsq", "meansq", "statistic", "p.value")
		
		# fix braycurtis
		pauv_plot_all_full$a_dist <- ifelse(pauv_plot_all_full$a_dist=="BrayCurtis", "Bray-Curtis", pauv_plot_all_full$a_dist)
		pauv_plot_all_full$a_dist <- factor(pauv_plot_all_full$a_dist, levels=unique(pauv_plot_all_full$a_dist))
		# export table
		write.table(pauv_plot_all_full, paste0(save_stats, "figure_5_stats.csv"), quote=F, sep="\t", row.names=F)

		# subset the data structure for plotting purposes
		similarity_scores_subset_by_lib <- similarity_scores[similarity_scores$Library_Size=="Large", ]
		# merge similarity scores
		subtable <- similarity_scores_subset_by_lib[, c("Diversity", "Library_Size", "Replicate", "Processing", "Method", "Treatment", "Level", "Sensitivity", "Specificity")]
		similarity_scores_subset_by_lib_long_fig_five <- cbind.data.frame(rbind.data.frame(subtable, subtable, subtable), similarity_score=c(similarity_scores_subset_by_lib$Aitchison, similarity_scores_subset_by_lib$Hellinger, similarity_scores_subset_by_lib$BrayCurtis), score_name=c(rep("Aitchison", nrow(similarity_scores_subset_by_lib)), rep("Hellinger", nrow(similarity_scores_subset_by_lib)), rep("BrayCurtis", nrow(similarity_scores_subset_by_lib))))
		
		# fix braycurtis
		similarity_scores_subset_by_lib_long_fig_five$score_name <- ifelse(similarity_scores_subset_by_lib_long_fig_five$score_name=="BrayCurtis", "Bray-Curtis", similarity_scores_subset_by_lib_long_fig_five$score_name)
		similarity_scores_subset_by_lib_long_fig_five$score_name <- factor(similarity_scores_subset_by_lib_long_fig_five$score_name, levels=c("Aitchison", "Bray-Curtis", "Hellinger"))

		# plot
		no_kingdom <- similarity_scores_subset_by_lib_long_fig_five[similarity_scores_subset_by_lib_long_fig_five$Level!="Kingdom", ]
		plot_sim_scores_pauv <- ggplot(no_kingdom, aes(x=Level, y=similarity_score, group=Treatment, colour=Treatment)) +
				stat_summary(position=position_dodge2(width=0.5), fun="mean", geom="line", lwd=1) +
				geom_point(data=no_kingdom[no_kingdom$Level != "Kingdom", ], aes(fill=Treatment, group=Treatment), colour="black", shape=21, size=2, alpha=0.6, position=position_dodge(width=0.5)) +
				scale_colour_manual(values=treatment_palette, labels=gsub("_", " ", levels(no_kingdom$Treatment))) +
				scale_fill_manual(values=treatment_palette, labels=gsub("_", " ", levels(no_kingdom$Treatment))) +
				ggtitle(NULL) +
				ggplot_theme(ang_le=45, leg_pos="bottom", fnt_sz=15, top=0.2, right=0.2, bottom=0, left=0.2) +
				guides(colour=guide_legend(override.aes=list(size=5), nrow=3)) +
				xlab(NULL) +
				ylab(NULL) +
				facet_grid(rows=vars(score_name), scales="free", space="free_x") +
				theme(panel.spacing = unit(0.7, "lines"))

		# export plot
		export_final_svg(paste0(root_path, "tables_final/fig_5"), plot_sim_scores_pauv, base=F, width=80*2, height=237)
		
	} else {

		print(a_dir)

		# define palette
		treatment_palette <- brewer.pal(7, "Dark2")
		names(treatment_palette) <- unique(similarity_scores$Treatment)

		# some stats
		specificity_stats_full_model <- summary(lme(Specificity ~ Processing + Method, random=list(~1|Diversity, ~1|Library_Size, ~1|Level), data=similarity_scores))$tTable
		sensitivity_stats_full_model <- summary(lme(Sensitivity ~ Processing + Method, random=list(~1|Diversity, ~1|Library_Size, ~1|Level), data=similarity_scores))$tTable

		# linear modelling
		write.table(cbind.data.frame(Factor=rownames(specificity_stats_full_model), specificity_stats_full_model), paste0(save_stats, a_dir, "_specificity_stats_full_MIXED_EFFECTS.csv"), quote=F, sep="\t", row.names=F)
		write.table(cbind.data.frame(Factor=rownames(sensitivity_stats_full_model), sensitivity_stats_full_model), paste0(save_stats, a_dir, "_sensitivity_stats_full_MIXED_EFFECTS.csv"), quote=F, sep="\t", row.names=F)

		# loop through diversities
		all_models <- vector()
		all_full_models <- vector()

		for (a_div in unique(similarity_scores$Diversity)) {

#			a_div <- 150

			# subset by diversity
			similarity_scores_subset_by_div <- similarity_scores[similarity_scores$Diversity == a_div, ]

			# define palette
			treatment_palette <- brewer.pal(7, "Dark2")
			names(treatment_palette) <- unique(similarity_scores$Treatment)

			# set factors
			similarity_scores_subset_by_div$Library_Size <- factor(similarity_scores_subset_by_div$Library_Size, levels=unique(similarity_scores_subset_by_div$Library_Size), ordered=F)
			similarity_scores_subset_by_div$Treatment <- factor(similarity_scores_subset_by_div$Treatment, levels=unique(similarity_scores_subset_by_div$Treatment), ordered=F)
			similarity_scores_subset_by_div$Level <- factor(similarity_scores_subset_by_div$Level, levels=unique(similarity_scores_subset_by_div$Level), ordered=F)
			colnames(similarity_scores_subset_by_div) <- c("Diversity", "Library_Size", "Replicate", "Processing", "Method", "Treatment", "Level", "Sensitivity", "Specificity", "Aitchison", "Hellinger", "BrayCurtis")
			# loop through dists
			# count for plot exporting
			count <- ifelse(a_dir=="final_analysis_results_ITS1f_ITS2", 20, 21)
			for (a_dist in c("Aitchison", "BrayCurtis", "Hellinger")) {

#				a_dist <- "Aitchison"

				print(paste("analysing", a_div, a_dist))

				# testing the existence of an effect due to stuff
				a_model <- broom::tidy(aov(get(a_dist)~Level+Treatment+Library_Size, data=similarity_scores_subset_by_div))
				# store models
				all_models <- rbind.data.frame(all_models, cbind.data.frame(a_div, a_dist, rsquare=a_model$sumsq/sum(a_model$sumsq), a_model))

				# full model
				full_model <- broom::tidy(aov(get(a_dist)~Level+Processing+Method+Library_Size, data=similarity_scores_subset_by_div))
				# store models
				all_full_models <- rbind.data.frame(all_full_models, cbind.data.frame(a_div, a_dist, rsquare=full_model$sumsq/sum(full_model$sumsq), full_model))

#####################################				# export stats
#####################################				write.table(a_model, paste0(save_stats, a_dir, "_", a_dist, "_short_model_.csv"), quote=F, sep="\t", row.names=F)
#####################################				write.table(full_model, paste0(save_stats, a_dir, "_", a_dist, "_full_model_.csv"), quote=F, sep="\t", row.names=F)

				# plot
				plot_dist <- ggplot(similarity_scores[similarity_scores$Level != "Kingdom", ], aes(x=Level, y=get(a_dist), group=Treatment, colour=Treatment)) +
							stat_summary(position=position_dodge2(width=0.5), fun="mean", geom="line", lwd=2) +
							stat_summary(position=position_dodge2(width=0.5), fun.data="mean_sdl", fun.args=list(mult=1), geom="errorbar", width=0.5, lwd=2) +
							scale_colour_manual(values=treatment_palette) +
							ggtitle(NULL) +
							facet_grid(rows = vars(Diversity), cols = vars(Library_Size)) +
							ggplot_theme(ang_le=45, leg_pos="bottom") +
							xlab(NULL)
				# fix ylab for braycurtis
				if (a_dist == "BrayCurtis") {
					plot_dist <- plot_dist + ylab("Bray-Curtis")
				} else {
					plot_dist <- plot_dist + ylab(a_dist)
				}
				# export treatment
				export_final_svg(paste0(supplementary_figs, "fig_S", count), plot_dist, base=F, width=168*6, height=168*6)
				count <- count + 2
			}
		}

		# plotting the stats
		plot_all <- all_models[all_models$term %in% c("Level", "Treatment", "Library_Size"), ]

		plot_all$term <- gsub("_", " ", plot_all$term)
		plot_all$term <- factor(plot_all$term, levels=unique(plot_all$term))
		colnames(plot_all) <- c("a_div", "a_dist", "rsquare", "Term", "df", "sumsq", "meansq", "statistic", "p.value")

		# set levels
		plot_all$a_div <- factor(plot_all$a_div, levels=c("50", "150", "250", "350", "450", "500", "600", "700", "800"))

		# plot
		plot_all_models <- ggplot(plot_all, aes(x=Term, y=rsquare)) +
					geom_boxplot(size=0.5, lwd=1, outlier.shape=NA) +
					ggnewscale::new_scale_fill() +
					geom_point(data=plot_all, aes(fill=a_div), colour="black", shape=21, size=2, alpha=0.8, position="jitter") +
					scale_fill_manual(name="Diversity", values=friendly_pal("muted_nine")[1:9], guide = "none") +
					ggplot_theme(ang_le=45, leg_pos="right", fnt_sz=15) +
					xlab(NULL) +
					scale_y_continuous(breaks=seq(0, 1, by=0.15)) +
					ggtitle("Short model") +
					ylab(bquote(paste("R", ""^2))) +
					facet_wrap(~a_dist) +
					theme(panel.spacing = unit(0.5, "lines"))
		# plotting the stats
		plot_all_full <- all_full_models[all_full_models$term %in% c("Level", "Processing", "Method", "Library_Size"), ]
		# fix labels and stuff
		plot_all_full$a_div <- factor(plot_all_full$a_div, levels=c("50", "150", "250", "350", "450", "500", "600", "700", "800"))
		plot_all_full$term <- gsub("_", " ", plot_all_full$term)
		plot_all_full$term <- factor(plot_all_full$term, levels=unique(plot_all_full$term))
		colnames(plot_all_full) <- c("a_div", "a_dist", "rsquare", "Term", "df", "sumsq", "meansq", "statistic", "p.value")
		
		# fix braycurtis
		plot_all_full$a_dist <- ifelse(plot_all_full$a_dist=="BrayCurtis", "Bray-Curtis", plot_all_full$a_dist)
		plot_all_full$a_dist <- factor(plot_all_full$a_dist, levels=unique(plot_all_full$a_dist))

		# plot
		plot_all_full_models <- ggplot(plot_all_full, aes(x=Term, y=rsquare)) +
						geom_boxplot(size=0.5, lwd=1, outlier.shape=NA) +
						ggnewscale::new_scale_fill() +
						geom_point(data=plot_all_full, aes(fill=a_div), colour="black", shape=21, size=2, alpha=0.8, position="jitter") +
						scale_fill_manual(name="Diversity", values=friendly_pal("muted_nine")[1:9]) +
						guides(fill=guide_legend(override.aes=list(shape=21, size=5), nrow=1)) +
						ggplot_theme(ang_le=45, leg_pos="right", fnt_sz=15, top=0.2, right=0.2, bottom=0.2, left=0.2) +
						scale_y_continuous(breaks=seq(0, 1, by=0.15)) +
						xlab(NULL) +
						ggtitle("Full model") +
						ylab(bquote(paste("R", ""^2))) +
						facet_wrap(~a_dist) +
						theme(panel.spacing = unit(0.5, "lines"))
		# get legend
		legend <- get_legend(plot_all_full_models)

		# final plot
#		final_plot_model <- ggarrange(plot_all_models, plot_all_full_models, ncol=2, labels=c('Short model', 'Full model'), common.legend=T, font.label=list(size=30, colour="black"), legend="bottom", )
		final_plot_model <- cowplot::plot_grid(
					cowplot::plot_grid(plot_all_models + theme(legend.position="none"), plot_all_full_models + theme(axis.title.y=element_blank(), legend.position="none"), ncol=2, align="vh", axis="tblr"),
					legend, nrow=2, rel_heights=c(1, 0.1))

		# export only for ITS2
		if (a_dir == "final_analysis_results_gITS7_ITS4") {
			# export treatment
			export_final_svg(paste0(root_path, "tables_final/fig_4"), final_plot_model, base=F, width=168*2, height=237)
			# export table
			write.table(plot_all_full, paste0(save_stats, "figure_4_stats.csv"), quote=F, sep="\t", row.names=F)
		} else if (a_dir == "final_analysis_results_ITS1f_ITS2") {
			# export treatment
			export_final_svg(paste0(supplementary_figs, "fig_S26"), final_plot_model, base=F, width=168*2, height=237)
			# export table
			write.table(plot_all_full, paste0(save_stats, "figure_S26_stats.csv"), quote=F, sep="\t", row.names=F)
		}
	}
}

