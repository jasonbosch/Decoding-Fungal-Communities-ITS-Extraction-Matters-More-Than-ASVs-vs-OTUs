# get best taxa table
barplot_supplementary <- function(phylo_obj, taxa_level, count) {

	# debug
	# phylo_obj <- true_its ; taxa_level <- "Genus" ; a_marker <- "ITS1"

	# empty table
	psmelted_phylo <- list()

	# loop through phylos
	for (a_marker in names(phylo_obj)){

		# get meta
		meta_phylo <- sample_data(phylo_obj[[a_marker]])

		# subset by Large
		pruned_phylo <- prune_samples(rownames(meta_phylo[meta_phylo$Library_Size=="Large", ]), phylo_obj[[a_marker]])

		# remove ASVs which are always zero for the subset under consideration
		phylo_clean <- prune_taxa(taxa_sums(pruned_phylo) > 0, pruned_phylo)

		# aggregate at order level
		aggregated_phylo <- tax_glom(transform_sample_counts(phylo_clean, function(x) x / sum(x)), taxrank=taxa_level, NArm=FALSE)
		meta_aggregated <- sample_data(aggregated_phylo)

		# export table
		abundance_table <- cbind.data.frame(tax_table(aggregated_phylo)[, taxa_level], rowMeans(otu_table(aggregated_phylo)))
		colnames(abundance_table) <- c("taxa_name", "avg_abundance")
		# rename empties
		abundance_table$taxa_name[abundance_table$taxa_name==""] <- "Unidentified"
		# get identified only
		identi_only <- abundance_table[abundance_table$taxa_name!="Unidentified", ]
		# merge unidenti
		back_together <- rbind.data.frame(identi_only, cbind.data.frame(taxa_name="Unidentified", avg_abundance=sum(abundance_table$avg_abundance[abundance_table$taxa_name=="Unidentified"])))
		write.table(back_together, paste0(save_abundances, "abundance_table_", taxa_level, "_", a_marker, ".csv"), quote=F, sep="\t", row.names=F)

		# # # # # # # # # # #  THIS IS THE PART WHERE THE BEST TAXA ARE RETRIEVED # # # # # # # # # # # 

		# remove unidentified from taxa_names
		unidentified_taxa <- back_together$taxa_name[grep("Unidentified", back_together$taxa_name)]
		identified_taxa <- back_together[grep("Unidentified", back_together$taxa_name, invert=T), ]

		# get best if greater than 0.01
		best_taxa <- identified_taxa$taxa_name[identified_taxa$avg_abundance > 0.01]
		# get others, otherwise
		other_orders <- identified_taxa$taxa_name[identified_taxa$avg_abundance <= 0.01]

		# psmelt phylo_sub
		psmelted_phylo[[a_marker]] <- psmelt(aggregated_phylo) %>%
					group_by(get(taxa_level), Diversity, Replicate) %>% # group by
					summarise(abundance=mean(Abundance)) %>% # and use the groups to compute abundances
					arrange(desc(abundance)) %>% # sort by
					data.frame() %>% # turn into data frame
					rename(t_level=get.taxa_level.) # rename weird column name

		# make new columns for plotting purposes
		psmelted_phylo[[a_marker]]$plot_order <- ifelse(psmelted_phylo[[a_marker]]$t_level %in% other_orders, "Other", ifelse(psmelted_phylo[[a_marker]]$t_level=="", "Unidentified", psmelted_phylo[[a_marker]]$t_level))
	}
	# get all taxa to plot
	all_taxa <- unique(unlist(lapply(psmelted_phylo, "[[", 5)))

	# set factors
	last_two <- as.character(c(all_taxa[!grepl("Other|Unidenti", all_taxa)], all_taxa[grepl("Other", all_taxa)], all_taxa[grepl("Unidenti", all_taxa)]))
	psmelted_phylo[["ITS1"]]$plot_order <- factor(psmelted_phylo[["ITS1"]]$plot_order, levels=unique(last_two))
	psmelted_phylo[["ITS2"]]$plot_order <- factor(psmelted_phylo[["ITS2"]]$plot_order, levels=unique(last_two))

	# set palette
	taxa_palette <- colorRampPalette(brewer.pal(8, "Paired"))(length(all_taxa))
	names(taxa_palette) <- unique(last_two)

	# loop through markers
	ribon_plot <- list()
	# get the plot ready
	ribon_plot[["ITS1"]] <- ggplot(psmelted_phylo[["ITS1"]], aes(x=Replicate)) +
					geom_bar(aes(weight=abundance, fill=plot_order), position="fill", color="black", show.legend=TRUE) +
					scale_fill_manual(values=taxa_palette, drop=FALSE) + # drop=FALSE forces ggplot to show all data in the legend
					guides(fill=guide_legend(ncol=1, title=taxa_level)) +
					ggplot_theme(leg_pos="right", ang_le=0, fnt_sz=50) +
					xlab(NULL) +
					ylab("Relative abundance") +
					ggtitle("ITS1") +
					facet_wrap(~Diversity, ncol=1) +
					theme(plot.margin=unit(c(1, 0.5, 0, 0.5), "cm"))
	# get the other plot ready
	ribon_plot[["ITS2"]] <- ggplot(psmelted_phylo[["ITS2"]], aes(x=Replicate)) +
					geom_bar(aes(weight=abundance, fill=plot_order), position="fill", color="black", show.legend=TRUE) +
					scale_fill_manual(values=taxa_palette, drop=FALSE) + # drop=FALSE forces ggplot to show all data in the legend
					ggplot_theme(leg_pos="right", ang_le=0, fnt_sz=50) +
					guides(fill=guide_legend(ncol=1, title=taxa_level)) +
					xlab(NULL) +
					ylab(NULL) +
					theme(axis.text.y=element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(1, 0.5, 0, 0.5), "cm")) +
					ggtitle("ITS2") +
					facet_wrap(~Diversity, ncol=1)
	# get legend
	legend <- get_legend(ribon_plot[["ITS1"]])

	# put plot together
	combo_plot <- cowplot::plot_grid(ribon_plot[["ITS1"]] + theme(legend.position="none"), ribon_plot[["ITS2"]] + theme(legend.position="none"), legend, ncol=3, align="vh", axis="tblr", rel_widths=c(1, 1, 0.5))
	export_final_svg(paste0(supplementary_figs, "fig_S", count), combo_plot, base=F, width=168*12, height=168*10)
}

# get best taxa table
barplot_main_phylum <- function(phylo_obj, taxa_level) {

	# debug
	# phylo_obj <- true_its ; taxa_level <- "Order"

	# merge phylos
	ITS1 <- subset_samples(phylo_obj[["ITS1"]], Library_Size == "Large")
	ITS2 <- subset_samples(phylo_obj[["ITS2"]], Library_Size == "Large")

	# Abundance tables
	Combined_Abundance_Table <- list()
	for (MARKER in c("ITS1","ITS2")) {
		#Tax glom to level of Phylum and get relative abundance
		Physeq_Current_Phylum <- tax_glom(get(MARKER), taxrank = "Phylum",NArm = FALSE)
		Physeq_Current_Phylum <- transform_sample_counts(Physeq_Current_Phylum, function (x) x/sum(x))
		Physeq_Current_Phylum_Table <- psmelt(Physeq_Current_Phylum)
		Physeq_Current_Phylum_Table$ID <- paste(Physeq_Current_Phylum_Table$Sample,Physeq_Current_Phylum_Table$Phylum,sep="_")
		colnames(Physeq_Current_Phylum_Table)[3] <- paste("Abundance",MARKER,sep = "_")
		
		#Tax glom to level of order and merge taxa with RA < 1%
		Physeq_Current_Order <- tax_glom(get(MARKER), taxrank = "Order",NArm = FALSE)
		Physeq_Current_Order <- transform_sample_counts(Physeq_Current_Order, function (x) x/sum(x))
		Physeq_Current_Order_Table <- psmelt(Physeq_Current_Order)
		Physeq_Current_Order_Table$ID <- paste(Physeq_Current_Order_Table$Sample,Physeq_Current_Order_Table$Order,sep="_")
		colnames(Physeq_Current_Order_Table)[3] <- paste("Abundance",MARKER,sep = "_")
		
		#Save to list
		Combined_Abundance_Table[paste(MARKER,"Phylum",sep = "")] <- list(Physeq_Current_Phylum_Table)
		Combined_Abundance_Table[paste(MARKER,"Order",sep = "")] <- list(Physeq_Current_Order_Table)
	}

	# Get Average Abundance
	Average_Abundance_Phylum <- merge.data.frame(Combined_Abundance_Table$ITS1Phylum,Combined_Abundance_Table$ITS2Phylum,by="ID")
	Average_Abundance_Phylum$Abundance <- rowMeans(Average_Abundance_Phylum[,c("Abundance_ITS1","Abundance_ITS2")])
	Average_Abundance_Order <- merge.data.frame(Combined_Abundance_Table$ITS1Order,Combined_Abundance_Table$ITS2Order,by="ID")
	Average_Abundance_Order$Abundance <- rowMeans(Average_Abundance_Order[,c("Abundance_ITS1","Abundance_ITS2")])

	# Get average abundance, mark rare phyla, summarise table and set phylum levels
	Sorted_Phylum <- aggregate(Abundance ~ Phylum.x,Average_Abundance_Phylum[,c("Phylum.x","Abundance")],mean)
	Sorted_Phylum[which(Sorted_Phylum$Abundance<0.01),"Phylum.x"] <- "Rare"
	Sorted_Phylum <- aggregate(Abundance ~ Phylum.x,Sorted_Phylum,sum)
	Sorted_Phylum <- Sorted_Phylum[order(-Sorted_Phylum$Abundance),]
	Sorted_Phylum$Phylum.x <- factor(Sorted_Phylum$Phylum.x, levels = Sorted_Phylum$Phylum.x, ordered = T)

	# Set the Phyla in the Order table to match which ones are not rare and find average abundance
	Average_Abundance_Order[!Average_Abundance_Order$Phylum.x%in%Sorted_Phylum$Phylum.x,"Phylum.x"] <- "Rare"
	Sorted_Order <- aggregate(Abundance ~ Phylum.x*Order.x,Average_Abundance_Order[,c("Phylum.x","Order.x","Abundance")],mean)
	# Change the order to rare (if rare phyla) or other X for other phyla
	Sorted_Order[which(Sorted_Order$Phylum.x=="Rare"),"Order.x"] <- "Rare"
	for (ROW in 1:nrow(Sorted_Order)) {
		if (Sorted_Order[ROW,"Abundance"]<0.01 & Sorted_Order[ROW,"Phylum.x"]!="Rare") {
			Sorted_Order[ROW,"Order.x"] <- paste("Other",Sorted_Order[ROW,"Phylum.x"],sep = "_")
		}
	}
	Sorted_Order <- aggregate(Abundance ~ Phylum.x*Order.x,Sorted_Order[,c("Phylum.x","Order.x","Abundance")],sum)
	# Get the order based on Phylum and Order abundances
	Sorted_Order_values <- NULL
	for (PHYLUM in levels(Sorted_Phylum$Phylum.x)) {
		Current_Sort <- Sorted_Order[Sorted_Order$Phylum.x==PHYLUM,]
		Current_Sort <- Current_Sort[order(-Current_Sort$Abundance),]
		if (paste("Other",PHYLUM,sep = "_")%in%Current_Sort$Order.x) {
			Sorted_Order_values <- c(Sorted_Order_values,Current_Sort$Order.x[which(Current_Sort$Order.x!=paste("Other",PHYLUM,sep = "_"))],paste("Other",PHYLUM,sep = "_"))
		} else {
			Sorted_Order_values <- c(Sorted_Order_values,Current_Sort$Order.x)
		}
	}
	Sorted_Order$Order.x <- factor(Sorted_Order$Order.x, levels = Sorted_Order_values, ordered = T)

	# Set palettes
	palette_Phylum <- friendly_pal("bright_seven")[1:(length(Sorted_Phylum$Phylum.x))]
	names(palette_Phylum) <- levels(Sorted_Phylum$Phylum.x)
	palette_order <- colorRampPalette(brewer.pal(12, "Paired"))(length(Sorted_Order$Order.x))
	names(palette_order) <- levels(Sorted_Order$Order.x)

	# Plot both Phylum and order on one plot
	Combined_Plot <- list()
	for (MARKER in c("ITS1","ITS2")) {

		#Tax glom to level of Phylum and name the rares
		Physeq_Current_Phylum <- tax_glom(get(MARKER), taxrank = "Phylum",NArm = FALSE)
		Physeq_Current_Phylum <- transform_sample_counts(Physeq_Current_Phylum, function (x) x/sum(x))
		tax_table(Physeq_Current_Phylum)[which(!tax_table(Physeq_Current_Phylum)[,"Phylum"]%in%levels(Sorted_Phylum$Phylum.x)),"Phylum"] <- "Rare"
		Physeq_Current_Phylum <- tax_glom(Physeq_Current_Phylum, taxrank = "Phylum",NArm = FALSE)
		Physeq_Current_Phylum_Table <- psmelt(Physeq_Current_Phylum)
		Physeq_Current_Phylum_Table$Phylum <- factor(Physeq_Current_Phylum_Table$Phylum, levels=levels(Sorted_Phylum$Phylum.x), ordered = T)

		#Tax glom to level of order and copy the names
		Physeq_Current_Order <- tax_glom(get(MARKER), taxrank = "Order",NArm = FALSE)
		Physeq_Current_Order <- transform_sample_counts(Physeq_Current_Order, function (x) x/sum(x))
		tax_table(Physeq_Current_Order)[which(!tax_table(Physeq_Current_Order)[,"Phylum"]%in%levels(Sorted_Phylum$Phylum.x)),"Phylum"] <- "Rare"
		for (ROW in 1:nrow(tax_table(Physeq_Current_Order))) {
			if (tax_table(Physeq_Current_Order)[ROW,"Phylum"]=="Rare") {
				tax_table(Physeq_Current_Order)[ROW,"Order"] <- "Rare"
			}
			if (!tax_table(Physeq_Current_Order)[ROW,"Order"]%in%levels(Sorted_Order$Order.x)) {
				tax_table(Physeq_Current_Order)[ROW,"Order"] <- paste("Other",tax_table(Physeq_Current_Order)[ROW,"Phylum"],sep = "_")
			}
		}
		tax_table(Physeq_Current_Order)[,"Class"] <- tax_table(Physeq_Current_Order)[,"Order"]
		Physeq_Current_Order <- tax_glom(Physeq_Current_Order, taxrank = "Order",NArm = FALSE)
		Physeq_Current_Order_Table <- psmelt(Physeq_Current_Order)
		Physeq_Current_Order_Table$Order <- factor(Physeq_Current_Order_Table$Order, levels=levels(Sorted_Order$Order.x), ordered = T)

		Plot_Combined <- 
			ggplot() +
			geom_col(data = Physeq_Current_Phylum_Table,aes(x = Replicate, y = Abundance, fill=Phylum), position = "fill", color="black", width=0.3, just=1.55) +
			scale_fill_manual(name="Phylum", values=palette_Phylum, drop=FALSE) +
			guides(fill=guide_legend(ncol=1,order = 1)) +
			new_scale_fill() +
			geom_col(data = Physeq_Current_Order_Table, aes(x = Replicate, y = Abundance, fill=Order), position = "fill", color="black", width=0.6, just=0.22) +
			scale_fill_manual(name="Order", values=palette_order, drop=FALSE) +
			guides(fill=guide_legend(ncol=1,order = 2)) +
			facet_wrap(~Diversity,ncol = 1, scales = "free_x") +
			labs(x= "Replicate", y = "Relative Abundance") +
			theme_bw() +
			theme(text=element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),legend.position = "right",axis.text = element_text(size=16))

		Combined_Plot[[MARKER]] <- Plot_Combined
	}


	# get legend
	combo_legend <- get_legend(Combined_Plot[["ITS1"]])
	# put plot together
	combo_plot <- plot_grid(Combined_Plot[["ITS1"]] + theme(legend.position="none"), 
			Combined_Plot[["ITS2"]] + theme(legend.position="none"), 
			combo_legend, ncol=3, align="vh", axis="tblr", rel_widths=c(1, 1, 0.8))
	# export
	export_final_svg(paste0(root_path, "tables_final/SUPER_ULTRA_NEW_fig_1"), combo_plot, base=F, width=168*2, height=237*2)

}

# get best taxa table
barplot_supplementary_pauvert <- function(phylo_obj, taxa_level, treatment, count, plotting=FALSE) {

	# debug
	# phylo_obj <- true_its[["Pauvert_combined"]] ; taxa_level <- "Genus" ; treatment <- "TRUE" ; plotting <- TRUE ;

	# get meta
	meta_phylo <- sample_data(phylo_obj)

	# remove ASVs which are always zero for the subset under consideration
	phylo_treat <- prune_samples(as.character(rownames(meta_phylo[which(meta_phylo$Treatment == treatment)])), phylo_obj)
	# remove zeroes
	phylo_clean <- prune_taxa(taxa_sums(phylo_treat) > 0, phylo_treat)

	# aggregate at order level
	aggregated_phylo <- tax_glom(transform_sample_counts(phylo_clean, function(x) x / sum(x)), taxrank=taxa_level, NArm=FALSE)
	meta_aggregated <- sample_data(aggregated_phylo)

	# export table
	abundance_table <- cbind.data.frame(tax_table(aggregated_phylo)[, taxa_level], rowMeans(otu_table(aggregated_phylo)))
	colnames(abundance_table) <- c("taxa_name", "avg_abundance")
	# rename empties
	abundance_table$taxa_name[abundance_table$taxa_name==""] <- "Unidentified"
	# get identified only
	identi_only <- abundance_table[abundance_table$taxa_name!="Unidentified", ]
	# merge unidenti
	back_together <- rbind.data.frame(identi_only, cbind.data.frame(taxa_name="Unidentified", avg_abundance=sum(abundance_table$avg_abundance[abundance_table$taxa_name=="Unidentified"])))

	# # # # # # # # # # #  THIS IS THE PART WHERE THE BEST TAXA ARE RETRIEVED # # # # # # # # # # # 

	if (plotting==TRUE) {
		# remove unidentified from taxa_names
		unidentified_taxa <- back_together$taxa_name[grep("Unidentified", back_together$taxa_name)]
		identified_taxa <- back_together[grep("Unidentified", back_together$taxa_name, invert=T), ]

		# get best if greater than 0.01
		best_taxa <- identified_taxa$taxa_name[identified_taxa$avg_abundance > 0.01]
		# get others, otherwise
		other_orders <- identified_taxa$taxa_name[identified_taxa$avg_abundance <= 0.01]

		# psmelt phylo_sub
		psmelted_phylo <- psmelt(aggregated_phylo) %>%
					group_by(get(taxa_level), Replicate) %>% # group by
					summarise(abundance=mean(Abundance)) %>% # and use the groups to compute abundances
					arrange(desc(abundance)) %>% # sort by
					data.frame() %>% # turn into data frame
					rename(t_level=get.taxa_level.) # rename weird column name

		# make new columns for plotting purposes
		psmelted_phylo$plot_order <- ifelse(psmelted_phylo$t_level %in% other_orders, "Other", ifelse(psmelted_phylo$t_level=="", "Unidentified", psmelted_phylo$t_level))

		# get all taxa to plot
		all_taxa_pauvert <- unique(psmelted_phylo$plot_order)

		# set factors
		last_two <- as.character(c(all_taxa_pauvert[!grepl("Other|Unidenti", all_taxa_pauvert)], all_taxa_pauvert[grepl("Other", all_taxa_pauvert)], all_taxa_pauvert[grepl("Unidentified", all_taxa_pauvert)]))
		psmelted_phylo$plot_order <- factor(psmelted_phylo$plot_order, levels=unique(last_two))

		# set palette
		taxa_palette <- colorRampPalette(brewer.pal(8, "Paired"))(length(all_taxa_pauvert))
		names(taxa_palette) <- unique(last_two)

		# get the plot ready
		ribon_plot <- ggplot(psmelted_phylo, aes(x=Replicate)) +
						geom_bar(aes(weight=abundance, fill=plot_order), position="fill", color="black") +
						scale_fill_manual(values=taxa_palette, drop=FALSE) + # drop=FALSE forces ggplot to show all data in the legend
						guides(names=taxa_level, fill=guide_legend(ncol=1)) +
						ggplot_theme(leg_pos="right", ang_le=0, fnt_sz=50) +
						xlab(NULL) +
						ylab("Relative abundance") +
						ggtitle("Pauvert") +
						labs(fill=taxa_level) +
						theme(plot.margin=unit(c(1, 0.5, 0, 0.5), "cm"))

		export_final_svg(paste0(supplementary_figs, "fig_S", count), ribon_plot, base=F, width=168*4, height=168*8)
	}
	# return table
	return(back_together)
}

# mm to inch
mm_to_inch <- function(mm) mm / 25.4

# export images for publication, A4 by default
export_final_svg <- function(filename, a_plot, width=420, height=394, base=F) {

	# if plot is made with ggplot2
	if (base==F) {
		# plot
		svg(paste0(filename, ".svg"), width=mm_to_inch(width), height=mm_to_inch(height))
			plot(a_plot)
		dev.off()
	# otherwise, it was made with base R
	} else {
		# export plot as svg
		svg(paste0(filename, ".svg"), width=mm_to_inch(width), height=mm_to_inch(height))
			a_plot
		dev.off()
	}
}

# create consistent theme for all ggplots
ggplot_theme <- function(leg_pos="right", ang_le=0, fnt_sz=50, leg_size=3, top=1, right=0.5, bottom=1, left=0.5) {

	if (ang_le == 45) {
		# set theme
		theme_bw() +
		theme(text=element_text(size=fnt_sz, color="black"), legend.position=leg_pos, legend.text=element_text(size=fnt_sz), legend.key.height=unit(1,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=ang_le, hjust=1, vjust=1, color="black", size=fnt_sz), axis.text.y=element_text(color="black", size=fnt_sz), strip.text=element_text(size=fnt_sz), panel.spacing=unit(3, "lines"), plot.margin=unit(c(top, right, bottom, left), "cm"), legend.key.size=unit(leg_size, "line"))
	} else if (ang_le == 90) {
		# set theme
		theme_bw() +
		theme(text=element_text(size=fnt_sz, color="black"), legend.position=leg_pos, legend.text=element_text(size=fnt_sz), legend.key.height=unit(1,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=ang_le, hjust=0.5, vjust=1, color="black", size=fnt_sz), axis.text.y=element_text(color="black", size=fnt_sz), strip.text=element_text(size=fnt_sz), panel.spacing=unit(3, "lines"), plot.margin=unit(c(top, right, bottom, left), "cm"), legend.key.size=unit(leg_size, "line"))
	} else {
		# set theme
		theme_bw() +
		theme(text=element_text(size=fnt_sz, color="black"), legend.position=leg_pos, legend.text=element_text(size=fnt_sz), legend.key.height=unit(1,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(color="black", size=fnt_sz), axis.text.y=element_text(color="black", size=fnt_sz), strip.text=element_text(size=fnt_sz), panel.spacing=unit(3, "lines"), plot.margin=unit(c(top, right, bottom, left), "cm"), legend.key.size=unit(leg_size, "line"))
	}
}
