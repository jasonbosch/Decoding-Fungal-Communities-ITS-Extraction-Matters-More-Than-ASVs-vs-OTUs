# source("X_Analyse_Data_gITS7_ITS4_GABRI_NEW.R")

# export images in case they are needed
export_svg <- function(filename, a_plot, base=F) {

	# if plot is made with ggplot2
	if (base==F) {
		# plot
		svg(paste0(filename, ".svg"), width=22, height=22)
			plot(a_plot)
		dev.off()

	# else with base R
	} else {
		# export plot as svg
		svg(paste0(filename, ".svg"), width=22, height=22)
			a_plot
		dev.off()
	}

	# export plot as Rds
	saveRDS(a_plot, paste0(filename, ".Rds"))
}

# create consistent theme for all ggplots
ggplot_theme <- function(leg_pos="right", ang_le=0, h_just_x=1) {
	# set theme
	theme_bw() +
	theme(text=element_text(size=15, colour="black"), legend.position=leg_pos, legend.text=element_text(size=15), legend.key.height=unit(0.7,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=ang_le, hjust=h_just_x, colour="black", size=15), axis.text.y=element_text(colour="black", size=15), strip.text=element_text(size=15)) +
	theme(legend.title=element_blank())
}

library(groundhog)

groundhog_day <- "2024-04-25"

meta.groundhog(groundhog_day)

###Libraries###

#CRAN libraries
groundhog.library(c("vegan",
                    "stringr",
                    "ggplot2",
                    "cowplot",
                    "ggsignif", #Is this needed?
                    "RColorBrewer"), #Is this needed?
                  groundhog_day)

#BioConductor Libraries
library(ggpubfigs)
library(phyloseq)
library(mirlyn) #IF DOING MULTIPLE RAREFACTION (Technically Github but installed through BioConductor)

###Directory###

root_path <- "/mnt/cinqueg/gabriele/work/microbiology/other_analyses/jason/data/"
final_res <- paste0(root_path, "final_analysis_results_gITS7_ITS4/")
save_img <- paste0(final_res, "figs/")

ifelse(dir.exists(final_res), "all cool", dir.create(final_res))
ifelse(dir.exists(save_img), "all cool", dir.create(save_img))

#combined_alphadiv_df <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/combined_alphadiv_df.Rds"))
#readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/all_phylos.Rds"))
#full_taxonomy <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/full_taxonomy_ready.Rds"))
#physeq_combined <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_combined.Rds"))
#physeq_FWD_ASV <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_FWD_ASV.Rds"))
#readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_FWD_ASV_Rolling.Rds"))
#readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_FWD_OTU.Rds"))
#readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_ITS_ASV.Rds"))
#readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_ITS_ASV_Rolling.Rds"))
#readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_ITS_OTU.Rds"))
#physeq_true <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/physeq_true.Rds"))
#Similarity_scores_2 <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/similarity_scores.Rds"))
#true_community <- readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/true_community_final_clean.Rds"))
#readRDS(paste0(root_path, "final_analysis_results_gITS7_ITS4/true_tax_table.Rds"))

#setwd("~/PostDoc/02_Projects/06_ASV_OTU_Comparison/04_Analysis_Results/gITS7_ITS4/")

###Functions###

#Import QIIME2 taxonomy function
import_qiime2_taxonomy <- function(qiime2_taxonomy) {
  taxonomy_cleaning <- as.data.frame(qiime2_taxonomy[,1])
  taxonomy_cleaning[,1] <- gsub("[a-z]__","",taxonomy_cleaning[,1])
  taxonomy_split <- as.data.frame(matrix(nrow = 1,ncol = 7))
  for (n in 1:nrow(taxonomy_cleaning)) {
    taxonomy_split[n,] <- str_split_fixed(taxonomy_cleaning[n,1],"; |;", n = 7)
  }
  colnames(taxonomy_split) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rownames(taxonomy_split) <- rownames(qiime2_taxonomy)
  tax_table(as.matrix(taxonomy_split))
}

#Filtering
filter_by_library <- function(otu_table, metadata_table) {
  #Seperate by library size
  otu_small <- otu_table[,grepl("small",colnames(otu_table))]
  otu_medium <- otu_table[,grepl("medium",colnames(otu_table))]
  otu_large <- otu_table[,grepl("large",colnames(otu_table))]
  #Remove any ASV with fewer 0.1% of the mean number of reads per sample
  otu_small_final <- otu_table(otu_small[(rowSums(otu_small)>round(mean(colSums(otu_small))*0.001)),],taxa_are_rows = TRUE)
  otu_medium_final <- otu_table(otu_medium[(rowSums(otu_medium)>round(mean(colSums(otu_medium))*0.001)),],taxa_are_rows = TRUE)
  otu_large_final <- otu_table(otu_large[(rowSums(otu_large)>round(mean(colSums(otu_large))*0.001)),],taxa_are_rows = TRUE)
  #Combine the filtered results
  otu_filtered <- merge_phyloseq(otu_small_final,otu_medium_final,otu_large_final)
  cat("Removed ",nrow(otu_table) - nrow(otu_filtered)," taxa.")
  return(otu_filtered)
}

#####LOAD DATA#####

#Shared
metadata_raw <- read.table(paste0(root_path, "Metadata.tsv"), 
                           colClasses = c("character","character","character","character","character","character","character"), header = TRUE, row.names = 1, sep = "\t",)
row.names(metadata_raw) <- gsub("-",".",row.names(metadata_raw))
metadata_raw$Diversity <- factor(metadata_raw$Diversity ,levels = sort(as.numeric(unique(metadata_raw$Diversity))),ordered = T)
metadata_raw$Library_Size <- factor(metadata_raw$Library_Size ,levels = rev(unique(metadata_raw$Library_Size)),ordered = T)
metadata_raw$Replicate <- factor(metadata_raw$Replicate ,levels = unique(metadata_raw$Replicate),ordered = T)
metadata_final <- sample_data(metadata_raw)

###
#FWD ASV Rolling

#Import
table_FWD_ASV_Rolling_raw <- read.table(paste0(root_path, "FWD_ASV_Rolling_gITS7_ITS4/export/OTU.tsv"), 
                                header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
taxonomy_FWD_ASV_Rolling <- import_qiime2_taxonomy(read.table(paste0(root_path, "FWD_ASV_Rolling_gITS7_ITS4/export/taxonomy.tsv"), header = TRUE, row.names = 1, sep = "\t"))

#Filtering
#Remove any ASV with fewer 0.1% of the mean number of reads per sample per library size
table_FWD_ASV_Rolling_final <- filter_by_library(table_FWD_ASV_Rolling_raw,metadata_final)
#407 taxa removed

#Create phyloseq object
physeq_FWD_ASV_Rolling <- phyloseq(table_FWD_ASV_Rolling_final, metadata_final, taxonomy_FWD_ASV_Rolling)
physeq_FWD_ASV_Rolling@sam_data$Processing <- rep("FWD",nrow(physeq_FWD_ASV_Rolling@sam_data))
physeq_FWD_ASV_Rolling@sam_data$Method <- rep("ASV_Rolling",nrow(physeq_FWD_ASV_Rolling@sam_data))
physeq_FWD_ASV_Rolling@sam_data$Treatment <- paste(physeq_FWD_ASV_Rolling@sam_data$Processing,physeq_FWD_ASV_Rolling@sam_data$Method,sep = "_")
sample_names(physeq_FWD_ASV_Rolling) <- gsub("$","_FWD_ASV_Rolling",sample_names(physeq_FWD_ASV_Rolling))

#Check taxonomic assignment
table(physeq_FWD_ASV_Rolling@tax_table[,"Kingdom"])
#4 Eukaryota_kgd_Incertae_sedis, 14 Stramenopila, 6 Unassigned and 21 Viridiplantae
physeq_FWD_ASV_Rolling <- subset_taxa(physeq_FWD_ASV_Rolling, tax_table(physeq_FWD_ASV_Rolling)[,"Kingdom"]=="Fungi")

# export object
saveRDS(physeq_FWD_ASV_Rolling, paste0(final_res, "physeq_FWD_ASV_Rolling.Rds"))

#FWD ASV

#Import
table_FWD_ASV_raw <- read.table(paste0(root_path, "FWD_ASV_gITS7_ITS4/export/OTU.tsv"), 
                                 header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
taxonomy_FWD_ASV <- import_qiime2_taxonomy(read.table(paste0(root_path, "FWD_ASV_gITS7_ITS4/export/taxonomy.tsv"), header = TRUE, row.names = 1, sep = "\t"))

#Filtering
#Remove any ASV with fewer 0.1% of the mean number of reads per sample per library size
table_FWD_ASV_final <- filter_by_library(table_FWD_ASV_raw,metadata_final)
#407 taxa removed

#Create phyloseq object
physeq_FWD_ASV <- phyloseq(table_FWD_ASV_final, metadata_final, taxonomy_FWD_ASV)
physeq_FWD_ASV@sam_data$Processing <- rep("FWD",nrow(physeq_FWD_ASV@sam_data))
physeq_FWD_ASV@sam_data$Method <- rep("ASV",nrow(physeq_FWD_ASV@sam_data))
physeq_FWD_ASV@sam_data$Treatment <- paste(physeq_FWD_ASV@sam_data$Processing,physeq_FWD_ASV@sam_data$Method,sep = "_")
sample_names(physeq_FWD_ASV) <- gsub("$","_FWD_ASV",sample_names(physeq_FWD_ASV))

#Check taxonomic assignment
table(physeq_FWD_ASV@tax_table[,"Kingdom"])
#4 Eukaryota_kgd_Incertae_sedis, 14 Stramenopila, 6 Unassigned and 21 Viridiplantae
physeq_FWD_ASV <- subset_taxa(physeq_FWD_ASV, tax_table(physeq_FWD_ASV)[,"Kingdom"]=="Fungi")

# export object
saveRDS(physeq_FWD_ASV, paste0(final_res, "physeq_FWD_ASV.Rds"))

###
#FWD OTU

#Import
table_FWD_OTU_raw <- read.table(paste0(root_path, "FWD_OTU_gITS7_ITS4/export/OTU-dn-97.tsv"), 
                                header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
taxonomy_FWD_OTU <- import_qiime2_taxonomy(read.table(paste0(root_path, "FWD_OTU_gITS7_ITS4/export/taxonomy.tsv"), header = TRUE, row.names = 1, sep = "\t"))

#Filtering
table_FWD_OTU_final <- filter_by_library(table_FWD_OTU_raw,metadata_final)
#18592 taxa removed

#Create phyloseq object
physeq_FWD_OTU <- phyloseq(table_FWD_OTU_final, metadata_final, taxonomy_FWD_OTU)
physeq_FWD_OTU@sam_data$Processing <- rep("FWD",nrow(physeq_FWD_OTU@sam_data))
physeq_FWD_OTU@sam_data$Method <- rep("OTU",nrow(physeq_FWD_OTU@sam_data))
physeq_FWD_OTU@sam_data$Treatment <- paste(physeq_FWD_OTU@sam_data$Processing,physeq_FWD_OTU@sam_data$Method,sep = "_")
sample_names(physeq_FWD_OTU) <- gsub("$","_FWD_OTU",sample_names(physeq_FWD_OTU))

#Check taxonomic assignment
table(physeq_FWD_OTU@tax_table[,"Kingdom"])
#2 Eukaryota_kgd_Incertae_sedis, 4 Stramenopila, 18 Unassigned and 13 Viridiplantae
physeq_FWD_OTU <- subset_taxa(physeq_FWD_OTU, tax_table(physeq_FWD_OTU)[,"Kingdom"]=="Fungi")

# export object
saveRDS(physeq_FWD_OTU, paste0(final_res, "physeq_FWD_OTU.Rds"))

###
#ITS ASV Rolling

#Import
table_ITS_ASV_Rolling_raw <- read.table(paste0(root_path, "ITS_ASV_Rolling_gITS7_ITS4/export/OTU.tsv"), 
                                header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
taxonomy_ITS_ASV_Rolling <- import_qiime2_taxonomy(read.table(paste0(root_path, "ITS_ASV_Rolling_gITS7_ITS4/export/taxonomy.tsv"), header = TRUE, row.names = 1, sep = "\t"))

#Filtering
table_ITS_ASV_Rolling_final <- filter_by_library(table_ITS_ASV_Rolling_raw,metadata_final)
#226 taxa removed

#Create phyloseq object
physeq_ITS_ASV_Rolling <- phyloseq(table_ITS_ASV_Rolling_final, metadata_final, taxonomy_ITS_ASV_Rolling)
physeq_ITS_ASV_Rolling@sam_data$Processing <- rep("ITS",nrow(physeq_ITS_ASV_Rolling@sam_data))
physeq_ITS_ASV_Rolling@sam_data$Method <- rep("ASV_Rolling",nrow(physeq_ITS_ASV_Rolling@sam_data))
physeq_ITS_ASV_Rolling@sam_data$Treatment <- paste(physeq_ITS_ASV_Rolling@sam_data$Processing,physeq_ITS_ASV_Rolling@sam_data$Method,sep = "_")
sample_names(physeq_ITS_ASV_Rolling) <- gsub("$","_ITS_ASV_Rolling",sample_names(physeq_ITS_ASV_Rolling))

#Check taxonomic assignment
table(physeq_ITS_ASV_Rolling@tax_table[,"Kingdom"])
#7 Eukaryota_kgd_Incertae_sedis, 6 Stramenopila, 7 Unassigned and 2 Viridiplantae
physeq_ITS_ASV_Rolling <- subset_taxa(physeq_ITS_ASV_Rolling, tax_table(physeq_ITS_ASV_Rolling)[,"Kingdom"]=="Fungi")

# export object
saveRDS(physeq_ITS_ASV_Rolling, paste0(final_res, "physeq_ITS_ASV_Rolling.Rds"))

#ITS ASV

#Import
table_ITS_ASV_raw <- read.table(paste0(root_path, "ITS_ASV_gITS7_ITS4/export/OTU.tsv"), 
                                header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
taxonomy_ITS_ASV <- import_qiime2_taxonomy(read.table(paste0(root_path, "ITS_ASV_gITS7_ITS4/export/taxonomy.tsv"), header = TRUE, row.names = 1, sep = "\t"))

#Filtering
table_ITS_ASV_final <- filter_by_library(table_ITS_ASV_raw,metadata_final)
#226 taxa removed

#Create phyloseq object
physeq_ITS_ASV <- phyloseq(table_ITS_ASV_final, metadata_final, taxonomy_ITS_ASV)
physeq_ITS_ASV@sam_data$Processing <- rep("ITS",nrow(physeq_ITS_ASV@sam_data))
physeq_ITS_ASV@sam_data$Method <- rep("ASV",nrow(physeq_ITS_ASV@sam_data))
physeq_ITS_ASV@sam_data$Treatment <- paste(physeq_ITS_ASV@sam_data$Processing,physeq_ITS_ASV@sam_data$Method,sep = "_")
sample_names(physeq_ITS_ASV) <- gsub("$","_ITS_ASV",sample_names(physeq_ITS_ASV))

#Check taxonomic assignment
table(physeq_ITS_ASV@tax_table[,"Kingdom"])
#7 Eukaryota_kgd_Incertae_sedis, 6 Stramenopila, 7 Unassigned and 2 Viridiplantae
physeq_ITS_ASV <- subset_taxa(physeq_ITS_ASV, tax_table(physeq_ITS_ASV)[,"Kingdom"]=="Fungi")

# export object
saveRDS(physeq_ITS_ASV, paste0(final_res, "physeq_ITS_ASV.Rds"))

###
#ITS OTU

#Import
table_ITS_OTU_raw <- read.table(paste0(root_path, "ITS_OTU_gITS7_ITS4/export/OTU-dn-97.tsv"), 
                                header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
taxonomy_ITS_OTU <- import_qiime2_taxonomy(read.table(paste0(root_path, "ITS_OTU_gITS7_ITS4/export/taxonomy.tsv"), header = TRUE, row.names = 1, sep = "\t"))

#Filtering
table_ITS_OTU_final <- filter_by_library(table_ITS_OTU_raw,metadata_final)
#25473 taxa removed

#Create phyloseq object
physeq_ITS_OTU <- phyloseq(table_ITS_OTU_final, metadata_final, taxonomy_ITS_OTU)
physeq_ITS_OTU@sam_data$Processing <- rep("ITS",nrow(physeq_ITS_OTU@sam_data))
physeq_ITS_OTU@sam_data$Method <- rep("OTU",nrow(physeq_ITS_OTU@sam_data))
physeq_ITS_OTU@sam_data$Treatment <- paste(physeq_ITS_OTU@sam_data$Processing,physeq_ITS_OTU@sam_data$Method,sep = "_")
sample_names(physeq_ITS_OTU) <- gsub("$","_ITS_OTU",sample_names(physeq_ITS_OTU))

#Check taxonomic assignment
table(physeq_ITS_OTU@tax_table[,"Kingdom"])
#1 Eukaryota_kgd_Incertae_sedis, 5 Stramenopila, 7 Unassigned and 4 Viridiplantae
physeq_ITS_OTU <- subset_taxa(physeq_ITS_OTU, tax_table(physeq_ITS_OTU)[,"Kingdom"]=="Fungi")

# export object
saveRDS(physeq_ITS_OTU, paste0(final_res, "physeq_ITS_OTU.Rds"))

###
#Bring in the true community
true_community <- data.frame()
true_communities <- as.character(system("ls /mnt/cinqueg/gabriele/work/microbiology/other_analyses/jason/data/ranks_gITS7_ITS4/*ranks.txt",intern = T))

for (FILE in true_communities) {
  tcommunity <- read.table(FILE)
  
  #Collapse reads from different amplicons
  tcommunity$V2 <- gsub("\\|.+","",tcommunity$V2)
  #tcommunity$V2 <- gsub(".dna.toplevel.+","",tcommunity$V2)
  tcommunity <- aggregate(.~V2, tcommunity[,c("V2","V3")],sum)
  
  #Fix underscore names
  tcommunity$V2 <- gsub("^_","",tcommunity$V2)
  
  #Count
  tcounts <- as.data.frame(tcommunity[,2])
  
  #Need to use half of the total reads (e.g. only FWD or merged FWD and REV)
  if (grepl(pattern = "large",x = FILE)) {LIBSIZE <- 0.5*50000} 
  if (grepl(pattern = "medium",x = FILE)) {LIBSIZE <- 0.5*25000}
  if (grepl(pattern = "small",x = FILE)) {LIBSIZE <- 0.5*10000}
  #Convert percentage community to idealised number of reads
  tcounts[,1] <- round(tcounts[,1]/100*LIBSIZE) 
  
  #Replace the percentage with the counts
  tcommunity$V3 <- tcounts[,1]
  colnames(tcommunity) <- c("V2",str_match(FILE,"library_[a-z]+_N_[0-9]+-[0-9]+"))

  #Create the true community
  if (nrow(true_community) == 0) {
    true_community <- tcommunity
  }
  else {
    true_community <- merge.data.frame(true_community,tcommunity, by="V2",all = T)
  }
}

#Replace all NA with 0
for (n in 1:ncol(true_community)) {
  true_community[is.na(true_community[,n]),n] <- 0
}

# RUN THIS INSTEAD
# cut -f2 sh_taxonomy_qiime_ver9_dynamic_all_25.07.2023_dev.txt | sed 's/[a-z]__//g' | tail -n +2 | cut -d";" -f1-7 | sort -u > taxonomy_qiime_gabri.csv

full_tax <- as(read.table("/mnt/cinqueg/gabriele/work/microbiology/other_analyses/jason/data/UNITE_ITS_Classifier/developer/taxonomy_qiime_gabri.csv", sep=";",header = T), "matrix")
colnames(full_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#find the taxonomic matches using the species or genus
true_tax <- as.data.frame(matrix(nrow = nrow(true_community),ncol = 8, dimnames = list(NULL,c("Original","Kingdom","Phylum","Class","Order","Family","Genus","Species")) ))
true_tax$Original <- str_extract(true_community$V2,"^([a-z]|[A-Z])+_([a-z]|[A-Z])+")
for (n in 1:nrow(true_tax)) {
  if (true_tax[n,"Original"]%in%full_tax@.Data[,"Species"]) {
    true_tax[n,2:8] <- full_tax@.Data[full_tax@.Data[,"Species"]==true_tax[n,"Original"]]
  }
  else if (str_extract(str_to_sentence(true_tax[n,"Original"]),"^([a-z]|[A-Z])+")%in%full_tax@.Data[,"Genus"]) {
    true_tax[n,2:8] <- c(unique(full_tax@.Data[full_tax@.Data[,"Genus"]==str_extract(str_to_sentence(true_tax[n,"Original"]),"^([a-z]|[A-Z])+"),1:6]),"")
  }
}
#SOME ARE KNOWN AT FAMILY LEVEL ONLY BUT DON'T GET INFO ADDED. IS THIS IMPORTANT ENOUGH TO FIX? ONLY A FEW TAXA AFFECTED.

#Make sure the true community and true taxa can be linked
row.names(true_community) <- paste("TRUE",row.names(true_community),sep = "")
row.names(true_tax) <- row.names(true_community)

#Clean up
true_tax <- true_tax[,-1]
true_community <- true_community[,-1]
colnames(true_community) <- gsub("-",".",colnames(true_community))

print("exporting true_community table")
saveRDS(true_community, paste0(final_res, "true_community.Rds"))

#Filtering
true_community_final <- filter_by_library(true_community,metadata_final)
#0 taxa removed

saveRDS(true_community_final, paste0(final_res, "true_community_final_clean.Rds"))

#Make QIIME tables
true_tax <- tax_table(as.matrix(true_tax))

#Create phyloseq
physeq_true <- phyloseq(true_community_final, metadata_final, true_tax)
physeq_true@sam_data$Processing <- rep("TRUE",nrow(physeq_true@sam_data))
physeq_true@sam_data$Method <- rep("TRUE",nrow(physeq_true@sam_data))
physeq_true@sam_data$Treatment <- rep("TRUE",nrow(physeq_true@sam_data))
sample_names(physeq_true) <- gsub("$","_T",sample_names(physeq_true))

#2 Stramenopila
physeq_true <- subset_taxa(physeq_true, tax_table(physeq_true)[,"Kingdom"]=="Fungi")

# export object
saveRDS(physeq_true, paste0(final_res, "physeq_true.Rds"))

#Combine the datasets
physeq_combined <- merge_phyloseq(physeq_true,physeq_FWD_ASV_Rolling,physeq_FWD_ASV,physeq_FWD_OTU,physeq_ITS_ASV_Rolling,physeq_ITS_ASV,physeq_ITS_OTU)
table(physeq_combined@tax_table[,"Kingdom"])

# export object
saveRDS(physeq_combined, paste0(final_res, "physeq_combined.Rds"))

# #Filter data to only keep fungal reads
# physeq_combined <- subset_taxa(physeq_combined, tax_table(physeq_combined)[,"Kingdom"]=="Fungi")

#####ANALYSIS#####

###Alpha diversity###
combined_alphadiv_df <- vector()
all_alpha_divs <- list()

###Alpha diversity###
for (libsize in c("Small", "Medium", "Large")) {
	#Rarefy communities
	#Using mirl to repeatedly rarefy to specified library size n times. This is how rarefying is supposed to work.
	#Separate by library size
	mirl_subset <- subset_samples(physeq_combined,Library_Size==libsize)
	#Using mirl to repeatedly rarefy to specified library size n times. This is how rarefying is supposed to work.
	mirl_object <- mirl(mirl_subset, libsize = min(colSums(mirl_subset@otu_table)), rep = 100, mc.cores = 6L, set.seed = 12345)
	# Generates dataframe of alpha-diversity metric from mirl_object
	#Can't use the default because we want different metrics
	alphadiv_small_df <- data.frame()
	for (REP in 1:length(mirl_object)) {
	  cat(paste("Calculating rep ",REP,"\n",sep = ""))
	  md <- sample_data(mirl_object[[REP]])
	  div <- estimate_richness(mirl_object[[REP]],measures = c("Observed","Simpson","Shannon"))
	  final <- cbind(md, div)
	  alphadiv_small_df <- rbind(alphadiv_small_df,final)
	}
	alphadiv_small_df$Treatment <- factor(alphadiv_small_df$Treatment, levels = c("FWD_ASV_Rolling","FWD_ASV","ITS_ASV_Rolling","ITS_ASV","TRUE","FWD_OTU","ITS_OTU"),ordered = T)
	#Join the diversity tables
	combined_alphadiv_df <- rbind.data.frame(combined_alphadiv_df, alphadiv_small_df)
}

# export all alpha divs
saveRDS(combined_alphadiv_df, paste0(final_res, "combined_alphadiv_df.Rds"))

# load data
#combined_alphadiv_df_list <- readRDS(paste0(final_res, "combined_alphadiv_df.Rds"))
#combined_alphadiv_df <- do.call(rbind.data.frame, combined_alphadiv_df_list)

#Draw richness plots
Alpha_observed_fungi <- 
  ggplot(combined_alphadiv_df, aes(x=Diversity, y=Observed)) +
  geom_point (aes(colour = Replicate),position = position_jitterdodge(jitter.width = 0.1),size = 0.3)  +
  geom_boxplot(fill=NA,size=0.3) +
  scale_y_continuous(limits = c(0,850),n.breaks = 6) +
  facet_wrap(~Library_Size*Treatment,ncol = 7) +
  labs(title = "Richness") +
  ggplot_theme(leg_pos="none", ang_le=90)
export_svg(paste0(save_img, "Alpha_observed_fungi"),Alpha_observed_fungi, base=F)

#Draw Simpsons plots
Alpha_simpson_fungi <- 
  ggplot(combined_alphadiv_df, aes(x=Diversity, y=Simpson)) +
  geom_point (aes(colour = Replicate),position = position_jitterdodge(jitter.width = 0.1),size = 0.3)  +
  geom_boxplot(fill=NA,size=0.3) +
  scale_y_continuous(limits = c(0,1), n.breaks = 6) +
  facet_wrap(~Library_Size*Treatment,ncol = 7) +
  labs(title = "Simpson (1-D)") +
  ggplot_theme(leg_pos="none", ang_le=90)
export_svg(paste0(save_img, "Alpha_simpson_fungi"),Alpha_simpson_fungi, base=F)

#Draw Shannon plots
Alpha_shannon_fungi <- 
  ggplot(combined_alphadiv_df, aes(x=Diversity, y=Shannon)) +
  geom_point (aes(colour = Replicate),position = position_jitterdodge(jitter.width = 0.1),size = 0.3)  +
  geom_boxplot(fill=NA,size=0.3) +
  scale_y_continuous(limits = c(0,5.25),n.breaks = 6) +
  facet_wrap(~Library_Size*Treatment,ncol = 7) +
  labs(title = "Shannon") +
  ggplot_theme(leg_pos="none", ang_le=90)
export_svg(paste0(save_img, "Alpha_shannon_fungi"),Alpha_shannon_fungi, base=F)

# #Combined Alpha Diversity
# Alpha <- plot_grid(Alpha_observed_fungi + labs(tag = "A"),
#                    Alpha_simpson_fungi + labs(tag = "B"),
#                    Alpha_shannon_fungi + labs(tag = "C"),
#           ncol = 1)
# #export_svg(paste0(save_img, "Alpha"),Alpha,width = 2*80, height = 2*80, units = "mm")

###Relative abundance###

#Find the differences between what was measured and what was expected.
combined_alphadiv_df$Observed_Diff <- as.numeric(combined_alphadiv_df$Observed)-as.numeric(as.character(combined_alphadiv_df$Diversity))

# get all alpha divs
all_alpha_divs[[libsize]] <- combined_alphadiv_df

Observed_Differences <- 
  ggplot(combined_alphadiv_df, aes(x=Diversity, y=Observed_Diff)) +
  geom_point (aes(colour = Replicate),position = position_jitterdodge(jitter.width = 0.1),size = 0.3)  +
  geom_boxplot(fill=NA,size=0.3) +
  scale_y_continuous(limits = c(-750,350),n.breaks = 6) +
  facet_wrap(~Library_Size*Treatment,ncol = 7) +
  labs(title = "Richness") +
  ylab("Difference between true diversity and observed diversity") +
  ggplot_theme(leg_pos="none", ang_le=90)
export_svg(paste0(save_img, "Observed_Differences"),Observed_Differences, base=F)

###Relative abundance###

#Normalisation
physeq_true_proportion <- transform_sample_counts(physeq_true, function (x) x/sum(x))
physeq_true_proportion@sam_data$Grouping <- paste(physeq_true_proportion@sam_data$Replicate,physeq_true_proportion@sam_data$Library_Size,sep = " ")
physeq_true_proportion@sam_data$Grouping <- factor(physeq_true_proportion@sam_data$Grouping, levels = sort(unique(physeq_true_proportion@sam_data$Grouping))[c(3,2,1,9,8,7,12,11,10,15,14,13,18,17,16,21,20,19,24,23,22,27,26,25,30,29,28,6,5,4)], ordered = T)

#Phylum level (Anything less than 1% is merged together)
physeq_phylum <- tax_glom(physeq_true_proportion, taxrank = "Phylum",NArm = FALSE)
physeq_phylum <- merge_taxa(physeq_phylum,taxa_names(filter_taxa(physeq_phylum, function(x) mean(x) < 0.01, TRUE)))
phylum_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_phylum)[,"Phylum"])))
names(phylum_colours) <- sort(unique(tax_table(physeq_phylum)[,"Phylum"]))
plot_true_phylum <-
  plot_bar(physeq_phylum, fill = "Phylum", x = "Grouping") +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x= element_blank(), y = "Abundance") +
  facet_wrap(~Diversity,ncol = 1, scales = "free_x") +
  scale_fill_manual(values = phylum_colours, labels = c(sort(unique(physeq_phylum@tax_table[,"Phylum"]),decreasing = F),"Rare_Taxa"), guide = guide_legend(ncol = 1)) +
  ggplot_theme(leg_pos="right", ang_le=45)
export_svg(paste0(save_img, "plot_true_phylum"),plot_true_phylum, base=F)

#Class level (Anything less than 1% is merged together)
physeq_class <- tax_glom(physeq_true_proportion, taxrank = "Class",NArm = FALSE)
physeq_class <- merge_taxa(physeq_class,taxa_names(filter_taxa(physeq_class, function(x) mean(x) < 0.01, TRUE)))
class_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_class)[,"Class"])))
names(class_colours) <- sort(unique(tax_table(physeq_class)[,"Class"]))
plot_true_class <-
  plot_bar(physeq_class, fill = "Class", x = "Grouping") +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x= element_blank(), y = "Abundance") +
  facet_wrap(~Diversity,ncol = 1, scales = "free_x") +
  scale_fill_manual(values = class_colours, labels = c(sort(unique(physeq_class@tax_table[,"Class"]),decreasing = F),"Rare_Taxa"), guide = guide_legend(ncol = 1)) +
  ggplot_theme(leg_pos="right", ang_le=45)
export_svg(paste0(save_img, "plot_true_class"),plot_true_class, base=F)

#Order level (Anything less than 1% is merged together)
physeq_order <- tax_glom(physeq_true_proportion, taxrank = "Order",NArm = FALSE)
physeq_order <- merge_taxa(physeq_order,taxa_names(filter_taxa(physeq_order, function(x) mean(x) < 0.01, TRUE)))
order_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_order)[,"Order"])))
names(order_colours) <- sort(unique(tax_table(physeq_order)[,"Order"]))
plot_true_order <-
  plot_bar(physeq_order, fill = "Order", x = "Grouping") +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x= element_blank(), y = "Abundance") +
  facet_wrap(~Diversity,ncol = 1, scales = "free_x") +
  scale_fill_manual(values = order_colours, labels = c(sort(unique(physeq_order@tax_table[,"Order"]),decreasing = F),"Rare_Taxa"), guide = guide_legend(ncol = 1)) +
  ggplot_theme(leg_pos="right", ang_le=45)
export_svg(paste0(save_img, "plot_true_order"),plot_true_order, base=F)

#Family level (Anything less than 1% is merged together)
physeq_family <- tax_glom(physeq_true_proportion, taxrank = "Family",NArm = FALSE)
physeq_family <- merge_taxa(physeq_family,taxa_names(filter_taxa(physeq_family, function(x) mean(x) < 0.01, TRUE)))
family_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_family)[,"Family"])))
names(family_colours) <- sort(unique(tax_table(physeq_family)[,"Family"]))
plot_true_family <-
  plot_bar(physeq_family, fill = "Family", x = "Grouping") +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x= element_blank(), y = "Abundance") +
  facet_wrap(~Diversity,ncol = 1, scales = "free_x") +
  scale_fill_manual(values = family_colours, labels = c(sort(unique(physeq_family@tax_table[,"Family"]),decreasing = F),"Rare_Taxa"), guide = guide_legend(ncol = 1)) +
  ggplot_theme(leg_pos="right", ang_le=45)
export_svg(paste0(save_img, "plot_true_family"),plot_true_family, base=F)

#Genus level (Anything less than 1% is merged together)
EvaluateByDiversity <- function(df) {
  res_table <- as.data.frame(matrix(nrow = nrow(df)),row.names = row.names(df))
  for (DIVERSITY in unique(physeq_combined@sam_data$Diversity)) {
    curr <- df[,sample_names(df)%in%row.names(physeq_combined@sam_data[physeq_combined@sam_data$Diversity==DIVERSITY,])]
    res_table[,DIVERSITY] <- rowMeans(curr)
  }
  res_table <- res_table[,-1]
  as.vector(apply(res_table, 1, FUN = max))
}

#Genus level (Anything less than 1% is merged together)
physeq_genus <- tax_glom(physeq_true_proportion, taxrank = "Genus",NArm = FALSE)
physeq_genus <- merge_taxa(physeq_genus,taxa_names(prune_taxa(EvaluateByDiversity(physeq_genus@otu_table) < 0.01,physeq_genus)))
genus_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_genus)[,"Genus"])))
names(genus_colours) <- sort(unique(tax_table(physeq_genus)[,"Genus"]))
plot_true_genus <-
  plot_bar(physeq_genus, fill = "Genus", x = "Grouping") +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x= element_blank(), y = "Abundance") +
  facet_wrap(~Diversity,ncol = 1, scales = "free_x") +
  scale_fill_manual(values = genus_colours, labels = c(sort(unique(physeq_genus@tax_table[,"Genus"]),decreasing = F),"Rare_Taxa"), guide = guide_legend(ncol = 1)) +
  ggplot_theme(leg_pos="right", ang_le=45)
export_svg(paste0(save_img, "plot_true_genus"),plot_true_genus, base=F)

#New loop
Similarity_scores_2 <- as.data.frame(matrix(ncol = 12, dimnames = list(NULL,c("Diversity","Library_Size","Replicate","Processing","Method","Treatment","Level","Sensitivity","Precision","Aitchison","Hellinger","BrayCurtis"))))
for (TAXLEVEL in c("Kingdom","Phylum","Class","Order","Family","Genus")) {
  Physeq_TAXLEVEL <- tax_glom(physeq_combined, taxrank = TAXLEVEL,NArm = FALSE)
  
  Physeq_TAXLEVEL_AITCHISON <- vegdist(t(Physeq_TAXLEVEL@otu_table),method = "robust.aitchison")
  Physeq_TAXLEVEL_HELLINGER <- distance(otu_table(decostand(Physeq_TAXLEVEL@otu_table,"hellinger",MARGIN = 2),taxa_are_rows = TRUE),method = "euclidean")
  Physeq_TAXLEVEL_BRAYCURTIS <- distance(Physeq_TAXLEVEL,method = "bray")
  
  for (TREATMENT in c("FWD_ASV_Rolling","FWD_ASV","FWD_OTU","ITS_ASV_Rolling","ITS_ASV","ITS_OTU")) {
    for (LIB_SIZE in c("Small","Medium","Large")) {
      for (DIVERSITY in unique(Physeq_TAXLEVEL@sam_data$Diversity)) {
        cat(paste("Calculating",TAXLEVEL,TREATMENT,LIB_SIZE,DIVERSITY,"\n",sep = " "))
        
        SAMPLES <- rownames(Physeq_TAXLEVEL@sam_data[Physeq_TAXLEVEL@sam_data$Treatment==TREATMENT&Physeq_TAXLEVEL@sam_data$Library_Size==LIB_SIZE&Physeq_TAXLEVEL@sam_data$Diversity==DIVERSITY])
        
        for (S in SAMPLES) {
          PT <- row.names(Physeq_TAXLEVEL@otu_table[Physeq_TAXLEVEL@otu_table[,gsub(paste(TREATMENT,sep = ""),"T",S)]>=1,]) #Present in true
          PS <- row.names(Physeq_TAXLEVEL@otu_table[Physeq_TAXLEVEL@otu_table[,S]>=1,]) #Present in sample
          TP <- sum(PS%in%PT)
          FN <- sum(!PT%in%PS)
          FP <- sum(!PS%in%PT)
          D_AITCH <- as.matrix(Physeq_TAXLEVEL_AITCHISON)[S,gsub(paste(TREATMENT,sep = ""),"T",S)]
          D_HELL <- as.matrix(Physeq_TAXLEVEL_HELLINGER)[S,gsub(paste(TREATMENT,sep = ""),"T",S)]
          D_BRAY <- as.matrix(Physeq_TAXLEVEL_BRAYCURTIS)[S,gsub(paste(TREATMENT,sep = ""),"T",S)]
          
          Similarity_scores_2[nrow(Similarity_scores_2)+1,1:12]<- c(DIVERSITY,
                                                                    LIB_SIZE,
                                                                    str_extract(str_extract(S,".[0-9]+_"),"[0-9]+"),
                                                                    str_extract(TREATMENT,"^[A-Z]{3}"),
                                                                    sub("^[A-Z]{3}_","",TREATMENT),
                                                                    TREATMENT,
                                                                    TAXLEVEL,
                                                                    TP/(TP+FN),
                                                                    TP/(TP+FP),
                                                                    D_AITCH,
                                                                    D_HELL,
                                                                    D_BRAY)
        }
      }
    }
  }
}
Similarity_scores_2 <- Similarity_scores_2[-1,]
Similarity_scores_2$Level <- factor(Similarity_scores_2$Level, levels = c("Kingdom","Phylum","Class","Order","Family","Genus"), ordered = T)
Similarity_scores_2$Library_Size <- factor(Similarity_scores_2$Library_Size, levels = c("Small","Medium","Large"), ordered = T)
Similarity_scores_2$Diversity <- factor(Similarity_scores_2$Diversity ,levels = sort(as.numeric(unique(Similarity_scores_2$Diversity))),ordered = T)
Similarity_scores_2$Sensitivity <- as.numeric(Similarity_scores_2$Sensitivity)
Similarity_scores_2$Precision <- as.numeric(Similarity_scores_2$Precision)
Similarity_scores_2$Aitchison <- as.numeric(Similarity_scores_2$Aitchison)
Similarity_scores_2$Hellinger <- as.numeric(Similarity_scores_2$Hellinger)
Similarity_scores_2$BrayCurtis <- as.numeric(Similarity_scores_2$BrayCurtis)

saveRDS(Similarity_scores_2, paste0(root_path, "final_analysis_results_gITS7_ITS4/similarity_scores.Rds"))

Similarity_scores_2$Diversity <- factor(Similarity_scores_2$Diversity, levels=unique(Similarity_scores_2$Diversity))
Similarity_scores_2$Library_Size <- factor(Similarity_scores_2$Library_Size, levels=unique(Similarity_scores_2$Library_Size))
Similarity_scores_2$Treatment <- factor(Similarity_scores_2$Treatment, levels=unique(Similarity_scores_2$Treatment))
Similarity_scores_2$Level <- factor(Similarity_scores_2$Level, levels=unique(Similarity_scores_2$Level))

#Can we see which is better at representing the true community?

#Find the differences between what was measured and what was expected.

# combined_alphadiv_df <- readRDS(paste0(final_res, "combined_alphadiv_df.Rds"))

Ranking_observed <- as.data.frame(aggregate(.~Diversity*Library_Size*Treatment, combined_alphadiv_df[,c("Diversity","Library_Size","Treatment","Observed_Diff")],mean))
Ranking_observed <- subset(Ranking_observed,Ranking_observed$Treatment!="TRUE")
Ranking_observed$Comparison <- paste(Ranking_observed$Diversity,Ranking_observed$Library_Size)
for (COMPARISON in unique(Ranking_observed$Comparison)) {
  Ranking_observed[Ranking_observed$Comparison==COMPARISON,"Rank"] <- rank(abs(Ranking_observed[Ranking_observed$Comparison==COMPARISON,"Observed_Diff"]),ties.method = "min")
}
ranks_obs_mean <- aggregate(.~Treatment, Ranking_observed[,c("Treatment","Rank")],mean)
ranks_obs_sum <- aggregate(.~Treatment, Ranking_observed[,c("Treatment","Rank")],sum)

Ranking_similarity <- as.data.frame(aggregate(.~Diversity*Library_Size*Treatment*Level, Similarity_scores_2[,c("Diversity","Library_Size","Treatment","Level","Sensitivity","Precision","Aitchison","Hellinger","BrayCurtis")],mean))
Ranking_similarity$Comparison <- paste(Ranking_similarity$Diversity,Ranking_similarity$Library_Size,Ranking_similarity$Level)
for (COMPARISON in unique(Ranking_similarity$Comparison)) {

  # COMPARISON <- unique(Ranking_similarity$Comparison)[1]

  Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Rank_sen"] <- rank(1-Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Sensitivity"],ties.method = "min")
  Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Rank_pre"] <- rank(1-Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Precision"],ties.method = "min")
  Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Rank_ait"] <- rank(Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Aitchison"],ties.method = "min")
  Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Rank_hel"] <- rank(Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Hellinger"],ties.method = "min")
  Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"Rank_bra"] <- rank(Ranking_similarity[Ranking_similarity$Comparison==COMPARISON,"BrayCurtis"],ties.method = "min")
}
ranks_sen_mean <- aggregate(.~Treatment, Ranking_similarity[,c("Treatment","Rank_sen","Rank_pre","Rank_ait","Rank_hel","Rank_bra")],mean)
ranks_sen_sum <- aggregate(.~Treatment, Ranking_similarity[,c("Treatment","Rank_sen","Rank_pre","Rank_ait","Rank_hel","Rank_bra")],sum)

Rankings_mean <- merge(ranks_obs_mean, ranks_sen_mean)
Rankings_mean$Average <- rowMeans(Rankings_mean[,c("Rank","Rank_sen","Rank_pre","Rank_ait","Rank_hel","Rank_bra")])
Rankings_mean$Final <- rank(Rankings_mean$Average)

Rankings_sum <- merge(ranks_obs_sum, ranks_sen_sum)
Rankings_sum$Sum <- rowSums(Rankings_sum[,c("Rank","Rank_sen","Rank_pre","Rank_ait","Rank_hel","Rank_bra")])
colnames(Rankings_sum) <- c("Treatment",paste(colnames(Rankings_sum),"sum",sep = "_")[-1])

Final_Rankings <- merge(Rankings_mean,Rankings_sum)
colnames(Final_Rankings) <- c("Treatment","Observed","Sensitivity","Precision","Aitchison","Hellinger","BrayCurtis","Average","Final","Observed_sum","Sensitivity_Sum","Precision_Sum","Aitchison_Sum","Hellinger_Sum","BrayCurtis_Sum","Sum")
Final_Rankings <- Final_Rankings[,c("Treatment","Observed","Observed_sum","Sensitivity","Sensitivity_Sum","Precision","Precision_Sum","Aitchison","Aitchison_Sum","Hellinger","Hellinger_Sum","BrayCurtis","BrayCurtis_Sum","Average","Sum","Final")]
write.csv(Final_Rankings, "Final_Rankings.csv",row.names = F,quote = T)

#Set colours#
Treatment_Colours <- friendly_pal("vibrant_seven",6)
names(Treatment_Colours) <- unique(Similarity_scores_2$Treatment)

#Plot the scoring
plot_sensitivity_2 <- 
  ggplot(Similarity_scores_2,aes(x=Level, y = Sensitivity, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x="", title = "Sensitivity (TP/(TP+FN))") + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = Treatment_Colours) +
  facet_grid(rows = vars(Diversity), cols = vars(Library_Size)) +
  ggplot_theme(leg_pos="bottom", ang_le=45)
export_svg(paste0(save_img, "plot_sensitivity_2"),plot_sensitivity_2, base=F)

plot_precision_2 <- 
  ggplot(Similarity_scores_2,aes(x=Level, y = Precision, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x="", title = "Precision (TP/(TP+FP))") + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = Treatment_Colours) +
  facet_grid(rows = vars(Diversity), cols = vars(Library_Size)) +
  ggplot_theme(leg_pos="bottom", ang_le=45)
export_svg(paste0(save_img, "plot_precision_2"),plot_precision_2, base=F)

plot_aitchison_2 <- 
  ggplot(Similarity_scores_2,aes(x=Level, y = Aitchison, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x="", title = "Aitchison distance") + 
  scale_y_continuous(n.breaks = 5) +
  coord_cartesian(ylim = c(0,max(Similarity_scores_2$Aitchison))) +
  scale_color_manual(values = Treatment_Colours) +
  facet_grid(rows = vars(Diversity), cols = vars(Library_Size)) +
  ggplot_theme(leg_pos="bottom", ang_le=45)
export_svg(paste0(save_img, "plot_aitchison_2"),plot_aitchison_2, base=F)

plot_hellinger_2 <- 
  ggplot(Similarity_scores_2,aes(x=Level, y = Hellinger, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x="", title = "Hellinger distance") + 
  scale_y_continuous(n.breaks = 5) +
  coord_cartesian(ylim = c(0,max(Similarity_scores_2$Hellinger))) +
  scale_color_manual(values = Treatment_Colours) +
  facet_grid(rows = vars(Diversity), cols = vars(Library_Size)) +
  ggplot_theme(leg_pos="bottom", ang_le=45)
export_svg(paste0(save_img, "plot_hellinger_2"),plot_hellinger_2, base=F)

plot_bray_2 <- 
  ggplot(Similarity_scores_2,aes(x=Level, y = BrayCurtis, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x="", title = "Bray Curtis dissimilarity") + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = Treatment_Colours) +
  facet_grid(rows = vars(Diversity), cols = vars(Library_Size)) +
  ggplot_theme(leg_pos="bottom", ang_le=45)
export_svg(paste0(save_img, "plot_bray_2"),plot_bray_2, base=F)

#Combine them but only showing Diversity = 800 for space reasons

plot_sensitivity_sub_2 <- 
  ggplot(Similarity_scores_2[Similarity_scores_2$Diversity==800,],aes(x=Level, y = Sensitivity, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x = element_blank()) + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = Treatment_Colours) +
  facet_grid(~Library_Size) +
  ggplot_theme(leg_pos="bottom")

plot_precision_sub_2 <- 
  ggplot(Similarity_scores_2[Similarity_scores_2$Diversity==800,],aes(x=Level, y = Precision, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x = element_blank()) + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = Treatment_Colours) +
  facet_grid(~Library_Size) +
  ggplot_theme(leg_pos="bottom")

plot_aitchison_sub_2 <- 
  ggplot(Similarity_scores_2[Similarity_scores_2$Diversity==800,],aes(x=Level, y = Aitchison, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x = element_blank()) + 
  coord_cartesian(ylim=c(0, max(Similarity_scores_2$Aitchison))) +
  scale_y_continuous(n.breaks = 6) +
  scale_color_manual(values = Treatment_Colours) +
  facet_wrap(~Library_Size,ncol = 3) +
  ggplot_theme(leg_pos="bottom")

plot_hellinger_sub_2 <- 
  ggplot(Similarity_scores_2[Similarity_scores_2$Diversity==800,],aes(x=Level, y = Hellinger, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x = element_blank()) + 
  coord_cartesian(ylim=c(0, max(Similarity_scores_2$Hellinger))) +
  scale_y_continuous(n.breaks = 6) +
  scale_color_manual(values = Treatment_Colours) +
  facet_wrap(~Library_Size,ncol = 3) +
  ggplot_theme(leg_pos="bottom")

plot_bray_sub_2 <- 
  ggplot(Similarity_scores_2[Similarity_scores_2$Diversity==800,],aes(x=Level, y = BrayCurtis, group = Treatment, colour = Treatment)) +
  stat_summary(position=position_dodge2(width=0.3), fun = "mean", geom="line") +
  stat_summary(position=position_dodge2(width=0.3), fun.data = "mean_sdl", fun.args = list(mult = 1), geom="errorbar", width = 0.3) +
  labs(x = element_blank()) + 
  coord_cartesian(ylim=c(0, 1)) +
  scale_y_continuous(n.breaks = 6) +
  scale_color_manual(values = Treatment_Colours) +
  facet_wrap(~Library_Size, ncol = 3) +
  ggplot_theme(leg_pos="bottom", ang_le=45)

#Combine them

#Sensitivity_legend_2<- get_legend(plot_sensitivity_sub_2)
Sensitivity_legend_2<- get_plot_component(plot_sensitivity_sub_2, 'guide-box-bottom', return_all = TRUE)
Sensitivity_plots_2 <-plot_grid(plot_grid(plot_sensitivity_sub_2 + theme(legend.position = "none"),
                                        plot_precision_sub_2 + theme(legend.position = "none"),
                                        plot_aitchison_sub_2 + theme(legend.position = "none"),
                                        plot_hellinger_sub_2 + theme(legend.position = "none"),
                                        plot_bray_sub_2 + theme(legend.position = "none"),
                                        ncol = 1, rel_heights = c(1,0.9,0.9,0.9,1.1)),
                                        Sensitivity_legend_2,
                                        ncol = 1, rel_heights = c(1,0.05))
export_svg(paste0(save_img, "Sensitivity_plots_2"),Sensitivity_plots_2, base=F)

