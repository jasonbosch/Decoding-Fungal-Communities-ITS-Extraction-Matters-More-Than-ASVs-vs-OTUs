#setwd("~/PostDoc/02_Projects/06_ASV_OTU_Comparison/06_Miscellaneous/Test_Run/Gabrielle_Test/")

library("stringr")

#Set arguments
args <- commandArgs(trailingOnly=TRUE) #args <- c("OUT=result","DUPS=TRUE","DIST=power","RICH=4","LIBS=5")
if (sum(grepl(pattern = "OUT=",args))) {
  outname <- sub("OUT=","",str_extract(args[grepl(pattern = "OUT=",args)],"OUT=[^[:blank:]]+"))
} else {
  outname <- "grinder_abundance_file" #OUT= filename
}
if (sum(grepl(pattern = "DUPS=",args))) {
  duplicate_species <- as.logical(str_extract(str_extract(args[grepl(pattern = "DUPS=",args)],"DUPS=[:alpha:]+"),"[:alpha:]+$"))
} else {
  duplicate_species <- FALSE #DUPS= TRUE or FALSE
}
if (sum(grepl(pattern = "DIST=",args))) {
  distribution <- str_extract(str_extract(args[grepl(pattern = "DIST=",args)],"DIST=[:alpha:]+"),"[:alpha:]+$")
} else {
  distribution <- "power" #DIST= power or uniform
}
if (sum(grepl(pattern = "RICH=",args))) {
  richness <- as.numeric(str_extract(str_extract(args[grepl(pattern = "RICH=",args)],"RICH=[0-9]+"),"[0-9]+$"))
} else {
  richness <- 150 #RICH= number
}
if (sum(grepl(pattern = "LIBS=",args))) {
  libraries <- as.numeric(str_extract(str_extract(args[grepl(pattern = "LIBS=",args)],"LIBS=[0-9]+"),"[0-9]+$"))
} else {
  libraries <- 5 #LIBS= number
}
if (sum(grepl(pattern = "SEED=",args))) {
  seed <- as.numeric(str_extract(str_extract(args[grepl(pattern = "SEED=",args)],"SEED=[0-9]+"),"[0-9]+$"))
} else {
  seed <- 123456 #SEED= number
}

#Uniform abundances 
uniform_abundances <- rep((1/richness),richness)

#powerlaw abundances
x <- 1:richness
y <- (x)^-1.08
powerlaw_abundances <- y/sum(y)

#Get species and assign abundances
amplicon_files <- as.data.frame(system("ls ./*.fa",intern = T))
colnames(amplicon_files) <- "file"
amplicon_files$Binomial <- gsub("./","",amplicon_files$file)
amplicon_files$Binomial <- str_extract(amplicon_files$Binomial, "[:alpha:]+_[:alpha:]+")

#If no duplicate species
if (!duplicate_species) {
  amplicon_files <- amplicon_files[!duplicated(amplicon_files$Binomial),]
}

#Check richness is compatible
if (nrow(amplicon_files)<richness) {
  stop("Richness exceeds the number of species.")
}

#File to store the communities
grinder_communities <- list()

#Do for each library
for (LIB in 1:libraries) {
  #Select the files that will be used and assign abundance
  set.seed(seed+LIB)
  selected_species <- amplicon_files[sample(nrow(amplicon_files),richness),]
  if (distribution == "power") {
    selected_species$RA <- powerlaw_abundances
  }
  if (distribution == "uniform") {
    selected_species$RA <- uniform_abundances
  }
  
  #Create an abundance file which is the fasta names and their abundances (even distribution of species abundance).
  library_abundance_file <- data.frame()
  for (n in 1:nrow(selected_species)) {
    amplicons <- as.data.frame(system(paste(r"(grep \> )",selected_species[n,"file"],sep = ""),intern = T))
    amplicons$RA <- selected_species[n,"RA"]/nrow(amplicons)
    library_abundance_file <- rbind(library_abundance_file,amplicons)
  }
  library_abundance_file[,1] <- gsub(">","",library_abundance_file[,1])
  colnames(library_abundance_file) <- c("Amplicon",paste("RA",LIB,sep = "_"))

  #Add it to the list
  grinder_communities[LIB] <- list(library_abundance_file)
}

#Collapse the list into a single file
reduce_merge <- function(x, y) {
  merge.data.frame(x,y,all = TRUE)
}
grinder_abundance_file <- Reduce(reduce_merge, grinder_communities[1:libraries])
grinder_abundance_file[is.na(grinder_abundance_file)] <- 0

#Save the abundance file
write.table(grinder_abundance_file, file = paste(outname,"N",richness,sep = "_"), quote = FALSE, col.names = FALSE, row.names = FALSE)




