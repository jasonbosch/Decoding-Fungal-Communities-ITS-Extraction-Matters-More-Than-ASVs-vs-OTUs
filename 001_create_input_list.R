# download this file
# https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/fungi/species_EnsemblFungi.txt

# add a first new column with numbers from 1-end so readcsv does not have issues in loading

# load libraries
library("stringr") 

# example URLS from the ensembl database
# https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/fungi/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz

# load the file with all the info
all_the_info <- read.csv("species_EnsemblFungi.txt", sep="\t", header=T, stringsAsFactors=F)

# starting string
beginning <- "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/fungi/fasta/"
# middle string
middle <- "/dna/"
# a dot
a_dot <- "."
# ending string
ending <- ".dna.toplevel.fa.gz"

# get info from the table containing all database info
f_one <- paste0(str_sub(sapply(strsplit(all_the_info$core_db, "core"), "[[", 1), end=-2), "/")
f_two <- all_the_info$species
f_two_mod <- paste0(toupper(str_sub(f_two, 1, 1)), substring(f_two, 2))
f_three <- paste0(".", all_the_info$assembly)

# put all together
all_downloadable_species <- cbind.data.frame(beginning, f_one, f_two, middle, f_two_mod, f_three, ending)

# export table to be used though a bash script
write.table(all_downloadable_species, "download_me_with_spaces.txt", sep="", col.names=F, row.names=F, quote=FALSE)

#library("Biostrings")

## gITS7 (5′-GTGAATCATCGAATCTTTG-3′)
#forward <- DNAString("GTGAATCATCGAATCTTTG")

#print(paste("Forward ", as.character(forward)))
#print(paste("Reverse forward ", as.character(reverse(forward))))
#print(paste("Reverse-complement forward ", as.character(reverseComplement(forward))))

## ITS4 (5′-TCCTCCGCTTATTGATATGC- 3′)
#reverse <- DNAString("TCCTCCGCTTATTGATATGC")

#print(paste("Reverse ", as.character(reverse)))
#print(paste("Reverse reverse ", as.character(reverse(reverse))))
#print(paste("Reverse-complement reverse ", as.character(reverseComplement(reverse))))

