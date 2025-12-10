# load library
library("stringr")
library("Biostrings")

# print out the script is running
print("substituting primers")

# get all files
all_files <- list.files("./")
# then filter by fasta
file_names <- all_files[grep("fasta", all_files)]

# get length of primers
fwd_primer <- "AATTYGAHTCAACRCGGG"
rev_primer <- "CGACRGGMGGTGTGBACA"

fwd_length <- str_length(fwd_primer)
rev_length <- str_length(rev_primer)

# loop through file names
for (flnm in file_names) {

#	flnm <- "Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa_reads.fasta"
	
	# read sequence file
	organism <- read.csv(flnm, header=F, stringsAsFactors=F)
	# read trimming file
	trim <- read.csv("./trim_here.csv", header=F, stringsAsFactors=F)
	
	# set colnames form trim file
	colnames(trim) <- c("start", "stop", "primer")
	
	# loop through sequences and substitute
	for (sqnc in seq(2, nrow(organism), by=2)) {
		
		# get sequence original start and stop to find which primer
		# was actually found by primersearch
		orig_start <- as.integer(strsplit(organism[sqnc-1, 1], "\\|")[[1]][2])
		orig_stop <- as.integer(strsplit(organism[sqnc-1, 1], "\\|")[[1]][3])
		# get primer from trim table
		primer <- trim$primer[which(trim$start==orig_start & trim$stop==orig_stop)]

		# get sequence length, start
		start_trim <- fwd_length + 1
		# get stop
		stop_trim <- str_length(organism[sqnc, 1]) - rev_length		

		# trimmed sequence
		trimmed <- str_sub(organism[sqnc, 1], start=start_trim, end=stop_trim)

		# which primer was found? if forward then
		if (primer == fwd_primer) {
			
			# substitute forward and reverse
			brand_new <- paste0(fwd_primer, trimmed, as.character(reverseComplement(DNAString(rev_primer))))
			
			# finally, substitute the brand new in the organism file
			organism[sqnc, 1] <- brand_new
		} else {
			# if reverse then
			
			# substitute forward and reverse
			brand_new <- paste0(rev_primer, trimmed, as.character(reverseComplement(DNAString(fwd_primer))))
			
			# finally, substitute the brand new in the organism file
			organism[sqnc, 1] <- brand_new
		}
	}
}

# export the new file
write.table(organism, paste0(gsub(".fasta", "_wobble_primers.fa", flnm)), row.names=F, col.names=F, quote=F)

# print out that the script ended
print("005 ends here")

