# source("substitute_scores.R")

# get scores from one fastq file

# head -n 200003 PRI-LANY-Jana-AD006_S1_L001_R1_001.fastq | grep -A1 "+" | grep -v "+" > hundred_thousand_scores.txt

# this files is used to substitute scores in existing fastq files

# list fastq files
fast_qs <- list.files()[grep("fastq", list.files())]

for (a_lib in fast_qs) {
	
	# a_lib <- fast_qs[1]
	
	print(paste0("analysing ", a_lib))
	
	# load a library
	a_library <- read.table(a_lib, header=F, sep="\n", stringsAsFactors=F)	

	print("lib loaded")

	# load scores
	scores <- read.table("/mnt/DATA01/ASV_OTU_Project/06_Miscellaneous/fastq_scores/hundred_thousands_scores.txt", header=F, sep="\n", stringsAsFactors=F)

	print("scores loaded")

	# loop through all sequences
	for (a_seq in seq(2, nrow(a_library), by=4)) {
	
		# a_seq <- seq(2, nrow(a_library), by=4)[3]

		# choose a random sequence to generate a score
		a_rand_seq <- scores[sample(100000, 1, replace=F)[1], 1]
		
		# length of current sequence
		current_seq_len <- nchar(a_library[a_seq, 1])
		
		# is score length enough to cover the sequence?
		if (nchar(a_rand_seq) >= current_seq_len) {
		
			# if they have same length
			if (nchar(a_rand_seq) == current_seq_len) {
	
				# update the score
				new_score <- a_rand_seq
			} else { # otherwise

				# cut the end of the score
				new_score <- substr(a_rand_seq, 1, current_seq_len)
			}
		} else { # otherwise get another random score
		
			# get score
			another_rand_seq <- scores[sample(100000, 1, replace=F), 1]
			
			# create a chimera score attacching the middle of one score
			# in the middle of the other score. this assumes that, to
			# keep the first bases with high quality score and last
			# bases with low quality, best would be to substitute the
			# central part with another central part of mixed quality.
			# to achieve this:
			# get missing number of scores
			missing_scores <- current_seq_len - nchar(a_rand_seq)
			
			# get missing from another_rand_seq
			half_of_another <- round(nchar(another_rand_seq)/2)
			
			# test if missing_scores is odd
			if (missing_scores %% 2 != 0) {
			
				# get the upper half and round to lower, if odd
				upper_half <- half_of_another + ceiling(missing_scores/2)
				# get the lower half and round to upper, if odd
				lower_half <- half_of_another - floor(missing_scores/2)
				# finally, get the missing score
				# remember the -1 to get a score that is long as
				# it should be
				middle_score <- substr(another_rand_seq, lower_half, upper_half-1)

			} else {
				# otherwise is even
			
				# get the upper half and round to lower, if odd
				upper_half <- half_of_another + ceiling(missing_scores/2)
				# get the lower half and round to upper, if odd
				lower_half <- (half_of_another + 1) - (floor(missing_scores/2))
				# finally, get the missing score
				# here, don't remove the -1
				middle_score <- substr(another_rand_seq, lower_half, upper_half)
			}
			
			# paste the score in the middle of the short sequence
			new_score <- paste0(substr(a_rand_seq, 1, floor(nchar(a_rand_seq)/2)), middle_score, substr(a_rand_seq, floor(nchar(a_rand_seq)/2)+1, nchar(a_rand_seq)))
		}
		
		# finally, substitute the new score
		a_library[a_seq+2, 1] <- new_score
		
	}
	
	# finally, export the updated library
	write.table(a_library, paste0(strsplit(a_lib, "\\.")[[1]][1], "_new_scores.fastq"), row.names=F, col.names=F, quote=F)
}

