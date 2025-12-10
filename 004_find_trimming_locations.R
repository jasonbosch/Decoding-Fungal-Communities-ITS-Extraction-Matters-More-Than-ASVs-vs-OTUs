# print out the script is running
print("find trimming locations")

# reorganising starts ans stops
# load positions
positions_raw <- read.csv("trimming_positions.txt", header=F, sep=",", stringsAsFactors=F)
colnames(positions_raw) <- c("start", "length")

# positions ok
positions <- cbind.data.frame(primer=gsub("\t", "", sapply(strsplit(positions_raw$start, " "), "[[", 1)), start=as.numeric(sapply(strsplit(positions_raw$start, " "), "[[", 6)), length=as.numeric(sapply(strsplit(positions_raw$length, " "), "[[", 3)))

# compute stop positions depending on where the fwd primer was found
positions$stop <- positions$start+(positions$length-1)

######################################################## switch start with stop if necessary

# reorganise start and stop
alreadyok <- positions[which(positions$start < positions$stop), ][, c("start", "stop", "primer")]
colnames(alreadyok) <- c("start", "stop", "primer")

# if start comes after stop
adjusted <- positions[which(positions$start > positions$stop), ][, c("start", "stop", "primer")]
# invert the column names
colnames(adjusted) <- c("stop", "start", "primer")

# switch columns to match the alreadyok table
adjusted <- adjusted[, c("start", "stop", "primer")]

# this table contains all starts and stops in proper order
all_good <- rbind.data.frame(alreadyok, adjusted)
colnames(all_good) <- c("start", "stop", "primer")

# get unique starts and stops combinations
tmp_extraction <- vector()
final_extraction <- vector()

# remove duplicate from starts
for (start_pos in unique(all_good$start)) {
	# get positions with same starting point
	subset <- all_good[which(all_good$start==start_pos), ]
	# get the maximum length from stop. start-max_stop will be retained
	max_ex <- max(subset$stop)
	# get primer from subset correspondin to max_ex
	primer <- subset$primer[which(subset$stop==max_ex)]
	# get partial table
	tmp_extraction <- rbind.data.frame(tmp_extraction, cbind.data.frame(start_pos, max_ex, primer))
}
# set colnames
colnames(tmp_extraction) <- c("start", "stop", "primer")

# remove duplicates from stop
for (stop_pos in unique(tmp_extraction$stop)) {
	# get positions with same starting point
	subset <- tmp_extraction[which(tmp_extraction$stop==stop_pos), ]
	# get the maximum length from stop. start-max_stop will be retained
	max_ex <- max(subset$start)
	# get primer from subset correspondin to max_ex
	primer <- subset$primer[which(subset$start==max_ex)]
	# get final table
	final_extraction <- rbind.data.frame(final_extraction, cbind.data.frame(max_ex, stop_pos, primer))
}
# set colnames
colnames(final_extraction) <- c("start", "stop", "primer")

# export
write.table(final_extraction, "trim_here.csv", quote=F, sep=",", row.names=F, col.names=F)

# print out that the script ended
print("004 ends here")
