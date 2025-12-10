# after downloading all genomes

# move to all_its_seqs folder and start looking for ITS regions
cd all_its_seqs

# example a_genome
# a_genome="Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa.gz"

# for each genome
while read a_genome; do

	# print out
	echo $a_genome
	# get its file name
	filegz=`echo $a_genome | cut -d"/" -f12`
	# unzip genome
	gzip -d -k "../all_genomes/"$filegz
	# remove gz from filename
	file=`echo $filegz | sed 's/.gz//g'`
	# parse file
	cat "../all_genomes/"$file | grep -v ">"  | tr -d "\n" > tmp_output.txt
	# find primers
	primersearch -seqall tmp_output.txt -infile ../primers.txt -mismatchpercent 20 -outfile found_primers_tmp
	# filter sequences by length as suggested by Jason
	grep -E -B 5 "(length: [0-9]{3} bp)|(length: 1[0-4][0-9]{2} bp)" found_primers_tmp > found_primers.txt
	# parse file, again
	cat tmp_output.txt | tr -d "\n" | tr ">" "\n" | sed "s/^//g" | sed "s/REF/REF\t/g" > output.txt
	# find trimming positions as suggested by Jason
	# this new version for detecting forward or reverse primer in forward reads
	paste -d"," <(grep forward found_primers.txt) <(grep length found_primers.txt) > trimming_positions.txt

	# test whether trimming_positions.txt were found
	# if yes, then go ahead with the analysis
	if [ -s trimming_positions.txt ]
		then
			echo "running 004 and 005"

			# run rscript to find sequences
			Rscript ../004_find_trimming_locations.R

			#get organism name
			orgn=`echo $file | cut -d"." -f1-4`
			# at this point trim sequences
			while read a_pos; do

				# get start and stop positions
				start_pos=`echo $a_pos | cut -d"," -f1`
				stop_pos=`echo $a_pos | cut -d"," -f2`

				# add header to file, so it becomes fasta-like
				echo ">"$orgn"|"$start_pos"|"$stop_pos >> $file"_reads.fasta"

				# get the substring
				cat output.txt | cut -c $start_pos-$stop_pos >> $file"_reads.fasta"
			done < trim_here.csv

			# run script to change the primers to have exact match
			Rscript ../005_substitute_primers.R

			# remove the old fasta, to keep only the substituted-primers version
			rm $file"_reads.fasta"
	fi
	# remove the original fasta, too. keep the .gz only
	rm "../all_genomes/"$file
# file to read from
done < ../all_downloaded_genomes.txt

# delete list with all genomes
rm tmp_output.txt found_primers.txt output.txt trimming_positions.txt trim_here.csv found_primers_tmp

echo "all dowloaded and extracted!"

