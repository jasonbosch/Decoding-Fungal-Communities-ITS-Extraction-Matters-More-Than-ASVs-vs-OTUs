# run example
# line="https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-57/fungi/fasta/fungi_basidiomycota1_collection/ustilago_bromivora_gca_900080155/dna/Ustilago_bromivora_gca_900080155.UBRO_v3.dna.toplevel.fa.gz"

# substitute spaces with underscore
sed 's/ /_/g' download_me_with_spaces.txt | sed 's/,//g' | sed 's/(/_/g' | sed 's/)/_/g' > download_me.txt
rm download_me_with_spaces.txt

# create necessary folders
mkdir all_its_seqs all_genomes all_genomes_bkp

# move to the all_genomes directory
cd all_genomes

# download all genomes using wget
while read download_me; do
	wget $download_me
done < ../download_me.txt

# store all the dowloaded genomes names in a file
ls *gz > ../all_downloaded_genomes.txt

# copy all the downloaded genomes in a backup in case they are needed
cp * ../all_genomes_bkp/

# move back one folder
cd ..

