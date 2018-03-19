#!/bin/bash
# Usage: sh nastybugs.sh sraDir ./hgDir/hg19 ./cardgene ./cardsnp 16 ./outDir 
# Required: id.txt ( List of SRA Accession Numbers - One Acession number per line)
# Required: Download Human Genome Sequence, e.g. ucsc.hg19.fasta, or other host sequence, and save in hgDir 
# cd hgDir
# makeblastdb -in ucsc.hg19.fasta -dbtype nucl -out hg19 #Create BLAST databases for host removal 

if [ "$#" -ne 6 ]; then
	echo "Not Enough Arguments. Please enter sh nastybugs.sh [List of SRA runs] [Host Genome Directory] [CARD Gene Database Directory] [CARD Variant Database Directory] [Cores] [Output Directory]"
	exit 1
fi

sraid= $(readlink -f "$1")
hostGen=$(readlink -f "$2")
cardgene=$(readlink -f "$3")
cardsnp=$(readlink -f "$4")
cores="$5"
outdir=$(readlink -f "$6")

## Getting CARD Gene and SNP databases and create BLAST db
## https://card.mcmaster.ca/download
printf "%s Getting CARD databases\n"
date
df
cd "$3"
echo $PWD
wget https://card.mcmaster.ca/download/0/broadstreet-v2.0.0.tar.gz
tar xvf broadstreet-v2.0.0.tar.gz
makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta -dbtype nucl -out cardgenedb
cp nucleotide_fasta_protein_variant_model.fasta $cardsnp/.
cd $cardsnp
makeblastdb -in nucleotide_fasta_protein_variant_model.fasta -dbtype nucl -out cardsnpdb
cd ..

for sra in $(cat "$1"); do

  printf "processing %s\n" "$sra"
  date
  # First we align to a host so we can subtract host reads 
  magicblast -sra $sra -db $hostGen -num_threads $cores -score 50 -penalty -3 -out $outdir/$sra.human.sam

  # Convert magicblast output sam file to FASTA
  # Getting unmapped reads (-f 4).  For mapped reads, use flag (-F 4)
  if [[ (-s $outdir/$sra.human.sam) ]]; then
      samtools fasta -f 4 $outdir/$sra.human.sam -1 $outdir/${sra}_unmapped_read_one -2 $outdir/${sra}_unmapped_read_two -0 $outdir/${sra}_unmapped_read_zero
  else
      printf "File not found!\n" && exit 0
  fi

  # Determine if SRA is PE or SE. then run magicblast on the nonhost reads. This could be done by looking at the meta data, too FASTA files

  read1_count=$(grep -c "^>" $outdir/${sra}_unmapped_read_one)
  read2_count=$(grep -c "^>" $outdir/${sra}_unmapped_read_two)
  read0_count=$(grep -c "^>" $outdir/${sra}_unmapped_read_zero)

  # Run magicblast using the CARD gene_homology  and CARD variant

  if [[ ( $read2_count > 0 ) && ( $read1_count > 0 ) ]]; then
    printf "Paired End - %d %d %d\n" "$read1_count" "$read2_count" "$read0_count"
    date
    fastx_clipper -i $outdir/${sra}_unmapped_read_one -o $outdir/${sra}_unmapped_read_one_trimmed # is this line still needed? 
    fastx_clipper -i $outdir/${sra}_unmapped_read_two -o $outdir/${sra}_unmapped_read_two_trimmed # is this line still needed? 
    
    magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_gene.sam -db $cardgene/cardgenedb
    magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_snp.sam -db $cardsnp/cardsnpdb

  elif [[ (-s $read1_count) ]]; then
    printf "Single End - %d %d %d\n" "$read1_count" "$read2_count" "$read0_count"
    date
    fastx_clipper -i $outdir/${sra}_unmapped_read_one -o $outdir/${sra}_unmapped_read_one_trimmed
    # Run magicblast using CARD Gene Homology nucleotide_fasta_protein_homolog_model.fasta 
    magicblast -num_threads $cores -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_gene.sam -db $cardgene/cardgenedb
    # Run magicblast using CARD variant  - nucleotide_fasta_protein_variant_model.fasta
    magicblast -num_threads $cores -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_snp.sam -db $cardgene/cardsnpdb 
  else 
    printf "ERROR: no reads to align to CARD databases\n"
  fi

  # Convert the gene_homology aligned SAM to BAM and sort. Delete SAM to save space
  samtools view -bS $outdir/$sra.CARD_gene.sam | samtools sort -o $outdir/$sra.CARD_gene.bam 
  # Convert the gene_variant aligned SAM to BAM and sort. Delete SAM to save space
  samtools view -bS $outdir/$sra.CARD_snp.sam | samtools sort -o $outdir/$sra.CARD_snp.bam # should we output the unaligned reads to a different bam?
  rm -f $outdir/$sra.CARD*.sam

  # Optional - Process the above bamfiles into a output-usable form 
  printf "Done! %s" "$sra"
  date
  printf "\n"
done > log

