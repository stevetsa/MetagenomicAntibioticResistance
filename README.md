# NastyBugs

### A Simple Method for Extracting Antimicrobial Resistance Information from Metagenomes
##### Hackathon team: Lead: Steve Tsang - SysAdmins: Greg Fedewa, Daniel Quang, Sherif Farag - Writers: Matthew Moss, Alexey V. Rakov

#### How to cite this work in a publication:  Tsang H, Moss M, Fedewa G et al. NastyBugs: A simple method for extracting antimicrobial resistance information from metagenomes [version 1; referees: 1 approved with reservations]. F1000Research 2017, 6:1971 
(doi: 10.12688/f1000research.12781.1)

Antibiotic resistance (AMR) of bacterial pathogens is a growing public health threat around the world. Fast and reliable extraction of antimicrobial resistance genomic signatures from large raw sequencing datasets obtained from human metagenomes is a key task for bioinformatics.

**NastyBugs** is a versatile workflow for fast extracting of antimicrobial resistance genomic signatures from metagenomic sequencing data.

*Objective*: Create a reusable, reproducible, scalable, and interoperable workflow 
to locate antimicrobial resistant genomic signatures in SRA shotgun sequencing (including metagenomics) datasets.

This project was part of the [Summer 2017 NCBI Hackathon](https://ncbi-hackathons.github.io/).

## Docker Image

A Docker image complete with all tools and dependencies in project is available.  
[Install Docker](https://docs.docker.com/install/#desktop)  
[How to use/run a Docker image from a previous hackathon](https://github.com/NCBI-Hackathons/Cancer_Epitopes_CSHL/blob/master/doc/Docker.md)  

## Pull and Run Docker Image
From Working Directory - 
Setting environment
```{sh}
touch id.txt #Add SRA Accession - one per line 
mkdir hgDir cardgene cardsnp outDir  #create directories or use existing directories
cp ucsc.hg19.fasta ./hgDir/.         #download host genome sequence and copy to the hgDir directory
cd hgDir
makeblastdb -in ucsc.hg19.fasta -dbtype nucl -out hg19 #Create BLAST databases for host removal
cd ..
```
Run nastybugs 
```{sh}
sh nastybugs.sh id.txt ./hgDir/hg19 ./cardgene ./cardsnp 16 ./outDir
docker run -v `pwd`:`pwd` -w `pwd` -i -t stevetsa/metagenomicantibioticresistance:latest
sh /MetagenomicAntibioticResistance/nastybugs.sh id.txt ./hgDir/hg19 ./cardgene ./cardsnp 16 ./outDir
```

## NastyBugs Workflow

![My image](https://github.com/NCBI-Hackathons/MetagenomicAntibioticResistance/blob/master/AbxResistanceMetagenomics.png)

## Workflow method

The pipeline use three databases that should be downloaded with the script:
1.	**GRCh37/hg19 human reference genome database** used for alignment and filtering reads of human origin from metagenomics samples.
2.	**CARD database** used for search of genomic signatures in the subset of reads unaligned to human genome.
3.	**RefSeq reference bacterial genomes database** used for search and assigning of 16S RNA taxonomic labels the subset of reads unaligned to human genome.

Step 1. 
Getting CARD Gene and SNP databases and create BLAST db - https://card.mcmaster.ca/download

```{sh}
wget https://card.mcmaster.ca/download/0/broadstreet-v2.0.0.tar.gz
tar xvf broadstreet-v2.0.0.tar.gz
makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta -dbtype nucl -out cardgenedb
makeblastdb -in nucleotide_fasta_protein_variant_model.fasta -dbtype nucl -out cardsnpdb
```

Step. 2

Loop Through SRA Accession List (id.txt - one accession per line)
First we align to a host so we can subtract host reads by mapping SRA to host genome using magicblast 
Extract unmapped read and convert to FASTA   #Getting unmapped reads (-f 4).  For mapped reads, use flag (-F 4)  

```{sh}
magicblast -sra $sra -db $hostGen -num_threads $cores -score 50 -penalty -3 -out $outdir/$sra.human.sam
samtools fasta -f 4 $outdir/$sra.human.sam -1 $outdir/${sra}_unmapped_read_one -2 $outdir/${sra}_unmapped_read_two -0 $outdir/${sra}_unmapped_read_zero
```
    
Step. 3
Trim FASTA files
```{sh}
fastx_clipper -i $outdir/${sra}_unmapped_read_one -o $outdir/${sra}_unmapped_read_one_trimmed
fastx_clipper -i $outdir/${sra}_unmapped_read_two -o $outdir/${sra}_unmapped_read_two_trimmed
```

Step. 4     
Map FASTA files to CARD gene and variant sequences

```{sh}
magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_gene.sa
m -db $cardgene/cardgenedb
magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_snp.sam
 -db $cardsnp/cardsnpdb
```

Step. 5
Convert the gene_homology aligned SAM to BAM and sort. Delete SAM to save space

```{sh}
samtools view -bS $outdir/$sra.CARD_gene.sam | samtools sort -o $outdir/$sra.CARD_gene.bam 
samtools view -bS $outdir/$sra.CARD_snp.sam | samtools sort -o $outdir/$sra.CARD_snp.bam # should we output the unaligned reads to a different bam?
rm -f $outdir/$sra.CARD*.sam
```

## Usage
From working directory - 
Required: id.txt ( List of SRA Accession Numbers - One Acession number per line)
Required: Download Human Genome Sequence, e.g. ucsc.hg19.fasta, or other host sequence, and save in ./hgDir
sh nastybugs.sh [List of SRA runs] [Host Genome Directory] [CARD Gene Database Directory] [CARD Variant Database Directory] [Cores] [Output Directory]

```{sh}
mkdir hgDir cardgene cardsnp outDir  #create directories or use existing directories
cp ucsc.hg19.fasta ./hgDir/.
makeblastdb -in ucsc.hg19.fasta -dbtype nucl -out hg19 #Create BLAST databases for host removal
cd ..
sh nastybugs.sh id.txt ./hgDir/hg19 ./cardgene ./cardsnp 16 ./outDir
```

## Input file format

Default - SRA accession numbers (ERR or SRR)
FASTQ files can be used by modifying magicblast steps.

## Validation

The NastyBugs workflow was validated using the next SRAs: ERR1600439 and ERR1600437

## Planned Features
1. Code optimization.
2. Improved more detailed output.
3. Prediction of novel resistance genes (using HMM).

## Dependencies:computer:

*Software:*

[Magic-BLAST 1.3](https://ncbi.github.io/magicblast/) is a tool for mapping large next-generation RNA or DNA sequencing runs against a whole genome or transcriptome.

[SAMtools 1.3.1](http://www.htslib.org/) is a suite of programs for interacting with high-throughput sequencing data.

[FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.

[Docker](https://www.docker.com/) is the leading software container platform.

*DBs used for BLAST databases:*

[NCBI GRCh37/UCSC hg19 human reference genome](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml)

[CARD (Comprehensive Antibiotic Resistance Database)](https://card.mcmaster.ca/)

[RefSeq Reference Bacterial Genomes](https://www.ncbi.nlm.nih.gov/refseq/)

## Optional - Output

1. Table (in CSV or TAB-delimited format) with the next columns:
* RefSeq accession number (Nucleotide)
* Genus
* Resistance gene
* ARO (Antibiotic Resistance Ontology) accession number
* Score (number of mapped reads per 1kb)

2. Dot plot showing relative abundance of antimicrobial resistance/bacterial species in metagenomic sample.

3. Pie chart vizualization of bacterial abundance in the given dataset using Krona ([Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. $

![My image](https://github.com/NCBI-Hackathons/MetagenomicAntibioticResistance/blob/master/MetagenomeVisualization.png)

## F.A.Q.
1. How to cite?

Tsang H, Moss M, Fedewa G et al. NastyBugs: A simple method for extracting antimicrobial resistance information from metagenomes [version 1; referees: awaiting peer review]. F1000Research 2017, 6:1971 [doi: 10.12688/f1000research.12781.1](https://f1000research.com/articles/6-1971/)

2. How to use?

Follow the instructions on this page.

3. What if I need a help?

Feel free to contact authors if you need help.

## Reference

Tsang H, Moss M, Fedewa G, Farag S, Quang D, Rakov AV, Busby B. NastyBugs: A simple method for extracting antimicrobial resistance information from metagenomes [version 1; referees: awaiting peer review]. F1000Research 2017, 6:1971 [doi: 10.12688/f1000research.12781.1](https://f1000research.com/articles/6-1971/)

## People/Team
* [Steve Tsang](https://github.com/stevetsa), NCI/NIH, Gaithersburg, MD, <tsang@mail.nih.gov>
* [Greg Fedewa](https://github.com/harper357), UCSF, San Francisco, CA, <greg.fedewa@gmail.com>
* [Sherif Farag](https://github.com/SWFarag), UNC, Chapel Hill, NC, <farags@email.unc.edu>
* [Matthew Moss](https://github.com/mmoss609), CSHL, Cold Spring Harbor, NY, <moss@cshl.edu>
* [Daniel Quang](https://github.com/daquang), UCI, Irvine, CA, <dxquang@uci.edu>
* [Alexey V. Rakov](https://github.com/alexeyrakov), UPenn, Philadelphia, PA, <rakovalexey@gmail.com>

