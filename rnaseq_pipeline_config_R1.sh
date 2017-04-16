# Config File for all input parameters nescessary to run the pipeline
# Change parameters according to your input

# Genome reference FASTA file (full path)
reference=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Master_Project/data/Human_genome/GRCh38/Gencode_v25/GRCh38.primary_assembly.genome.fa

# Genome annotation File (full path)
annotation=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Master_Project/data/Human_genome/GRCh38/Gencode_v25/gencode.v25.primary_assembly.annotation.gtf
#Fastq input Files (do not need to be trimmed), if single-end reads second file is ignored (full path)
in_fq_1=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Master_Project/data/HEK293/RNA_seq/Raw_sequences/SRR2549078_1.fastq
in_fq_2=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Master_Project/data/HEK293/RNA_seq/Raw_sequences/SRR2549078_2.fastq

#---------------------------------------------------


#Trimmomatic inputs
	# Set number of threads trimmomatic should use (24 seems to work well)
	tr_threads=24

	#Define adapter properties (adapter fasta is distributed with trimmomatic under http://www.usadellab.org/cms/?page=trimmomatic)
	#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold> ;seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed;palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.;simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read. 
	ILLUMINACLIP=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Master_Project/data/Trimmomatic_adapters/TruSeq3-SE.fa:2:30:10
	
	#LEADING=<quality>; quality: Specifies the minimum quality required to keep a base
	LEADING=3 

	#TRAILING=<quality>; quality: Specifies the minimum quality required to keep a base
	TRAILING=3
 
	#SLIDINGWINDOW=<windowSize>:<requiredQuality>; windowSize: specifies the number of bases to average across ; requiredQuality: specifies the average quality required.
	SLIDINGWINDOW=4:15

	#MINLEN=<length> ; length: Specifies the minimum length of reads to be kept.
	MINLEN=36


#-------------------------------------------------


# STAR inputs

	# Set number of cores STAR will use for indexing the genome and read mapping
	st_threads=48

	# speciﬁes path to the directory (henceforth called ”genome directory” where thevgenome indices are stored. This directory has to be created (with mkdir) before STAR run and needs to writing permissions
	genomeDir=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Master_Project/data/Human_genome/GRCh38/Gencode_v25/STAR_genome_index

	# specifies one or more FASTA ﬁles with the genome reference sequences
	genomeFastaFiles=$reference

	#  speciﬁes the path to the ﬁle with annotated transcripts in the standard GTF format
	sjdbGTFfile=$annotation

	# speciﬁes the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads
	sjdbOverhang=99

	# Remove non-canonical splice junctions (Recommended in STAR manual for cufflinks compatibility)
	outFilterIntronMotifs=RemoveNoncanonicalUnannotated


#-----------------------------------------------


# Samtools inputs

	# Sets number of cores Samtools uses for sorting bam files, recommended: 8
	sa_threads=8


#----------------------------------------------


# Cufflinks inputs

	# Sets the number of cores used for transcript assembly
	cl_threads=24

	# Specifies the library preparation. Chosse: fr-unstranded (e.g. unstranded Illumina TruSeq), fr_firststrand (dUTP, NSR, NNSR) , etc.
	librarytype=fr-firststrand




