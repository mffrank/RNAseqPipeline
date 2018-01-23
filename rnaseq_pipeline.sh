#!/bin/bash

# Input arguments to script
PAIRED=$1
STAGE=$2
RUN_NAME=$3
OUT_DIR=$4
PARAMS=$5
in_fq_1=$6
in_fq_2=$7
# show usage if we don't have the last command-line argument

if [ "$in_fq_1" == "" ]; then
	echo "-----------------------------------------------------------------";
	echo "Trims, aligns and assembles transcript reads from illumina sequencing runs to a transcriptome fasta on the Euler cluster. Parameters of the individual modules can be set in the rnaseq_pipeline_config.sh file. Outputs can be found in the specified output folder and consist of .sam and .bam alignments from STAR, transcripts in .gtf and .fa format from cufflinks. ";
	echo "";
	echo "USAGE INFORMATION:";
	echo "";
	echo "rna-seq-pipeline PAIRED[0|1] STAGE[0:4]  RUN_NAME OUT_DIR Path/to/parameter_file Path/to/input/fasta";
	echo "";
	echo "PAIRED	1=Paired end sequencing; 0=Single end sequencing";
	echo "STAGE	Start with 0:FASTQC, 1: TRIMMOMATIC, 2: STAR, 3: SAMTOOLS, 4:CUFFLINKS";
	echo "RUN_NAME	Name of the Sequencing run, output files will be saved under that name"
	echo "OUT_DIR	Path to desired output directory";
	echo "PARAMS Path/to/parameter_file Path to shell script containing parameter variables";
	echo "in_fq_1/2 Fastq input Files (do not need to be trimmed), if single-end reads second file is ignored (full path)"
	echo "";
	exit 1;
fi
#Function that exits after printing its text argument
#   in a standard format which can be easily grep'd.
err() {
  echo "$1...exiting";
  exit 1; # any non-0 exit code signals an error
}
# function to check return code of programs.
# exits with standard message if code is non-zero;
# otherwise displays completiong message and date.
#   arg 1 is the return code (usually $?)
#   arg2 is text describing what ran
ckRes() {
  if [ "$1" == "0" ]; then
    echo "..Done $2 `date`";
  else
    err "$2 returned non-0 exit code $1";
  fi
}
# function that checks if a file exists
#   arg 1 is the file name
#   arg2 is text describing the file (optional)
ckFile() {
  if [ ! -e "$1" ]; then
    err "$2 File '$1' not found";
  fi
}
# function that checks if a file exists and
#   that it has non-0 length. needed because
#   programs don't always return non-0 return
#   codes, and worse, they also create their
#   output file with 0 length so that just
#   checking for its existence is not enough
#   to ensure the program ran properly
ckFileSz() {
  ckFile $1 $2;
  SZ=`ls -l $1 | awk '{print $5}'`;
  if [ "$SZ" == "0" ]; then
    err "$2 file '$1' is zero length";
  else
    echo "$2 file '$1' checked";
  fi
}

# Load Parameters from config File
ckFile $PARAMS "Config File";
source $PARAMS

#Check if Paired is set to 0 or 1

if [ $PAIRED -ne 1 -a $PAIRED -ne 0 ]; then err "PAIRED must be 0 or 1"; fi

# Check if input files are present

echo "";
ckFileSz $in_fq_1 "1st FASTQ Input";
if [ $PAIRED == 1 ]; then  ckFileSz $in_fq_2 "2nd FASTQ Input"; fi
ckFileSz $reference "Reference Fasta";
ckFileSz $annotation "Genome annotation";
echo "";

# Set up folder structure for output files

if [ ! -d ${OUT_DIR} ]; then mkdir ${OUT_DIR}; fi #OUT_DIR s a parameter defined in the function input
if [ ! -d ${OUT_DIR}/fastqc ]; then mkdir ${OUT_DIR}/fastqc; fi
if [ ! -d ${OUT_DIR}/trimmomatic ]; then mkdir ${OUT_DIR}/trimmomatic; fi
if [ ! -d ${OUT_DIR}/star ]; then mkdir ${OUT_DIR}/star; fi
if [ ! -d ${OUT_DIR}/samtools ]; then mkdir ${OUT_DIR}/samtools; fi
if [ ! -d ${OUT_DIR}/cufflinks ]; then mkdir ${OUT_DIR}/cufflinks; fi
if [ ! -d ${OUT_DIR}/cufflinks ]; then err "OUT_DIR variable seems not to have been set correctly"; fi

#Loading requred modules

module load java
module load boost/1.55.0
module load gdc
module load fastqc
module load trimmomatic
module load star
module load samtools
module load cufflinks

# Generate an identifier to mark jobs
id=$RANDOM

# Start a Test job if starting at later stages

case $STAGE in
	1) echo "Starting at Stage 1 - Adapter trimming with Trimmomatic";;
	2) bsub -J job2 echo "Starting at Stage 2 - Genomic mapping with Star";echo "Starting at Stage 2";;
	3) bsub -J job4 echo "Starting at Stage 3 - SAM to BAM conversion and BAM sorting with Samtools";echo "Starting at Stage 3";;
	4) bsub -J job6 echo "Starting at Stage 4";echo "Starting at Stage 4 - Transcriptome assembly with cufflinks";;
esac


#--------------
# ACTUAL WORK
#--------------

if [ $PAIRED -eq 1 ]; then

	## First QC Step with fastqc
		if [ $STAGE -eq 0 ]; then
		echo "Submitting FastQC-raw read quality control job"
		bsub -J job1 -o $OUT_DIR/fastqc/${RUN_NAME}_untrimmed_fastqc_euler_log.txt -e $OUT_DIR/fastqc/${RUN_NAME}_untrimmed_fastqc_log.txt fastqc $in_fq_1 $in_fq_2 -o $OUT_DIR/fastqc

	# needs only one core, saves report as zip file in output directory

	fi
	if [ $STAGE -le 1 ]; then

	#Trimmomatic base command
		echo "Submitting Trimmomatic job"
		bsub -n $tr_threads -J job2 -o $OUT_DIR/trimmomatic/${RUN_NAME}_trimmomatic_euler_log.txt -e $OUT_DIR/trimmomatic/${RUN_NAME}_trimmomatic_log.txt trimmomatic PE -threads $tr_threads -trimlog $OUT_DIR/trimmomatic/${RUN_NAME}_trimming_log.txt $in_fq_1 $in_fq_2 -baseout $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed.fastq ILLUMINACLIP:$ILLUMINACLIP LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$SLIDINGWINDOW MINLEN:$MINLEN
		#Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided).  Initially
		#Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches.
		#These seeds will be extended and clipped if in the case of paired end reads a score of
		#30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases).

		#Remove leading low quality or N bases (below quality 3)

		#Remove trailing low quality or N bases (below quality 3)

		#Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15

		#Drop reads which are less than 36bases long after these steps

	fi
	if [ $STAGE -le 2 ]; then

		# Second FASTQC run after trimming
		echo "Submitting FastQC-trimmed read quality control job"
		bsub -J job2b -w "done(job2)" -o $OUT_DIR/fastqc/${RUN_NAME}_trimmed_fastqc_euler_log.txt -e $OUT_DIR/fastqc/${RUN_NAME}_trimmed_fastqc_log.txt fastqc $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed_1P.fastq $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed_2P.fastq -o $OUT_DIR/fastqc

		# needs only one core, saves report as zip file in output directory

		#STAR Genome Index generation

		if [ ! -s ${genomeDir}/Genome ]; then

			GI=0
			echo "No indexed Genome found in '$genomeDir'. Creating one from scratch!"
			bsub -n $st_threads -J job3 -o $genomeDir/STAR_gi__euler_log.txt -e $genomeDir/STAR_gi_log.txt STAR --runThreadN $st_threads --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --sjdbOverhang $sjdbOverhang
		else
			GI=1
			echo "Found indexed genome in '$genomeDir'. Starting to map reads with Star..."
		fi

		#STAR Mapping
		if [ $GI -eq 1 ]; then
			# Genome index exists - Wait only for trimming to finish
			bsub -n $st_threads -J job4 -w "done(job2)" -o $OUT_DIR/star/${RUN_NAME}_STAR_euler_log.txt -e $OUT_DIR/star/${RUN_NAME}_STAR_log.txt STAR --runThreadN $st_threads --genomeDir $genomeDir --readFilesIn $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed_1P.fastq $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed_2P.fastq --outFileNamePrefix $OUT_DIR/star/$RUN_NAME --outFilterIntronMotifs $outFilterIntronMotifs
		else
			# Genome index doesnt exist - Wait for trimming AND Genome Indexiing to finish
			bsub -n $st_threads -J job4 -w "done(job2) && done(job3)" -o $OUT_DIR/star/${RUN_NAME}_STAR_euler_log.txt -e $OUT_DIR/star/${RUN_NAME}_STAR_log.txt STAR --runThreadN $st_threads --genomeDir $genomeDir --readFilesIn $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed_1P.fastq $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed_2P.fastq --outFileNamePrefix $OUT_DIR/star/$RUN_NAME --outFilterIntronMotifs $outFilterIntronMotifs

		fi

	fi
	if [ $STAGE -le 3 ]; then

		# SAMTOOLS Conversion
		echo "Submitting Samtools SAM/BAM conversion job"
		bsub  -J job5 -w "done(job4)" -o $OUT_DIR/samtools/${RUN_NAME}_STconversion_euler_log.txt -e $OUT_DIR/samtools/${RUN_NAME}_STconversion_log.txt samtools view -b -u -o  $OUT_DIR/samtools/${RUN_NAME}Aligned.out.bam $OUT_DIR/star/${RUN_NAME}Aligned.out.sam

		echo "Submitting Samtools BAM sort job"

		bsub -n $sa_threads -J job6 -w "done(job5)" -o $OUT_DIR/samtools/${RUN_NAME}_STsort_euler_log.txt -e $OUT_DIR/samtools/${RUN_NAME}_STsort_log.txt samtools sort -@ $sa_threads -m 1G -o $OUT_DIR/samtools/${RUN_NAME}Aligned.out_sorted.bam $OUT_DIR/samtools/${RUN_NAME}Aligned.out.bam

	fi
	if [ $STAGE -le 4 ]; then

		# Cufflinks assembly
		echo "Submitting Cufflinks Aglinment job"
		bsub -n $cl_threads -J job7 -w "done(job6)" -o $OUT_DIR/cufflinks/${RUN_NAME}_cufflinks_assembly_euler_log.txt -e $OUT_DIR/cufflinks/${RUN_NAME}_cufflinks_assembly_log.txt cufflinks -p $cl_threads -G $annotation -b $reference --library-type $librarytype $OUT_DIR/samtools/${RUN_NAME}Aligned.out_sorted.bam --output-dir $OUT_DIR/cufflinks

		# Conversion to transcriptome FASTA
		echo "Submitting Transcriptome-GTF to Fasta conversion job"
		bsub  -J job8 -w "done(job7)" -N -o $OUT_DIR/cufflinks/${RUN_NAME}_gffread_fasta_conversion_euler_log.txt -e $OUT_DIR/cufflinks/${RUN_NAME}_gffread_fasta_conversion_log.txt gffread $OUT_DIR/cufflinks/transcripts.gtf -g $reference -w $OUT_DIR/cufflinks/transcriptome.fa

		echo "All jobs submitted successfully. You can monitor them with bjobs. An email will be sent to your NETHZ adress when the last job has finished."


	else
		 echo "Please Specify a STAGE value between 1 and 4: 0:FASTQC, 1: TRIMMOMATIC, 2: STAR, 3: SAMTOOLS, 4:CUFFLINKS"
	fi


else
	echo "Running Single-end mode..."
	#err "Single end read Processing is not yet supported";
	## First QC Step with fastqc
	if [ $STAGE -eq 0 ]; then
		echo "Submitting FastQC-raw read quality control job"
		bsub -J job1 -o $OUT_DIR/fastqc/${RUN_NAME}_untrimmed_fastqc_euler_log.txt -e $OUT_DIR/fastqc/${RUN_NAME}_untrimmed_fastqc_log.txt fastqc $in_fq_1 -o $OUT_DIR/fastqc

	# needs only one core, saves report as zip file in output directory
	fi

	if [ $STAGE -le 1 ]; then

	#Trimmomatic base command
		echo "Submitting Trimmomatic job"
		bsub -n $tr_threads -J job2 -o $OUT_DIR/trimmomatic/${RUN_NAME}_trimmomatic_euler_log.txt -e $OUT_DIR/trimmomatic/${RUN_NAME}_trimmomatic_log.txt trimmomatic SE -threads $tr_threads -trimlog $OUT_DIR/trimmomatic/${RUN_NAME}_trimming_log.txt $in_fq_1 $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed.fastq ILLUMINACLIP:$ILLUMINACLIP LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$SLIDINGWINDOW MINLEN:$MINLEN
		#Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided).  Initially
		#Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches.
		#These seeds will be extended and clipped if in the case of paired end reads a score of
		#30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases).

		#Remove leading low quality or N bases (below quality 3)

		#Remove trailing low quality or N bases (below quality 3)

		#Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15

		#Drop reads which are less than 36bases long after these steps
	fi

	if [ $STAGE -le 2 ]; then

		# Second FASTQC run after trimming
		echo "Submitting FastQC-trimmed read quality control job"
		bsub -J job2b -w "done(job2)" -o $OUT_DIR/fastqc/${RUN_NAME}_trimmed_fastqc_euler_log.txt -e $OUT_DIR/fastqc/${RUN_NAME}_trimmed_fastqc_log.txt fastqc $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed.fastq -o $OUT_DIR/fastqc

		# needs only one core, saves report as zip file in output directory

		#STAR Genome Index generation

		if [ ! -s ${genomeDir}/Genome ]; then

			GI=0
			echo "No indexed Genome found in '$genomeDir'. Creating one from scratch!"
			bsub -n $st_threads -J job3 -o $genomeDir/STAR_gi_euler_log.txt -e $genomeDir/STAR_gi_log.txt STAR --runThreadN $st_threads --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --sjdbOverhang $sjdbOverhang
		else
			GI=1
			echo "Found indexed genome in '$genomeDir'. Starting to map reads with Star..."
		fi

		#STAR Mapping
		if [ $GI -eq 1 ]; then
			# Genome index exists - Wait only for trimming to finish
			bsub -n $st_threads -R "rusage[scratch=30000]" -J job4 -w "done(job2)" -o $OUT_DIR/star/${RUN_NAME}_STAR_euler_log.txt -e $OUT_DIR/star/${RUN_NAME}_STAR_log.txt STAR --runThreadN $st_threads --genomeDir $genomeDir --readFilesIn $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed.fastq  --outFileNamePrefix $OUT_DIR/star/$RUN_NAME --outFilterIntronMotifs $outFilterIntronMotifs  --outSAMstrandField intronMotif
		else
			# Genome index doesnt exist - Wait for trimming AND Genome Indexiing to finish
			bsub -n $st_threads -J job4 -w "done(job2) && done(job3)" -o $OUT_DIR/star/${RUN_NAME}_STAR_euler_log.txt -e $OUT_DIR/star/${RUN_NAME}_STAR_log.txt STAR --runThreadN $st_threads --genomeDir $genomeDir --readFilesIn $OUT_DIR/trimmomatic/${RUN_NAME}_trimmed.fastq  --outFileNamePrefix $OUT_DIR/star/$RUN_NAME --outFilterIntronMotifs $outFilterIntronMotifs --outSAMstrandField intronMotif

		fi
	fi

	if [ $STAGE -le 3 ]; then

		# SAMTOOLS Conversion
		echo "Submitting Samtools SAM/BAM conversion job"
		bsub  -J job5 -w "done(job4)" -o $OUT_DIR/samtools/${RUN_NAME}_STconversion_euler_log.txt -e $OUT_DIR/samtools/${RUN_NAME}_STconversion_log.txt samtools view -b -u -o  $OUT_DIR/samtools/${RUN_NAME}Aligned.out.bam $OUT_DIR/star/${RUN_NAME}Aligned.out.sam

		echo "Submitting Samtools BAM sort job"

		bsub -n $sa_threads -J job6 -w "done(job5)" -o $OUT_DIR/samtools/${RUN_NAME}_STsort_euler_log.txt -e $OUT_DIR/samtools/${RUN_NAME}_STsort_log.txt samtools sort -@ $sa_threads -m 1G -o $OUT_DIR/samtools/${RUN_NAME}Aligned.out_sorted.bam $OUT_DIR/samtools/${RUN_NAME}Aligned.out.bam

        fi

	if [ $STAGE -le 4 ]; then

		# Cufflinks assembly
		echo "Submitting Cufflinks Aglinment job"
		bsub -n $cl_threads -J job7 -w "done(job6)" -o $OUT_DIR/cufflinks/${RUN_NAME}_cufflinks_assembly_euler_log.txt -e $OUT_DIR/cufflinks/${RUN_NAME}_cufflinks_assembly_log.txt cufflinks -p $cl_threads -G $annotation -b $reference --library-type $librarytype $OUT_DIR/samtools/${RUN_NAME}Aligned.out_sorted.bam --output-dir $OUT_DIR/cufflinks

		# Conversion to transcriptome FASTA
		echo "Submitting Transcriptome-GTF to Fasta conversion job"
		bsub  -J job8 -w "done(job7)" -N -o $OUT_DIR/cufflinks/${RUN_NAME}_gffread_fasta_conversion_euler_log.txt -e $OUT_DIR/cufflinks/${RUN_NAME}_gffread_fasta_conversion_log.txt gffread $OUT_DIR/cufflinks/transcripts.gtf -g $reference -w $OUT_DIR/cufflinks/transcriptome.fa

		echo "All jobs submitted successfully. You can monitor them with bjobs. An email will be sent to your NETHZ adress when the last job has finished."



	else
		 echo "Please Specify a STAGE value between 1 and 4: 0:FASTQC, 1: TRIMMOMATIC, 2: STAR, 3: SAMTOOLS, 4:CUFFLINKS"
	fi

fi
