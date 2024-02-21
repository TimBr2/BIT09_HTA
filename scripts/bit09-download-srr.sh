#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: ./bit09-download-srr.sh 
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="This download script requires 3 parameters:\n
		1. File with SRR accession numbers.\n
        2. Path of the output folder.\n
        3. SE for single end or PE for paired end reads.\n
        Run fastq-dump.";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 3 ]; then
	usage
fi
################################################################################
# INPUT FILE CONTAINING SRR ACCESSION NUMBERS
# Note: txt files must end with newline or data after last newline not read
################################################################################
inputFile=$1;
################################################################################
# OUTPUT FOLDER (CREATE IF NOT EXISTS)
################################################################################
outFolder=$2;
# Remove trailing slash if this is last char
len=${#outFolder};
lastPos=$(expr $len - 1);
lastChar=${outFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	outFolder=${outFolder:0:$lastPos};
fi
mkdir -p ${outFolder}
################################################################################
# RUN FASTQ-DUMP
################################################################################
# Read SRR numbers from inputFile and download SE or PE fastq.gz files
while read line
do
    # Remove \r \n from line and echo SRR
    line="$(echo "$line"|tr -d '\r')"
    line="$(echo "$line"|tr -d '\n')"
    echo "### $line ###"
    # Compose command
    fastqdumpCommand="fastq-dump --gzip"; #do -X 100 at the end to only get the first 100 reads, i deleted this to get all the reads
    if [[ $3 == "PE" ]]; then
        # Add --split-3 for paired end reads
        fastqdumpCommand="$fastqdumpCommand --split-3";
    fi
    # Add output folder and SRR
    fastqdumpCommand="$fastqdumpCommand -O ${outFolder} ${line}";
    # Show command
    echo -e "$fastqdumpCommand";
    # Execute command
    output=$(eval $fastqdumpCommand);
    # Show output
    echo -e "$output\n";
done < $inputFile
################################################################################
