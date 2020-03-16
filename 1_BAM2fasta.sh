#!/usr/bin/bash
# 23 Oct 2018
# A. Vershinina
# Goal: convert BAM files to FASTA files using ANGSD. 
# If your genomes are low coverage, you may want to exclude transitions. The foloowing script does not do that.

HOME=/avershinina
FILES=${HOME}/horse_genomes/*.bam
ANGSD=${HOME}/tools/angsd/angsd
FILTERS='-minQ 25 -minMapQ 25 -uniqueOnly -setMinDepth 5 -setMaxDepth 100 -iupacRatio 0.35 -nThreads 10 -howOften 1000'
REF=${HOME}/horse_genomes/ecab2.fasta

echo "Starting conversion, transitions included"

for SAMPLE in ${FILES}; do
	NAME=$(echo ${SAMPLE} | sed 's/[.].*//')
	echo "Making w_Ti fasta for ${NAME}"
	$ANGSD -doFasta 4 -doCounts 1 $FILTERS -i $SAMPLE -out ${NAME}_w_Ti
done
