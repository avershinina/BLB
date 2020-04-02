#!/usr/bin/bash
# 23 Oct 2018
# A. Vershinina
# Goal: convert BED regions to FASTA files using bedtools. 

HOME=/private/groups/shapirolab/avershinina
FILES=${HOME}/horse_genomes/fastas_w_Ti/*.noXUn.fasta
BED=/private/groups/shapirolab/avershinina/gphocs_genomewide/putative_regions_2.extracted_regions.final.final.bed
OUT=${HOME}/gphocs_genomewide

for SAMPLE in ${FILES}; do
	NAME=$(echo ${SAMPLE} | sed 's/[.].*//' | sed 's:.*/::')
	echo $NAME
	bedtools getfasta -fi $SAMPLE -bed $BED -fo ${OUT}/${NAME}.loci.fa
wait
done
