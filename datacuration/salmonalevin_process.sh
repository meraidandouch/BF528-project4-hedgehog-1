#!/bin/bash
r1_1=singlecell_pancreas/data/SRR3879604_1.fastq
r1_2=singlecell_pancreas/data/SRR3879605_1.fastq
r1_3=singlecell_pancreas/data/SRR3879606_1.fastq
r2_1=singlecell_pancreas/data/SRR3879604_2.fastq
r2_2=singlecell_pancreas/data/SRR3879605_2.fastq
r2_3=singlecell_pancreas/data/SRR3879606_2.fastq

cat WL_just_BC.txt > WL_just_BC_file.txt
cat tg_map.tsv > tg_map_file.tsv

salmon alevin -l ISR \
	-1 $r1_1 $r1_2 $r1_3 \
	-2 $r2_1 $r2_2 $r2_3 \
	--index salmon_index \
	--whitelist WL_just_BC_file.txt \
	-o GSM2230760__salmon_quant \
	--tgMap tg_map_file.tsv \
    --end 5 --barcodeLength 19 --umiLength 6 \
	--dumpMtx
echo "Job finished"
