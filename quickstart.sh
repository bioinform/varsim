#!/bin/bash

set -x

# Download varsim
wget https://github.com/bioinform/varsim/releases/download/v0.4.1/varsim-0.4.1.tar.gz
tar xfz varsim-0.5.tar.gz

# Download reference and variant databases 
wget http://goo.gl/lgT18V
wget http://web.stanford.edu/group/wonglab/varsim/insert_seq.txt
wget http://web.stanford.edu/group/wonglab/varsim/GRCh37_hg19_supportingvariants_2013-07-23.txt
wget http://goo.gl/NUG0dy
gunzip hs37d5.fa.gz

# Download samtools and index reference
mkdir samtools
cd samtools
wget http://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
tar xfj samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
cd ..
samtools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools faidx hs37d5.fa

# Download ART
mkdir ART
cd ART
wget http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz
tar xfz artbinvanillaicecream031114linux64tgz.tgz
cd ..

# Test run varsim to generate 1x coverage data
./varsim.py --vc_in_vcf All.vcf.gz --sv_insert_seq insert_seq.txt \
--sv_dgv GRCh37_hg19_supportingvariants_2013-07-23.txt \
--reference hs37d5.fa --id simu --read_length 100 --vc_num_snp 3000000 --vc_num_ins 100000 \
--vc_num_del 100000 --vc_num_mnp 50000 --vc_num_complex 50000 --sv_num_ins 2000 \
--sv_num_del 2000 --sv_num_dup 200 --sv_num_inv 1000 --sv_percent_novel 0.01 \
--vc_percent_novel 0.01 --mean_fragment_size 350 --sd_fragment_size 50 \
--vc_min_length_lim 0 --vc_max_length_lim 49 --sv_min_length_lim 50 \
--sv_max_length_lim 1000000 --nlanes 3 --total_coverage 1 \
--simulator_executable ART/art_bin_VanillaIceCream/art_illumina --out_dir out --log_dir log --work_dir work \
--simulator art
