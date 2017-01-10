#!/bin/bash

set -x

b37_source="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
dbsnp_source="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz"

# Download varsim
version=$(java -jar VarSim.jar -version)
wget https://github.com/bioinform/varsim/releases/download/v$version/varsim-$version.tar.gz
tar xfz varsim-$version.tar.gz

# Download reference and variant databases 
wget $b37_source -O - | gunzip -c > hs37d5.fa
wget http://web.stanford.edu/group/wonglab/varsim/insert_seq.txt
wget http://web.stanford.edu/group/wonglab/varsim/GRCh37_hg19_supportingvariants_2013-07-23.txt
wget $dbsnp_source -O All.vcf.gz

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
