#!/bin/bash

num_snp=3000000
num_ins=100000
num_del=100000
num_mnp=50000
num_complex=50000
num_sv_ins=2000
num_sv_del=2000
num_inv=1000
num_dup=1000
novel_pc=0.01
nlanes=3
coverage=30
dbsnp=/net/kodiak/volumes/river/shared/users/johnmu/sv_sim/marghoob/dbsnp_138.b37.vcf
reference=/net/kodiak/volumes/lake/shared/references/human/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta
dgv=/net/kodiak/volumes/river/shared/users/johnmu/sv_sim/marghoob/GRCh37_hg19_supportingvariants_2013-07-23.txt
insert_seq=/net/kodiak/volumes/river/shared/users/johnmu/sv_sim/marghoob/insert_seq.txt
art=/home/marghoob/lake/opt/ART/art_illumina
dwgsim=/net/kodiak/volumes/lake/shared/opt/dwgsim/dwgsim

./varsim.py --vc_num_snp $num_snp --vc_num_ins $num_ins --vc_num_del $num_del --vc_num_mnp $num_mnp --vc_num_complex $num_complex --vc_percent_novel $novel_pc --vc_in_vcf $dbsnp \
            --reference $reference --id sv \
            --sv_num_ins $num_sv_ins --sv_num_del $num_sv_del --sv_num_inv $num_inv --sv_num_dup $num_dup --sv_percent_novel $novel_pc --sv_insert_seq $insert_seq --sv_dgv $dgv \
            --dwgsim $dwgsim --nlanes $nlanes --total_coverage $coverage --art $art --out_dir art/out --log_dir art/log --work_dir art/work --use_art \
            --disable_rand_vcf --disable_rand_dgv --disable_vcf2diploid
