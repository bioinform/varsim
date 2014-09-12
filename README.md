VarSim: A high-fidelity simulation validation framework for high-throughput genome sequencing with cancer applications (BETA)
===================
# Website
For access to webapps http://bioinform.github.io/varsim/

# System Requirements
<p>
For variant simulation and read generation:
<ul>
<li>32GB free RAM</li>
<li>Enough free disk space to store twice the number of reads generated</li>
<li>ART or dwgsim installed</li>
<li>python 2.x</li>
</ul>
</p>
<p>
For alignment and variant validation:
<ul>
<li>8GB free RAM</li>
</ul>
</p>

# Building
<p>VarSim uses maven to build some of the tools. Please make sure maven is installed on your system.</p>

<p>To build, execute `mvn package`</p>

<p>Contact us if you would like pre-built JARs</p>

# Quick Start
<p>
This quick start guide will provide steps for generating a random genome with pre-specified and random variants. Then generate reads from this genome with ART. Finally, results of analysis on the output of secondary analysis is plotted.  
</p>

<p>
Create a directory `varsim_run`
</p>

Download varsim to `varsim_run`

```
wget http://web.stanford.edu/group/wonglab/varsim/varsim.tar.gz
tar xfz varsim.tar.gz 
```

<p>
Download the following files to `varsim_run`:
<ol>
<li>Reference genome in FASTA format <a href='http://goo.gl/lgT18V'>[hs37d5.fa.gz]</a></li>
<li>Concatented insert sequences as <a href='http://web.stanford.edu/group/wonglab/varsim/insert_seq.txt'>[insert_seq.txt]</a></li>
<li>DGV annotations <a href='http://web.stanford.edu/group/wonglab/varsim/GRCh37_hg19_supportingvariants_2013-07-23.txt'>[GRCh37_hg19_supportingvariants_2013-07-23.txt]</a></li>
<li>dbSNP 141 file <a href='http://goo.gl/NUG0dy'>[All.vcf.gz]</a></li>
<li>(optional) Extra variants to include in simulation variants.vcf</li>
</ol>
</p>

by running the following commands

```
wget http://goo.gl/lgT18V
wget http://web.stanford.edu/group/wonglab/varsim/insert_seq.txt
wget http://web.stanford.edu/group/wonglab/varsim/GRCh37_hg19_supportingvariants_2013-07-23.txt
wget http://goo.gl/NUG0dy
gunzip hs37d5.fa.gz
```

Download Samtools to index the reference genome

```
mkdir samtools
cd samtools
wget http://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
tar xfj samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
cd ..
samtools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools faidx hs37d5.fa
```

Dowload ART to `varsim_run` under an the directory ART by running

```
mkdir ART
cd ART
wget http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz
tar xfz artbinvanillaicecream031114linux64tgz.tgz
cd ..
```

The structure of `varsim_run` should look like

```
\varsim_run
- All.vcf.gz  
- \ART  
- - \art_bin_VanillaIceCream  
- - - ART files ...
- - artbinvanillaicecream031114linux64tgz.tgz
- GRCh37_hg19_supportingvariants_2013-07-23.txt	
- hs37d5.fa  
- hs37d5.fa.fai  
- insert_seq.txt  
- \samtools  
- - samtools files...
- \varsim  
- - varsim files ...
- varsim.tar.gz
```

<p>
Run the following command to generate the simulated genome and reads to 30x depth. Replace the values in square brackets with the appropriate values. This will take a few hours to run. The last --vcfs option is optional and only required if you want to add additional variants to the simluation. 
</p>

```
varsim/varsim.py --vc_in_vcf All.vcf.gz --sv_insert_seq insert_seq.txt \
--sv_dgv GRCh37_hg19_supportingvariants_2013-07-23.txt \
--reference hs37d5.fa --id simu --read_length 100 --vc_num_snp 3000000 --vc_num_ins 100000 \
--vc_num_del 100000 --vc_num_mnp 50000 --vc_num_complex 50000 --sv_num_ins 2000 \
--sv_num_del 2000 --sv_num_dup 200 --sv_num_inv 1000 --sv_percent_novel 0.01 \
--vc_percent_novel 0.01 --mean_fragment_size 350 --sd_fragment_size 50 \
--vc_min_length_lim 0 --vc_max_length_lim 49 --sv_min_length_lim 50 \
--sv_max_length_lim 1000000 --nlanes 3 --total_coverage 30 \
--art ART/art_bin_VanillaIceCream/art_illumina --out_dir out --log_dir log --work_dir work \
--dwgsim ART/art_bin_VanillaIceCream/art_illumina --simulator art  \
--vcfs [Optional VCF file to include, variants.vcf]
```

<p>
The reads will be generated in the out directory.  Run the secondary analysis tools (alignment and variant calling) on those.
The ground truth VCF file is also in the out directory called simu.truth.vcf 
</p>

<p>
After running the alignment and variant calling we can evaluate the results. In order to validate the variants run the following command:
</p>

```
java -jar target/build/vcfcompare.jar simu.truth.vcf [VCF from result of secondary analysis] simu 0.8 20
```

<p>
This will output a JSON file that can be used as input to the VCF Compare webapp [http://bioinform.github.io/varsim/webapp/variant_compare.html]
</p>

<p>
In order to validate the alignments run the following command:
</p>

```
java -jar target/build/samcompare.jar [BAM files from result of secondary analysis] 20 simu
```

<p>
This will output a JSON file that can be used as input to the Alignment Compare webapp [http://bioinform.github.io/varsim/webapp/alignment_compare.html]
</p>

# Testing
<p>example.sh is a sample run of VarSim</p>

# Running
Type `varsim.py -h` for help.

# References/Tools

* vcf2diploid: Rozowsky J, Abyzov A, Wang J, Alves P, Raha D, Harmanci A, Leng J, Bjornson R, Kong Y, Kitabayashi N, Bhardwaj N, Rubin M, Snyder M, Gerstein M. <a href="http://msb.embopress.org/content/7/1/522">AlleleSeq: analysis of allele-specific expression and binding in a network framework</a>. Mol Syst Biol. 2011 Aug 2;7:522. doi: 10.1038/msb.2011.54. <a href="http://alleleseq.gersteinlab.org/tools.html">Download link</a>

* DWGSIM: Nils Homer. <a href="https://github.com/nh13/DWGSIM">Github repository</a>

* ART: Huang W1, Li L, Myers JR, Marth GT. <a href="http://www.ncbi.nlm.nih.gov/pubmed/22199392">ART: a next-generation sequencing read simulator.</a> Bioinformatics. 2012 Feb 15;28(4):593-4. doi: 10.1093/bioinformatics/btr708. Epub 2011 Dec 23. <a href="http://www.niehs.nih.gov/research/resources/software/biostatistics/art/">Download link</a>
