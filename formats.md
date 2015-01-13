Formats for different intermediate and final output files generated/used by VarSim
===================

# Map file
Our modified vcf2diploid generates a file indicating the mapping from blocks of the host genome (the simulated genome) to the reference (b37, hg19, etc). The map file is used by the FASTQ-liftover tool to convert the alignment coordinates to the reference. The map file will contain one line per block. For each block, we have the line in the following format
```
<size of the block> <host chr> <host loc> <ref chr> <ref loc> <direction of block> <feature name> <variant_id>
```
Direction of a block can be '+' or '-'. The feature name indicates which feature the block comes from-could be 'INS', 'DEL', 'DUP', 'DUP_TANDEM', 'SEQ', 'INV'. When the feature is 'INS', the 'ref chr' and 'ref loc' fields indicate the location before which the insertion happened. This is vice-versa for 'DEL'. When the direction is '-', then the mapping for the host block is done after reversing the locations within the block-for '+', the mapping is sequential. Also, all the locations are 1-indexed.  

Variant_id is to keep track of variants that have multiple features. 

