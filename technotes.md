# VCF format for complex structural variations
## Introduction
VarSim follows the latest VCF format specification [4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf).
## Details
* In VCF Specification 4.3, the only values allowed for `SVTYPE` in `INFO` are `DEL`, `INS`, `DUP`, `INV`, `CNV` and `BND`.
* In `ALT` column, a subtype can be specified, e.g. `<DUP:TANDEM>`.
* In `ALT` column, values allowed for first level symbolic alleles are `DEL`, `INS`, `DUP`, `INV`, `CNV`. Reserved subtypes include `DUP:TANDEM`, `DEL:ME`, `INS:ME`.
* A cut-and-paste translocation is represented in 2 lines, one as `<DEL:TRA>`, the other as `<DUP:TRA>`. That is, they must appear in pairs.
Evaluation will be performed independently on these two lines.
 The two lines can be reported as a single event (translocation) or 2 unrelated events (deletion + duplication).
 If reported as a single event, VarSim uses `TRAID` in `INFO` field to differentiate between different translocations.
* For a duplication (copy-and-paste), it comes with additional INFO fields (`CHR2`, `POS2`, `END2`) to indicate source of duplicated sequences.
If these additional fields are absent, we assume it is tandem duplication otherwise `<DUP>` and `<DUP:ISP>` are essentially the same.
 The first level symbolic allele in `ALT` must be `DUP` to indicate that this line denotes a duplication.
 `SVLEN` will be used to determine length of duplication regardless of `CHR2`, `POS2` and `END2`.
* An additional tag is required in `INFO` field to indicate the direction of inserted duplicate. VarSim recognizes `ISINV`. If it is absent, then we assume it is on positive strand.
* Some programs may report translocations in a different way, their output will be converted to VarSim format before complex variants can be analyzed.
* Subtypes for `DUP` may include `TRA`, `TANDEM` and `ISP` (interspersed). Other subtypes are ignored (treated as if they were not there).
* `INS` must represent novel sequence insertion (i.e., cannot be duplication). VarSim takes this as a truth without enforcing it.
## Examples
### Header
The following fields are honored by VarSim. Other fields may be present, but will be ignored. Note here the type of `CN` is `String` not `Integer`. Its format is similar to `GT`.
```
##fileformat=VCFv4.3
##reference=referenceFileName.fasta
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome of source sequence">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="1-based start position of source sequence">
##INFO=<ID=END2,Number=1,Type=Integer,Description="1-based end position of source sequence">
##INFO=<ID=ISINV,Number=1,Type=Flag,Description="whether a duplication is inverted">
##INFO=<ID=TRAID,Number=1,Type=String,Description="translocation ID">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=String,Description="Copy number for each genotype">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DEL:TRA,Description="Deletion in translocation">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##ALT=<ID=DUP:ISP,Description="Interspersed duplication">
##ALT=<ID=DUP:TRA,Description="Duplication in translocation">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INV,Description="Inversion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT  SAMPLEX SAMPLEY ...
```
### Example complex variants
#### Deletion
```
1	4	.	C	<DEL>	.	.	SVTYPE=DEL;SVLEN=-3	GT	1|0
```
It means reference sequence `1:5-7` is deleted.

#### Tandem duplication
```
1	3	.	T	<DUP:TANDEM>	.	PASS	SVTYPE=DUP;SVLEN=4	GT	1|1
```
The following line is equivalent.
```
1	3	.	T	<DUP>	.	PASS	SVTYPE=DUP;SVLEN=4	GT	1|1
```
It means 4bp of reference sequence (`1:4-7`) is copied and inserted right after `1:3`. Copy number is assumed to be 1 by default.

#### Cut-and-paste translocation
```
1	1	.	A	<DUP:TRA>	.	PASS	SVTYPE=DUP;SVLEN=3;TRAID=1;CHR2=2;POS2=3;END2=5	GT	1/1
2	2	.	T	<DEL:TRA>	.	PASS	SVTYPE=DEL;SVLEN=-3;TRAID=1	GT	1|1
```
It means `2:3-5` is cut and pasted right after `1:1`.

#### Interspersed duplication
```
1	1	.	A	<DUP:ISP>	.	PASS	SVTYPE=DUP;SVLEN=3;CHR2=2;POS2=3;END2=5	GT	1/1
```

The following line is equivalent.
```
1	1	.	A	<DUP>	.	PASS	SVTYPE=DUP;SVLEN=3;CHR2=2;POS2=3;END2=5	GT	1/1
```
It means 3bp of reference sequence (`2:3-5`) is copied and inserted right after `1:1`.

#### Cut-and-paste translocation with inverted duplication
```
1	1	.	A	<DUP:TRA>	.	PASS	SVTYPE=DUP;SVLEN=3;ISINV;TRAID=1;CHR2=2;POS2=3;END2=5	GT	1/1
2	2	.	T	<DEL:TRA>	.	PASS	SVTYPE=DEL;SVLEN=-3;TRAID=1	GT	1|1
```
It means `2:3-5` is cut and pasted right after `1:1` in reverse-complement orientation.
#### Insertion
```
1	3	.	T	<INS>	.	.	SVTYPE=INS;SVLEN=3	GT	1/1
```
It means 3 bp of unknown sequence is inserted after `1:3`.
#### Inversion
```
1	3	.	T	<INV>	.	PASS	SVTYPE=INV;SVLEN=4	GT	0|1
```
It means 4 bp of reference sequence (`1:4-7`) is inverted.

## Implementation details

### `<DUP:TRA>`, `<DUP:ISP>`

As of 12/12/2016，VarSim validates translocation duplication via breakends. Each `<DUP:TRA>` or `<DUP:ISP>` is broken into two nonadjacent novel breakends, i.e. the two breakends at the outermost positions of duplication insertion.

### validation of breakend

Breakends are validated by two criteria:
1. the locus, plus/minus `wiggle`.
2. the oritentation

It's possible that there are inserted sequences before/after breakends, for now, we ignore them.

### `<DEL:TRA>`

As of 12/12/2016, we treat translocation deletion same as regular deletion, i.e. as intervals whose overlapping is subject to control of `overlapRatio` and `wiggle`.

### `<TRAID>`

Translocation ID is used to link `<DUP:TRA>` and `<DEL:TRA>`, each combination acts as if they are one translocation, and will be valided and reported as a single event. Partial matching will have to exceed `overlapRatio` to make this translocation match with another.
