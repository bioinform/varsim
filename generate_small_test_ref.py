import os
import subprocess
import sys
import numpy
import argparse
import pysam
import vcf
import pybedtools
import logging
from intervaltree import Interval, IntervalTree
from collections import defaultdict, OrderedDict

def uint(value):
  if not value.isdigit(): raise argparse.ArgumentTypeError("%s is not digit-only" % value)
  ret = int(value)
  if ret < 0: raise argparse.ArgumentTypeError("%s is negative" % value)
  return ret

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Generate restricted FASTAs and VCFs given a BED file. The contigs are the sequences for each genomic region in the BED file and the name of the contigs reflects that. The VCFs use the coordinates on the new contigs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--reference", help="Reference FASTA", required=True)
parser.add_argument("--regions", help="Regions BED", required=True)
parser.add_argument("--vcfs", nargs="+", required=True, default=[])
parser.add_argument("--outdir", required=True)
parser.add_argument("--flank", type=uint, default=0, help="ignore variants this close to the edges of a region")

args = parser.parse_args()

reference = args.reference
regions = args.regions
invcfs = args.vcfs
outdir = args.outdir

if not os.path.exists(outdir):
  os.makedirs(outdir)


reference_fasta = pysam.Fastafile(reference)
regions_bedtool = pybedtools.BedTool(regions)
regions_intervaltree = defaultdict(IntervalTree)

contigs = []
with open(os.path.join(outdir, "ref.fa"), "w") as out_fasta:
  first = True
  for region in regions_bedtool:
    sequence = reference_fasta.fetch(reference=str(region.chrom), start=region.start, end=region.end)
    region_name = "%s_%d_%d" % (str(region.chrom), region.start, region.end)
    if first:
      out_fasta.write(">{}\n{}".format(region_name, sequence))
      first = False
    else: out_fasta.write("\n>{}\n{}".format(region_name, sequence))
    contigs.append((region_name, len(sequence)))
    regions_intervaltree[str(region.chrom)].add(Interval(region.start, region.end, data=region))
subprocess.check_output("samtools faidx {}/ref.fa".format(outdir), shell=True)
logger.info("Lifter over the reference to %s/ref.fa" % outdir)

# This will liftover the VCFs to the new reference
for invcf in invcfs:
  if not os.path.isfile(invcf):
    logger.error("%s not found" % invcf)
    continue
  # get the base name and use it in the output
  outvcf = os.path.join(outdir, os.path.splitext(os.path.basename(invcf))[0])
  vcf_template_reader = vcf.Reader(open(invcf, "r"))
  vcf_template_reader.metadata["reference"] = os.path.join(outdir, "ref.fa")
  vcf_template_reader.contigs = OrderedDict([(contig_name, vcf.parser._Contig(contig_name, contig_length)) for (contig_name, contig_length) in contigs])
  vcf_writer = vcf.Writer(open(outvcf, "w"), vcf_template_reader)

  tabix_vcf = pysam.TabixFile(invcf, parser=pysam.asVCF())
  info_warned = False
  for region in regions_bedtool:
    records = None
    try: records = tabix_vcf.fetch(reference=str(region.chrom), start=region.start, end=region.end)
    except ValueError: logger.error("Failed to retrieve %s from %s" % (str(region).strip(), invcf))
    if records is None: continue
    for record in records:
      if record.pos <= region.start + args.flank or record.pos + len(record.ref) + args.flank - 1 >= region.end: continue
      record.contig = "%s_%d_%d" % (str(region.chrom), region.start, region.end)
      # record.pos seems to be zero-based, at least in the infinite wisdom of my version of pysam
      record.pos = record.pos - region.start
      fmt = "GT" if len(record) else None
      sample_indexes = [0] if len(record) else []
      info = None
      try: info = vcf_template_reader._parse_info(record.info)
      except ValueError:
        if not info_warned:
          logger.error("Failed to process INFO %s" % str(record.info))
          logger.error("more warnings will be suppressed")
          info_warned = True
      vcfrecord = vcf.model._Record(record.contig, record.pos + 1, record.id, record.ref, [record.alt,], record.qual, record.filter, info, fmt, sample_indexes)
      if len(record):
        vcfrecord.samples = vcf_template_reader._parse_samples([record[0]], "GT", vcfrecord)
      vcf_writer.write_record(vcfrecord)
  vcf_writer.close()
  pysam.tabix_index(outvcf, force=True, preset='vcf')
  logger.info("Lifter over the VCF %s to %s" % (invcf, outvcf))
