import os
import sys
import argparse
import vcf
import pysam
import logging
import copy
from collections import defaultdict, OrderedDict

def lift_vcfs(vcfs, out_vcf, reference):
  logger = logging.getLogger(lift_vcfs.__name__)

  if not vcfs:
    logger.info("No VCFs to lift")
    return

  contigs = []

  if reference:
    fasta_handle = pysam.FastaFile(reference)
    contigs = [(contig, fasta_handle.get_reference_length(contig)) for contig in fasta_handle.references]
    fasta_handle.close()

  vcf_template_reader = vcf.Reader(open(vcfs[0], "r"))

  if reference:
    vcf_template_reader.metadata["reference"] = reference
    vcf_template_reader.contigs = OrderedDict([(contig_name, vcf.parser._Contig(contig_name, contig_length)) for (contig_name, contig_length) in contigs])
  vcf_writer = vcf.Writer(open(out_vcf, "w"), vcf_template_reader)
  for index, vcf_ in enumerate(vcfs):
    with open(vcf_) as vcf_fd:
      vcf_reader = vcf.Reader(vcf_fd)
      for record in vcf_reader:
        original_fields = record.CHROM.rsplit("_", 2)
        original_chrom = original_fields[0]
        original_start = int(original_fields[1])
        info = None
        try:
          #info = vcf_template_reader._parse_info(record.INFO)
          info = copy.deepcopy(record.INFO)
          if "END" in info:
            info["END"] += original_start
          if "CHR2" in info:
            chr2_fields = info["CHR2"].rsplit("_", 2)
            chr2_start = int(chr2_fields[1])
            info["CHR2"] = chr2_fields[0]
            info["POS2"] += chr2_start
            info["END2"] += chr2_start
        except ValueError:
          logger.error("Failed to process INFO %s" % str(record.INFO))

        new_record = vcf.model._Record(original_chrom, record.POS + original_start, record.ID, record.REF, record.ALT, record.QUAL, record.FILTER, info, record.FORMAT, record._sample_indexes, record.samples)
        vcf_writer.write_record(new_record)
  vcf_writer.close()

  logger.info("Finished liftover of VCF to original reference")


def lift_maps(maps, out_map):
  logger = logging.getLogger(lift_maps.__name__)

  if not maps:
    logger.info("No MAP files to lift")
    return

  with open(out_map, "w") as out_fd:
    for map_f in maps:
      logger.info("Lifting over " + map_f)
      with open(map_f) as map_fd:
        for line in map_fd:
          if not line.strip():
            continue
          fields = line.strip().split("\t")
          ref_fields = fields[3].split("_")
          offset = int(ref_fields[1])
          out_fd.write("\t".join(fields[:3] + [ref_fields[0], str(int(fields[4]) + offset)] + fields[5:]) + "\n")


def main():
  FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
  logging.basicConfig(level=logging.INFO, format=FORMAT)
  logger = logging.getLogger(main.__name__)

  parser = argparse.ArgumentParser(description="Lift over truth VCF and MAP files to a proper reference in case of a reference derived by restricting to BED regions", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("--vcfs", nargs="+", help="VCFs to liftover", default=[])
  parser.add_argument("--maps", nargs="+", help="MAPs to liftover", default=[])
  parser.add_argument("--reference", help="Original reference", default=None)
  parser.add_argument("--out_dir", help="Output directory", default=".")

  args = parser.parse_args()

  lift_vcfs(args.vcfs, os.path.join(args.out_dir, "truth.vcf"), args.reference)
  lift_maps(args.maps, os.path.join(args.out_dir, "truth.map"))

if __name__ == "__main__":
  main()
