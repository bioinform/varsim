import os
import subprocess
import sys
import numpy
import argparse
import pysam
import vcf
import pybedtools
import logging
from collections import defaultdict, OrderedDict

def uint(value):
  if not value.isdigit(): raise argparse.ArgumentTypeError("%s is not digit-only" % value)
  ret = int(value)
  if ret < 0: raise argparse.ArgumentTypeError("%s is negative" % value)
  return ret

def gen_restricted_reference(reference, regions_bed, out_reference, use_short_contigs_names=False):
    logger = logging.getLogger(gen_restricted_reference.__name__)

    reference_handle = pysam.Fastafile(reference)
    regions_bedtool = pybedtools.BedTool(regions_bed)

    with open(out_reference, "w") as out_fasta:
        for region_index, region in enumerate(regions_bedtool, start=1):
            sequence = reference_handle.fetch(reference=str(region.chrom), start=region.start, end=region.end)
            region_name = str(region_index) if use_short_contigs_names else ("%s_%d_%d" % (str(region.chrom), region.start, region.end) )
            if region_index == 1:
                out_fasta.write(">{}\n{}".format(region_name, sequence))
            else: out_fasta.write("\n>{}\n{}".format(region_name, sequence))
    pysam.faidx(out_reference)
    logger.info("Lifted over the reference to {}".format(out_reference))

    reference_handle.close()

    return out_reference


def gen_restricted_vcf(in_vcf, regions_bed, out_vcf, restricted_reference, targeted_samples, flank=0, use_short_contig_names=False):
    logger = logging.getLogger(gen_restricted_vcf.__name__)

    if not in_vcf:
        return None

    if not os.path.isfile(in_vcf):
        logger.error("%s not found" % in_vcf)
        return None

    reference_handle = pysam.Fastafile(restricted_reference)
    contigs = list(zip(reference_handle.references, reference_handle.lengths))
    reference_handle.close()

    # get the base name and use it in the output
    vcf_template_reader = vcf.Reader(open(in_vcf, "r"))
    vcf_template_reader.metadata["reference"] = restricted_reference
    vcf_template_reader.contigs = OrderedDict([(contig_name, vcf.parser._Contig(contig_name, contig_length)) for (contig_name, contig_length) in contigs])

    new_samples = []
    if targeted_samples:
        for k,v in sorted(vcf_template_reader._sample_indexes.iteritems()):
            if k in targeted_samples:
                new_samples.append(k)
        vcf_template_reader.samples = new_samples

    vcf_writer = vcf.Writer(open(out_vcf, "w"), vcf_template_reader)

    if targeted_samples:
        vcf_template_reader = vcf.Reader(open(in_vcf, "r"))

    #tabix_vcf = pysam.TabixFile(invcf, parser=pysam.asVCF())
    regions_bedtool = pybedtools.BedTool(regions_bed)

    for region_index, region in enumerate(regions_bedtool, start=1):
        records = None
        try: records = vcf_template_reader.fetch(chrom=str(region.chrom), start=region.start, end=region.end)
        except ValueError: logger.error("Failed to retrieve %s from %s" % (str(region).strip(), in_vcf))
        if records is None: continue
        for record in records:
            if record.POS <= region.start + flank or record.POS + len(record.REF) + flank - 1 >= region.end: continue
            record.CHROM = str(region_index) if use_short_contig_names else ("%s_%d_%d" % (str(region.chrom), region.start, region.end))
            # record.POS seems to be zero-based, at least in the infinite wisdom of my version of pysam
            record.POS = record.POS - region.start
            if not new_samples:
                vcf_writer.write_record(record)
                continue
            else:
                snames = []
                sindexes = {}
            for s in new_samples:
                for i in xrange(len(record.samples)):
                    if s == record.samples[i].sample:
                        sindexes[s] = i
                        snames.append(record.samples[i])
            vcfrecord = vcf.model._Record(record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, record.FILTER, record.INFO, record.FORMAT, sindexes, snames)
            vcf_writer.write_record(vcfrecord)
    vcf_writer.close()
    pysam.tabix_index(out_vcf, force=True, preset='vcf')
    logger.info("Lifted over the VCF %s to %s" % (in_vcf, out_vcf))

    return "{}.gz".format(out_vcf)


def gen_restricted_ref_and_vcfs(reference, invcfs, regions, samples, outdir, flank=0, short_contig_names=False):
    restricted_fasta = reference
    outvcfs = invcfs

    if regions:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        restricted_fasta = os.path.join(outdir, os.path.basename(reference))
        gen_restricted_reference(reference, regions, restricted_fasta, short_contig_names)

        if outvcfs:
            outvcfs = map(lambda x: os.path.join(outdir, os.path.splitext(os.path.basename(x))[0]) if x else None, invcfs)
            generated_vcfs = []
            for invcf, outvcf in zip(invcfs, outvcfs):
                generated_vcfs.append(gen_restricted_vcf(invcf, regions, outvcf, restricted_fasta, samples, flank, short_contig_names))
            outvcfs = generated_vcfs

    return (restricted_fasta, outvcfs)


def main():
    logger = logging.getLogger(main.__name__)

    parser = argparse.ArgumentParser(description="Generate restricted FASTAs and VCFs given a BED file. The contigs are the sequences for each genomic region in the BED file and the name of the contigs reflects that. The VCFs use the coordinates on the new contigs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--reference", help="Reference FASTA", required=True)
    parser.add_argument("--regions", help="Regions BED", required=True)
    parser.add_argument("--vcfs", nargs="+", required=True, default=[])
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--flank", type=uint, default=0, help="Ignore variants this close to the edges of a region")
    parser.add_argument("--short_contig_names", action="store_true", help="Generate short contig names instead of the chr_start_end naming")
    parser.add_argument("--samples", nargs="+", default=[], help="Select specific samples. Select all samples if leave empty")

    args = parser.parse_args()

    gen_restricted_ref_and_vcfs(args.reference, args.vcfs, args.regions, args.samples, args.outdir, args.flank, args.short_contig_names)


if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    
    main()
