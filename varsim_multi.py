#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import logging
import shutil
import time
import signal
import itertools
import glob
import tempfile
import re
from distutils.version import LooseVersion
from liftover_restricted_vcf_map import lift_vcfs, lift_maps
from generate_small_test_ref import gen_restricted_ref_and_vcfs 
from varsim import varsim_main, get_version, check_java, get_loglevel, makedirs, RandVCFOptions, RandDGVOptions, run_randvcf
import pybedtools

MY_DIR = os.path.dirname(os.path.realpath(__file__))
VARSIMJAR = os.path.realpath(os.path.join(MY_DIR, "VarSim.jar"))
DEFAULT_VARSIMJAR = os.path.join(MY_DIR, "VarSim.jar")
REQUIRE_VARSIMJAR = not os.path.isfile(DEFAULT_VARSIMJAR)
if REQUIRE_VARSIMJAR: DEFAULT_VARSIMJAR = None


def varsim_multi(reference,
                 simulator,
                 simulator_exe,
                 total_coverage,
                 variant_vcfs=[],
                 sampling_vcf=None,
                 dgv_file=None,
                 regions=None,
                 randvcf_options=None,
                 randdgv_options=None,
                 nlanes=1,
                 simulator_options="",
                 samples=[],
                 out_dir="out",
                 sv_insert_seq=None,
                 seed=0,
                 sex="MALE",
                 remove_filtered=False,
                 keep_temp=False,
                 force_five_base_encoding=False,
                 lift_ref=False,
                 disable_vcf2diploid=False):
    logger = logging.getLogger(varsim_multi.__name__)

    makedirs([out_dir])

    restricted_dir = os.path.join(out_dir, "restricted")

    restricted_reference, restricted_vcfs = gen_restricted_ref_and_vcfs(reference, variant_vcfs, regions, samples, restricted_dir, flank=0, short_contig_names=False)

    merged_bedtool = pybedtools.BedTool(regions).merge() if regions else None

    if regions and sampling_vcf:
        merged_bed = os.path.join(out_dir, "merged.bed")
        merged_bedtool = pybedtools.BedTool(regions).merge().saveas(merged_bed)
        _, [restricted_sampling_vcf] = gen_restricted_ref_and_vcfs(reference, [sampling_vcf], merged_bed, [], os.path.join(out_dir, "region_restricted"), flank=0)
        # Now lift over the restricted_sampling_vcf to get the region-limited VCF
        sampling_vcf = lift_vcfs([restricted_sampling_vcf], os.path.join(out_dir, "region_restricted", "region-restricted-sampling.vcf"), reference)
        

    #_, [restricted_sampling_vcf] = gen_restricted_ref_and_vcfs(reference, [sampling_vcf], regions, [], os.path.join(out_dir, "restricted_sampling"), flank=0, short_contig_names=False)

    for index, sample in enumerate(samples):
        sample_dir = os.path.join(out_dir, sample)
        logger.info("Simulating sample {} in {}".format(sample, sample_dir))

        # Run RandVCF first to get the sampled variants for the sample
        if randvcf_options:
            sampled_vcf = os.path.join(sample_dir, "randvcf.vcf")
            with open(sampled_vcf, "w") as randvcf_out, open(os.path.join(sample_dir, "randvcf.err"), "w") as randvcf_log:
                run_randvcf(sampling_vcf, randvcf_out, randvcf_log, seed + 1000 * index).wait()
            # Now generate the restricted sampled VCF for the sample
            _, [restricted_sampled_vcf] = gen_restricted_ref_and_vcfs(reference, [sampled_vcf], regions, [], os.path.join(sample_dir, "restricted_randvcf.vcf"), flank=0)
            sample_variant_vcfs = restricted_vcfs + [restricted_sampled_vcf]
        else:
            sample_variant_vcfs = restricted_vcfs

        varsim_main(restricted_reference,
                    simulator,
                    simulator_exe,
                    total_coverage,
                    restricted_vcfs,
                    None,
                    dgv_file,
                    None,
                    randdgv_options,
                    nlanes,
                    simulator_options,
                    sample,
                    os.path.join(sample_dir, "log"),
                    os.path.join(sample_dir, "out"),
                    sv_insert_seq,
                    seed + 1000 * index,
                    sex,
                    remove_filtered,
                    keep_temp,
                    force_five_base_encoding,
                    lift_ref,
                    disable_vcf2diploid)


if __name__ == "__main__":
    check_java()

    main_parser = argparse.ArgumentParser(description="VarSim: A high-fidelity simulation validation framework",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("--out_dir", metavar="DIR",
                             help="Output directory for the simulated genome, reads and variants", required=False,
                             default="out")
    main_parser.add_argument("--reference", metavar="FASTA", help="Reference genome that variants will be inserted into", required=True)
    main_parser.add_argument("--regions", help="Regions of interest for simulation. Skip for whole genome simulation")
    main_parser.add_argument("--seed", metavar="seed", help="Random number seed for reproducibility", type=int, default=0)
    main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (MALE/FEMALE)", required=False, type=str,
                             choices=["MALE", "FEMALE"], default="MALE")
    main_parser.add_argument("--samples", help="Samples to be simulated", required=True, nargs="+")
    main_parser.add_argument("--simulator", metavar="SIMULATOR", help="Read simulator to use", choices=["art", "dwgsim", "longislnd"], default="art")
    main_parser.add_argument("--simulator_executable", metavar="PATH",
                             help="Path to the executable of the read simulator chosen")
    main_parser.add_argument("--simulator_options", help="Simulator options", required=False, default="")
    main_parser.add_argument("--nlanes", metavar="INTEGER",
                             help="Number of lanes to generate, coverage will be divided evenly over the lanes. Simulation is parallized over lanes. Each lane will have its own pair of files",
                             default=1, type=int)
    main_parser.add_argument("--total_coverage", metavar="FLOAT", help="Total coverage to simulate", default=1.0,
                             type=float)
    main_parser.add_argument("--vcfs", metavar="VCF",
                             help="Addtional list of VCFs to insert into genome, priority is lowest ... highest", nargs="+",
                             default=[])
    main_parser.add_argument("--force_five_base_encoding", action="store_true", help="Force output bases to be only ACTGN")
    main_parser.add_argument("--filter", action="store_true", help="Only use PASS variants for simulation")
    main_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files after simulation")
    main_parser.add_argument("--lift_ref", action="store_true", help="Liftover chromosome names from restricted reference")
    main_parser.add_argument('--version', action='version', version=get_version())
    main_parser.add_argument('--log_to_stderr', action='store_true', help='Output log to stderr instead of log_dir/varsim.log')
    main_parser.add_argument("--loglevel", help="Set logging level", choices=["debug", "warn", "info"], default="info")

    pipeline_control_group = main_parser.add_argument_group("Pipeline control options. Disable parts of the pipeline.")
    pipeline_control_group.add_argument("--disable_rand_vcf", action="store_true",
                                        help="Disable sampling from the provided small variant VCF")
    pipeline_control_group.add_argument("--disable_rand_dgv", action="store_true",
                                        help="Disable sampline from the provided DGV file")
    pipeline_control_group.add_argument("--disable_vcf2diploid", action="store_true",
                                        help="Disable diploid genome simulation")
    pipeline_control_group.add_argument("--disable_sim", action="store_true", help="Disable read simulation")

    # RandVCF2VCF seed num_SNP num_INS num_DEL num_MNP num_COMPLEX percent_novel min_length_lim max_length_lim reference_file file.vcf
    rand_vcf_group = main_parser.add_argument_group("Small variant simulation options")
    rand_vcf_group.add_argument("--vc_num_snp", metavar="INTEGER", help="Number of SNPs to sample from small variant VCF",
                                default=0, type=int)
    rand_vcf_group.add_argument("--vc_num_ins", metavar="INTEGER",
                                help="Number of insertions to sample from small variant VCF", default=0, type=int)
    rand_vcf_group.add_argument("--vc_num_del", metavar="INTEGER",
                                help="Number of deletions to sample from small variant VCF", default=0, type=int)
    rand_vcf_group.add_argument("--vc_num_mnp", metavar="INTEGER", help="Number of MNPs to sample from small variant VCF",
                                default=0, type=int)
    rand_vcf_group.add_argument("--vc_num_complex", metavar="INTEGER",
                                help="Number of complex variants to sample from small variant VCF", default=0,
                                type=int)
    rand_vcf_group.add_argument("--vc_percent_novel", metavar="FLOAT",
                                help="Percent variants sampled from small variant VCF that will be moved to novel positions",
                                default=0, type=float)
    rand_vcf_group.add_argument("--vc_min_length_lim", metavar="INTEGER",
                                help="Min length of small variant to accept [inclusive]", default=0, type=int)
    rand_vcf_group.add_argument("--vc_max_length_lim", metavar="INTEGER",
                                help="Max length of small variant to accept [inclusive]", default=99,
                                type=int)
    rand_vcf_group.add_argument("--sampling_vcf", metavar="VCF", help="Input small variant sampling VCF, usually dbSNP")
    rand_vcf_group.add_argument("--vc_prop_het", metavar="FLOAT", help="Proportion of heterozygous small variants",
                                default=0.6,
                                type=float)

    # RandDGV2VCF seed num_INS num_DEL num_DUP num_INV percent_novel min_length_lim max_length_lim reference_file insert_seq.txt dgv_file.txt
    rand_dgv_group = main_parser.add_argument_group("Structural variant simulation options")
    rand_dgv_group.add_argument("--sv_num_ins", metavar="INTEGER", help="Number of insertions to sample from DGV",
                                default=20, type=int)
    rand_dgv_group.add_argument("--sv_num_del", metavar="INTEGER", help="Number of deletions to sample from DGV",
                                default=20, type=int)
    rand_dgv_group.add_argument("--sv_num_dup", metavar="INTEGER", help="Number of duplications to sample from DGV",
                                default=20, type=int)
    rand_dgv_group.add_argument("--sv_num_inv", metavar="INTEGER", help="Number of inversions to sample from DGV",
                                default=20, type=int)
    rand_dgv_group.add_argument("--sv_percent_novel", metavar="FLOAT",
                                help="Percent variants sampled from DGV that will be moved to novel positions", default=0,
                                type=float)
    rand_dgv_group.add_argument("--sv_min_length_lim", metavar="min_length_lim",
                                help="Min length of structural variant to accept [inclusive]", default=100,
                                type=int)
    rand_dgv_group.add_argument("--sv_max_length_lim", metavar="max_length_lim",
                                help="Max length of structural variant to accept [inclusive]", default=1000000,
                                type=int)
    rand_dgv_group.add_argument("--sv_insert_seq", metavar="FILE",
                                help="Path to file containing concatenation of real insertion sequences",
                                required=False)
    rand_dgv_group.add_argument("--sv_dgv", metavar="DGV_FILE", help="DGV file containing structural variants",
                                required=False)

    args = main_parser.parse_args()

    makedirs([args.out_dir])

    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = get_loglevel(args.loglevel)
    if not args.log_to_stderr:
        logging.basicConfig(filename=os.path.join(args.out_dir, "varsim.log"), filemode="w", level=loglevel, format=FORMAT)
    else:
        logging.basicConfig(level=loglevel, format=FORMAT)

    simulator = None if args.disable_sim else args.simulator
    randvcf_options = None if args.disable_rand_vcf else RandVCFOptions(args.vc_num_snp, args.vc_num_ins, args.vc_num_del, args.vc_num_mnp, args.vc_num_complex, args.vc_percent_novel, args.vc_min_length_lim, args.vc_max_length_lim, args.vc_prop_het)
    randdgv_options = None if args.disable_rand_dgv else RandDGVOptions(args.sv_num_ins, args.sv_num_del, args.sv_num_dup, args.sv_num_inv, args.sv_percent_novel, args.sv_min_length_lim, args.sv_max_length_lim)
    varsim_multi(args.reference,
                 simulator,
                 args.simulator_executable,
                 args.total_coverage,
                 variant_vcfs=args.vcfs,
                 sampling_vcf=args.sampling_vcf,
                 dgv_file=args.sv_dgv,
                 regions=args.regions,
                 randvcf_options=randvcf_options,
                 randdgv_options=randdgv_options,
                 nlanes=args.nlanes,
                 simulator_options=args.simulator_options,
                 samples=args.samples,
                 out_dir=args.out_dir,
                 sv_insert_seq=args.sv_insert_seq,
                 seed=args.seed,
                 sex=args.sex,
                 remove_filtered=args.filter,
                 keep_temp=args.keep_temp,
                 force_five_base_encoding=args.force_five_base_encoding,
                 lift_ref=args.lift_ref,
                 disable_vcf2diploid=args.disable_vcf2diploid)
