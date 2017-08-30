#!/usr/bin/env python

# this is a temporary attempt to get a workable somatic workflow

import argparse
import os
import sys
import logging
import time
from varsim import SUPPORTED_SIMULATORS
from varsim import check_java, makedirs, monitor_processes, run_vcfstats, run_randvcf, get_version, RandVCFOptions
from varsim import varsim_main, add_varsim_id


def varsim_somatic_main():

    check_java()

    main_parser = argparse.ArgumentParser(description="VarSim: somatic workflow",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("--out_dir", metavar="Out directory", help="Output directory",
                             default="somatic_out")
    main_parser.add_argument("--work_dir", metavar="Work directory", help="Work directory",
                             default="somatic_work")
    main_parser.add_argument("--log_dir", metavar="Log directory", help="Directory to log to",
                             default="somatic_log")
    main_parser.add_argument("--reference", metavar="FASTA", help="Reference genome", required=True)
    main_parser.add_argument("--seed", metavar="INT", help="Random number seed", type=int, default=0)
    main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (MALE/FEMALE)",
                             choices=["MALE", "FEMALE"], default="MALE")
    main_parser.add_argument("--id", metavar="id", help="Sample ID", required=True)
    main_parser.add_argument("--simulator", metavar="simulator", help="Read simulator to use",
                             choices=SUPPORTED_SIMULATORS, default="art")
    main_parser.add_argument("--simulator_executable", metavar="PATH",
                             help="Path to the executable of the read simulator chosen"
                             , required=True)
    main_parser.add_argument("--simulator_options", help="Options to be passed to read simulator", default="")
    main_parser.add_argument("--regions", help="Restrict simulation to regions in the BED file")
    main_parser.add_argument("--varsim_jar", help="Path to VarSim.jar (deprecated)")
    main_parser.add_argument("--read_length", help="Length of read to simulate (deprecated)", default=100, type=int)
    main_parser.add_argument("--nlanes", metavar="INT",
                             help="Number of lanes to generate, coverage will be divided evenly over the lanes. Simulation is parallized over lanes. Each lane will have its own pair of files",
                             default=3, type=int)
    main_parser.add_argument("--total_coverage", metavar="FLOAT", help="Total coverage to simulate", default=1.0,
                             type=float)
    main_parser.add_argument("--mean_fragment_size", help="Mean fragment size (deprecated)", default=350, type=int)
    main_parser.add_argument("--sd_fragment_size", help="Standard deviation of fragment size (deprecated)",
                             default=50, type=int)

    main_parser.add_argument("--force_five_base_encoding", action="store_true", help="Force bases to be ACTGN")
    main_parser.add_argument("--filter", action="store_true", help="Only use PASS variants")
    main_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
    main_parser.add_argument('--version', action='version', version=get_version())

    input_vcf_group = main_parser.add_argument_group("Input VCFs options")
    input_vcf_group.add_argument("--cosmic_vcf", metavar="VCF", help="COSMIC database VCF. Need to specify when random COSMIC sampling is enabled. (deprecated)")
    input_vcf_group.add_argument("--somatic_sampling_vcf", help="Somatic database VCF for sampling.")
    input_vcf_group.add_argument("--normal_vcf", metavar="VCF", help="Normal VCF from previous VarSim run", required=True)
    input_vcf_group.add_argument("--somatic_vcfs", metavar="VCF", nargs="+", help="Somatic VCF", default=[])
    input_vcf_group.add_argument("--merge_priority", choices=["sn", "ns"], help="Priority of merging (lowest first) somatic (s) and normal truth (n).", default="sn")

    pipeline_control_group = main_parser.add_argument_group("Pipeline control options. Disable parts of the pipeline.")
    pipeline_control_group.add_argument("--disable_rand_vcf", action="store_true", help="Disable RandVCF2VCF somatic")
    pipeline_control_group.add_argument("--disable_vcf2diploid", action="store_true", help="Disable vcf2diploid")
    pipeline_control_group.add_argument("--disable_sim", action="store_true", help="Disable read simulation")

    # RandVCF2VCF seed num_SNP num_INS num_DEL num_MNP num_COMPLEX percent_novel min_length_lim max_length_lim reference_file file.vcf
    rand_vcf_group = main_parser.add_argument_group("RandVCF2VCF somatic options")
    rand_vcf_group.add_argument("--som_num_snp", help="Number of somatic SNPs", default=9000, type=int)
    rand_vcf_group.add_argument("--som_num_ins", help="Number of somatic insertions", default=1000, type=int)
    rand_vcf_group.add_argument("--som_num_del", help="Number of somatic deletions", default=1000, type=int)
    rand_vcf_group.add_argument("--som_num_mnp", help="Number of somatic MNPs", default=100, type=int)
    rand_vcf_group.add_argument("--som_num_complex", help="Number of somatic complex variants", default=100, type=int)
    rand_vcf_group.add_argument("--som_min_length_lim", help="Min length lim", default=0, type=int)
    rand_vcf_group.add_argument("--som_max_length_lim", help="Max length lim", default=49, type=int)
    rand_vcf_group.add_argument("--som_prop_het", help="Proportion of somatic heterozygous variants", default=1.0, type=float)
    rand_vcf_group.add_argument("--sv_insert_seq",
                                help="Path to file containing concatenation of real insertion sequences (deprecated)", type=file,
                                required=False)

    dwgsim_group = main_parser.add_argument_group("DWGSIM options (deprecated)")
    dwgsim_group.add_argument("--dwgsim_start_e", help="Error rate on the first base (deprecated)", default=0.0001, type=float)
    dwgsim_group.add_argument("--dwgsim_end_e", help="Error rate on the last base (deprecated)", default=0.0015, type=float)
    dwgsim_group.add_argument("--dwgsim_options", help="DWGSIM command-line options (deprecated)", default="")

    art_group = main_parser.add_argument_group("ART options (deprecated)")
    art_group.add_argument("--profile_1", help="Profile for first end (deprecated)")
    art_group.add_argument("--profile_2", help="Profile for second end (deprecated)")
    art_group.add_argument("--art_options", help="ART command-line options (deprecated)", default="")

    args = main_parser.parse_args()

    makedirs([args.log_dir, args.out_dir])

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(filename=os.path.join(args.log_dir, "varsim.log"), filemode="w", level=logging.DEBUG, format=FORMAT)
    logger = logging.getLogger(varsim_somatic_main.__name__)

    if not args.disable_sim:
        if not args.simulator_executable:
            logger.error("Please specify %s binary with --simulator_executable option" % args.simulator)
            sys.exit(os.EX_USAGE)

    t_s = time.time()

    somatic_sampled_vcfs = []
    if not args.disable_rand_vcf:
        if not args.somatic_sampling_vcf:
            logger.error("COSMIC database VCF not specified using --somatic_sampling_vcf")
            sys.exit(os.EX_USAGE)
        randvcf_options = RandVCFOptions(args.som_num_snp, args.som_num_ins, args.som_num_del, args.som_num_mnp, args.som_num_complex, 0, args.som_min_length_lim, args.som_max_length_lim, args.som_prop_het)

        with open(os.path.join(args.out_dir, "random.somatic.vcf"), "w") as rand_vcf_stdout, \
                open(os.path.join(args.log_dir, "random.somatic.err"), "w") as rand_vcf_stderr:
            somatic_sampled_vcfs = [rand_vcf_stdout.name]

            # Not able to support novel yet for somatic variants
            monitor_processes([run_randvcf(os.path.realpath(args.somatic_sampling_vcf), rand_vcf_stdout, rand_vcf_stderr, args.seed, args.sex, randvcf_options, args.reference)])

    normal_vcfs = [args.normal_vcf]
    somatic_vcfs = somatic_sampled_vcfs + args.somatic_vcfs
    fixed_somatic_vcfs = []
    if somatic_vcfs:
        vcfs_dir = os.path.join(args.out_dir, "somatic_vcfs")
        makedirs([vcfs_dir])
        for index, src_vcf in enumerate(somatic_vcfs):
            dst_vcf = os.path.join(vcfs_dir, "%d.vcf" % index)
            add_varsim_id(src_vcf, dst_vcf, sindex=0, id_prefix="VarSimSomatic%d" % index)
            fixed_somatic_vcfs.append(dst_vcf)

    vcf_files = (fixed_somatic_vcfs + normal_vcfs) if args.merge_priority == "sn" else (normal_vcfs + fixed_somatic_vcfs)
    vcf_files = map(os.path.realpath, filter(None, vcf_files))

    monitor_processes(run_vcfstats(vcf_files, args.out_dir, args.log_dir))

    varsim_main(args.reference,
                args.simulator if not args.disable_sim else None,
                args.simulator_executable,
                args.total_coverage,
                variant_vcfs=vcf_files,
                sampling_vcf=None,
                dgv_file=None,
                randvcf_options=None,
                randdgv_options=None,
                nlanes=args.nlanes,
                simulator_options=args.simulator_options,
                sample_id=args.id,
                log_dir=args.log_dir,
                out_dir=args.work_dir,
                sv_insert_seq=None,
                seed=args.seed,
                sex=args.sex,
                remove_filtered=args.filter,
                keep_temp=args.keep_temp,
                force_five_base_encoding=args.force_five_base_encoding,
                lift_ref=args.lift_ref,
                disable_vcf2diploid=args.disable_vcf2diploid)

    # Split the tumor truth VCF into normal variants and somatic variants
    tumor_vcf = os.path.realpath(os.path.join(args.out_dir, "%s.truth.vcf" % args.id))
    normal_vcf = os.path.join(args.out_dir, "%s_norm.vcf" % args.id)
    somatic_vcf = os.path.join(args.out_dir, "%s_somatic.vcf" % args.id)
    logger.info("Splitting the truth VCF %s into normal and somatic VCFs" % tumor_vcf)
    with open(tumor_vcf, "r") as tumor_truth_fd, \
        open(normal_vcf, "w") as normal_vcf_fd, \
        open(somatic_vcf, "w") as somatic_vcf_fd:
        for line in tumor_truth_fd:
            if line.startswith("#"):
                somatic_vcf_fd.write(line)
                normal_vcf_fd.write(line)
                continue
            if line.find("VarSimSomatic") >= 0:
                somatic_vcf_fd.write(line)
            else:
                normal_vcf_fd.write(line)

    monitor_processes(run_vcfstats([normal_vcf, somatic_vcf], args.out_dir, args.log_dir))

    logger.info("Done! (%g hours)" % ((time.time() - t_s) / 3600.0))

if __name__ == "__main__":
    varsim_somatic_main()
