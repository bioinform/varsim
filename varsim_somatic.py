#!/usr/bin/python

# this is a temporary attempt to get a workable somatic workflow

import argparse
import os
import sys
import subprocess
import logging
import time
from varsim import VERSION, MY_DIR, VARSIMJAR, DEFAULT_VARSIMJAR, REQUIRE_VARSIMJAR
from varsim import check_java, makedirs, monitor_processes, check_executable, run_vcfstats, run_randvcf

VARSIM_PY = os.path.join(MY_DIR, "varsim.py")

if __name__ == "__main__":

    main_parser = argparse.ArgumentParser(description="VarSim: somatic workflow",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("--out_dir", metavar="Out directory", help="Output directory",
                             default="somatic_out")
    main_parser.add_argument("--work_dir", metavar="Work directory", help="Work directory",
                             default="somatic_work")
    main_parser.add_argument("--log_dir", metavar="Log directory", help="Directory to log to",
                             default="somatic_log")
    main_parser.add_argument("--reference", metavar="FASTA", help="Reference genome", required=True, type=file)
    main_parser.add_argument("--seed", metavar="INT", help="Random number seed", type=int, default=0)
    main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (MALE/FEMALE)", required=False, type=str,
                             choices=["MALE", "FEMALE"], default="MALE")
    main_parser.add_argument("--id", metavar="id", help="Sample ID", required=True)
    main_parser.add_argument("--simulator", metavar="simulator", help="Read simulator to use", required=False, type=str,
                             choices=["art", "dwgsim"], default="art")
    main_parser.add_argument("--simulator_executable", metavar="PATH",
                             help="Path to the executable of the read simulator chosen"
                             , required=True, type=file)
    main_parser.add_argument("--varsim_jar", metavar="PATH", help="Path to VarSim.jar (deprecated)", type=file,
                             default=DEFAULT_VARSIMJAR,
                             required=False)
    main_parser.add_argument("--read_length", metavar="INT", help="Length of read to simulate", default=100, type=int)
    main_parser.add_argument("--nlanes", metavar="INT",
                             help="Number of lanes to generate, coverage will be divided evenly over the lanes. Simulation is parallized over lanes. Each lane will have its own pair of files",
                             default=3, type=int)
    main_parser.add_argument("--total_coverage", metavar="FLOAT", help="Total coverage to simulate", default=1.0,
                             type=float)
    main_parser.add_argument("--mean_fragment_size", metavar="INT", help="Mean fragment size", default=350,
                             type=int)
    main_parser.add_argument("--sd_fragment_size", metavar="INT", help="Standard deviation of fragment size",
                             default=50, type=int)

    main_parser.add_argument("--force_five_base_encoding", action="store_true", help="Force bases to be ACTGN")
    main_parser.add_argument("--filter", action="store_true", help="Only use PASS variants")
    main_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
    main_parser.add_argument('--version', action='version', version='VarSim: %(prog)s ' + VERSION)


    input_vcf_group = main_parser.add_argument_group("Input VCFs options")
    input_vcf_group.add_argument("--cosmic_vcf", metavar="VCF", help="COSMIC database VCF. Need to specify when random COSMIC sampling is enabled.")
    input_vcf_group.add_argument("--normal_vcf", metavar="VCF", help="Normal VCF from previous VarSim run", required=True)
    input_vcf_group.add_argument("--somatic_vcfs", metavar="VCF", nargs="+", help="Somatic VCF", default=[])
    input_vcf_group.add_argument("--merge_priority", choices=["sn", "ns"], help="Priority of merging (lowest first) somatic (s) and normal truth (n).", default="sn")

    pipeline_control_group = main_parser.add_argument_group("Pipeline control options. Disable parts of the pipeline.")
    pipeline_control_group.add_argument("--disable_rand_vcf", action="store_true", help="Disable RandVCF2VCF somatic")
    pipeline_control_group.add_argument("--disable_vcf2diploid", action="store_true", help="Disable vcf2diploid")
    pipeline_control_group.add_argument("--disable_sim", action="store_true", help="Disable read simulation")

    # RandVCF2VCF seed num_SNP num_INS num_DEL num_MNP num_COMPLEX percent_novel min_length_lim max_length_lim reference_file file.vcf
    rand_vcf_group = main_parser.add_argument_group("RandVCF2VCF somatic options")
    rand_vcf_group.add_argument("--som_num_snp", metavar="INT", help="Number of somatic SNPs", default=9000, type=int)
    rand_vcf_group.add_argument("--som_num_ins", metavar="INT", help="Number of somatic insertions", default=1000,
                                type=int)
    rand_vcf_group.add_argument("--som_num_del", metavar="INT", help="Number of somatic deletions", default=1000,
                                type=int)
    rand_vcf_group.add_argument("--som_num_mnp", metavar="INT", help="Number of somatic MNPs", default=100, type=int)
    rand_vcf_group.add_argument("--som_num_complex", metavar="INT", help="Number of somatic complex variants",
                                default=100, type=int)
    # rand_vcf_group.add_argument("--som_percent_novel", metavar="percent_novel", help="Percent novel", default=0, type=float)
    rand_vcf_group.add_argument("--som_min_length_lim", metavar="INT", help="Min length lim", default=0,
                                type=int)
    rand_vcf_group.add_argument("--som_max_length_lim", metavar="INT", help="Max length lim", default=49,
                                type=int)
    # rand_vcf_group.add_argument("--som_vcf", metavar="in_vcf", help="Input somatic variant database VCF", type=file, required=False)
    rand_vcf_group.add_argument("--som_prop_het", metavar="FLOAT", help="Proportion of somatic heterozygous variants",
                                default=1.0, type=float)

    dwgsim_group = main_parser.add_argument_group("DWGSIM options")
    dwgsim_group.add_argument("--dwgsim_start_e", metavar="first_base_error_rate", help="Error rate on the first base",
                              default=0.0001, type=float)
    dwgsim_group.add_argument("--dwgsim_end_e", metavar="last_base_error_rate", help="Error rate on the last base",
                              default=0.0015, type=float)
    dwgsim_group.add_argument("--dwgsim_options", help="DWGSIM command-line options", default="", required=False)

    art_group = main_parser.add_argument_group("ART options")
    art_group.add_argument("--profile_1", metavar="profile_file1", help="Profile for first end", default=None, type=file)
    art_group.add_argument("--profile_2", metavar="profile_file2", help="Profile for second end", default=None, type=file)
    art_group.add_argument("--art_options", help="ART command-line options", default="", required=False)

    args = main_parser.parse_args()

    makedirs([args.log_dir, args.out_dir])

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(filename=os.path.join(args.log_dir, "varsim.log"), filemode="w", level=logging.DEBUG, format=FORMAT)
    logger = logging.getLogger(__name__)

    check_java()

    if not args.disable_sim:
        if not args.simulator_executable:
            logger.error("Please specify %s binary with --simulator_executable option" % args.simulator)
            sys.exit(os.EX_USAGE)
        check_executable(args.simulator_executable.name)

    t_s = time.time()

    cosmic_sampled_vcfs = []
    if not args.disable_rand_vcf:
        if not args.cosmic_vcf:
            logger.error("COSMIC database VCF not specified using --cosmic_vcf")
            sys.exit(os.EX_USAGE)
        rand_vcf_stdout = open(os.path.join(args.out_dir, "random.cosmic.vcf"), "w")
        rand_vcf_stderr = open(os.path.join(args.log_dir, "random.cosmic.err"), "w")
        cosmic_sampled_vcfs = [rand_vcf_stdout.name]

        # Not able to support novel yet for COSMIC variants
        monitor_processes([run_randvcf(os.path.realpath(args.cosmic_vcf), rand_vcf_stdout, rand_vcf_stderr,
                         args.seed, args.sex, args.som_num_snp, args.som_num_ins, args.som_num_del, args.som_num_mnp,
                         args.som_num_complex, 0, args.som_min_length_lim, args.som_max_length_lim,
                         os.path.realpath(args.reference.name), args.som_prop_het)])

    normal_vcfs = [args.normal_vcf]
    somatic_vcfs = cosmic_sampled_vcfs + args.somatic_vcfs
    fixed_somatic_vcfs = []
    if somatic_vcfs:
        vcfs_dir = os.path.join(args.out_dir, "somatic_vcfs")
        makedirs([vcfs_dir])
        count = 0
        for index, vcf in enumerate(somatic_vcfs):
            copied_vcf = os.path.join(vcfs_dir, "%d.vcf" % index)
            logger.info("Copying somatic VCF %s to %s and adding VARSIMSOMATIC id to entries if missing" % (vcf, copied_vcf))
            with open(vcf, "r") as vcf_fd, open(copied_vcf, "w") as copied_vcf_fd:
                for line in vcf_fd:
                    if line.startswith("#"):
                        copied_vcf_fd.write(line)
                    else:
                        line_fields = line.split("\t")
                        line_fields[2] = ("VARSIMSOMATIC%d" % count) if line_fields[2] == "." else ("%s,VARSIMSOMATIC%d" % (line_fields[2], count))
                        copied_vcf_fd.write("\t".join(line_fields))
                        count += 1
            fixed_somatic_vcfs.append(copied_vcf)

    vcf_files = (fixed_somatic_vcfs + normal_vcfs) if args.merge_priority == "sn" else (normal_vcfs + fixed_somatic_vcfs)
    vcf_files = map(os.path.realpath, filter(None, vcf_files))

    processes = run_vcfstats(vcf_files, args.out_dir, args.log_dir)

    # Run VarSim
    varsim_stdout = open(os.path.join(args.log_dir, "som_varsim.out"), "w")
    varsim_stderr = open(os.path.join(args.log_dir, "som_varsim.log"), "w")

    vcf_arg_list = ["--vcfs"] + vcf_files

    # need to fix the store true ones
    filter_arg_list = ["--filter"] if args.filter else []
    disable_sim_arg_list = ["--disable_sim"] if args.disable_sim else []
    force_five_base_encoding_arg_list = ["--force_five_base_encoding"] if args.force_five_base_encoding else []
    keep_temp_arg_list = ["--keep_temp"] if args.keep_temp else []
    profile_1_arg_list = ["--profile_1", args.profile_1.name] if args.profile_1 is not None else []
    profile_2_arg_list = ["--profile_2", args.profile_2.name] if args.profile_2 is not None else []
    other_varsim_opts = []
    if args.simulator == "dwgsim":
        other_varsim_opts = ["--dwgsim_start_e", str(args.dwgsim_start_e), "--dwgsim_end_e", str(args.dwgsim_end_e)]
        if args.dwgsim_options: other_varsim_opts += ["--dwgsim_options", str(args.dwgsim_options)]
    elif args.simulator == "art" and args.art_options:
        other_varsim_opts += ["--art_options", args.art_options]

    varsim_command = ["python", os.path.realpath(VARSIM_PY),
                      "--out_dir", str(os.path.realpath(args.out_dir)),
                      "--work_dir", str(os.path.realpath(args.work_dir)),
                      "--log_dir", str(os.path.realpath(os.path.join(args.log_dir, "varsim"))),
                      "--reference", str(os.path.realpath(args.reference.name)),
                      "--seed", str(args.seed),
                      "--sex", str(args.sex),
                      "--id", str(args.id),
                      "--simulator", str(args.simulator),
                      "--simulator_executable", str(args.simulator_executable.name),
                      "--read_length", str(args.read_length),
                      "--nlanes", str(args.nlanes),
                      "--total_coverage", str(args.total_coverage),
                      "--mean_fragment_size", str(args.mean_fragment_size),
                      "--sd_fragment_size", str(args.sd_fragment_size),
                      "--disable_rand_vcf",
                      "--disable_rand_dgv"] + other_varsim_opts + vcf_arg_list + filter_arg_list + disable_sim_arg_list \
                     + force_five_base_encoding_arg_list + keep_temp_arg_list + profile_1_arg_list + profile_2_arg_list
    varsim_command = " ".join(varsim_command)
    p_varsim = subprocess.Popen(varsim_command, stdout=varsim_stdout, stderr=varsim_stderr, shell=True)
    logger.info("Executing command " + varsim_command + " with pid " + str(p_varsim.pid))
    processes.append(p_varsim)

    processes = monitor_processes(processes)

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
            if line.find("VARSIMSOMATIC") >= 0:
                somatic_vcf_fd.write(line)
            else:
                normal_vcf_fd.write(line)

    monitor_processes(run_vcfstats([normal_vcf, somatic_vcf], args.out_dir, args.log_dir))

    logger.info("Done! (%g hours)" % ((time.time() - t_s) / 3600.0))
