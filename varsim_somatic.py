#!/usr/bin/python

# this is a temporary attempt to get a workable somatic workflow

import argparse
import os
import sys
import re
import subprocess
import logging
import time
import signal


def get_contigs_list(reference):
    with open("%s.fai" % (reference)) as fai_file:
        contigs = [line.strip().split()[0] for line in fai_file.readlines()]
    return contigs


def makedirs(dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)


my_dir = os.path.dirname(os.path.realpath(__file__))

default_varsim_jar = os.path.join(my_dir, "VarSim.jar")
require_varsim_jar = not os.path.isfile(default_varsim_jar)

if not os.path.isfile(default_varsim_jar):
    default_varsim_jar = None

default_varsim = os.path.join(my_dir, "varsim.py")
require_varsim = not os.path.isfile(default_varsim)

if not os.path.isfile(default_varsim):
    default_varsim = None

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
main_parser.add_argument("--varsim_jar", metavar="PATH", help="Path to VarSim.jar", type=file,
                         default=default_varsim_jar,
                         required=require_varsim_jar)
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
main_parser.add_argument("--cosmic_vcf", metavar="VCF", help="COSMIC database VCF")
main_parser.add_argument("--normal_vcf", metavar="VCF", help="Normal VCF from previous VarSim run", required=True)
main_parser.add_argument("--somatic_vcfs", metavar="VCF", nargs="+", help="Somatic VCF. Make sure that the VCF entries have the SOMATIC flag", default=[])
main_parser.add_argument("--force_five_base_encoding", action="store_true", help="Force bases to be ACTGN")
main_parser.add_argument("--filter", action="store_true", help="Only use PASS variants")
main_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
main_parser.add_argument("--varsim_py", metavar="PATH", help="Path to VarSim.py", type=file,
                         default=default_varsim, required=require_varsim)

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


def run_shell_command(cmd, cmd_stdout, cmd_stderr, cmd_dir="."):
    subproc = subprocess.Popen(cmd, stdout=cmd_stdout, stderr=cmd_stderr, cwd=cmd_dir, shell=True, preexec_fn=os.setsid)
    retcode = subproc.wait()
    sys.exit(retcode)


def monitor_processes(processes):
    logger = logging.getLogger(monitor_processes.__name__)
    while processes:
        time.sleep(1)
        kill_all = False
        processes_running = []
        for p in processes:
            status = p.poll()
            if status is not None:
                logger.info("Process %s exited with code %d" % (p.pid, status))
                if status != 0:
                    kill_all = True
                    logger.error("Process %s failed. Will kill the remaining processes." % (p.pid))
            else:
                processes_running.append(p)
        if kill_all:
            for p in processes:
                status = p.poll()
                if status is None:
                    try:
                        os.killpg(p.pid, signal.SIGTERM)
                    except OSError, e:
                        try:
                            p.terminate()
                        except OSError, ex:
                            logger.write("Could not kill the process %d\n" % p.pid)
            sys.exit(1)
        processes = processes_running
    return []


def check_executable(fpath):
    logger = logging.getLogger(check_executable.__name__)
    if not os.path.isfile(fpath):
        logger.error("File %s does not exist" % fpath)
        sys.exit(os.EX_NOINPUT)
    if not os.access(fpath, os.X_OK):
        logger.error("File %s is not executable" % fpath)
        sys.exit(os.EX_NOINPUT)


def run_vcfstats(vcfs, varsim_jar, out_dir, log_dir):
    logger = logging.getLogger(run_vcfstats.__name__)
    processes = []
    for in_vcf in vcfs:
        out_prefix = os.path.basename(in_vcf)
        vcfstats_stdout = open(os.path.join(out_dir, "%s.stats" % (out_prefix)), "w")
        vcfstats_stderr = open(os.path.join(log_dir, "%s.vcfstats.err" % (out_prefix)), "w")
        vcfstats_command = ["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(varsim_jar), "vcfstats", "-vcf",
                        in_vcf]
        p_vcfstats = subprocess.Popen(vcfstats_command, stdout=vcfstats_stdout, stderr=vcfstats_stderr)
        logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
        processes.append(p_vcfstats)
    return processes


if not args.disable_sim:
    if not args.simulator_executable:
        logger.error("Please specify %s binary with --simulator_executable option" % args.simulator)
        sys.exit(os.EX_USAGE)
    check_executable(args.simulator_executable.name)

processes = []

t_s = time.time()

cosmic_sampled_vcf = []
if not args.disable_rand_vcf:
    if not args.cosmic_vcf:
        raise Exception("COSMIC database VCF not specified using --cosmic_vcf")
    rand_vcf_stdout = open(os.path.join(args.out_dir, "random.cosmic.vcf"), "w")
    rand_vcf_stderr = open(os.path.join(args.log_dir, "random.cosmic.err"), "w")
    cosmic_sampled_vcf = [rand_vcf_stdout.name]

    rand_vcf_command = ["java", "-jar", os.path.realpath(args.varsim_jar.name), "randvcf2vcf", "-seed", str(args.seed),
                        "-num_snp", str(args.som_num_snp),
                        "-num_ins", str(args.som_num_ins),
                        "-num_del", str(args.som_num_del),
                        "-num_mnp", str(args.som_num_mnp),
                        "-num_complex", str(args.som_num_complex),
                        # "-novel", str(args.som_percent_novel), # Not able to support novel yet (COS)
                        "-min_len", str(args.som_min_length_lim),
                        "-max_len", str(args.som_max_length_lim),
                        "-ref", os.path.realpath(args.reference.name),
                        "-prop_het", str(args.som_prop_het),
                        "-vcf", os.path.realpath(args.cosmic_vcf)]

    p_rand_vcf = subprocess.Popen(rand_vcf_command, stdout=rand_vcf_stdout, stderr=rand_vcf_stderr)
    logger.info("Executing command " + " ".join(rand_vcf_command) + " with pid " + str(p_rand_vcf.pid))
    processes.append(p_rand_vcf)

processes = monitor_processes(processes)

somatic_vcfs = cosmic_sampled_vcf + args.somatic_vcfs
fixed_somatic_vcfs = []
if somatic_vcfs:
    somatic_vcfs_dir = os.path.join(args.out_dir, "somatic_vcfs")
    makedirs([somatic_vcfs_dir])
    logger.info("Copying somatic VCFs over and adding SOMATIC flag to entries if missing")
    for index, vcf in enumerate(somatic_vcfs):
        copied_vcf = os.path.join(somatic_vcfs_dir, "%d.vcf" % index)
        with open(vcf, "r") as vcf_fd, open(copied_vcf, "w") as copied_vcf_fd:
            for line in vcf_fd:
                if line.startswith("#"):
                    copied_vcf_fd.write(line)
                else:
                    line_fields = line.split("\t")
                    if line_fields[7].find("SOMATIC") < 0:
                        line_fields[7] += ";SOMATIC"
                    copied_vcf_fd.write("\t".join(line_fields))
        fixed_somatic_vcfs.append(copied_vcf)
            

vcf_files = map(os.path.realpath, filter(None, fixed_somatic_vcfs + [args.normal_vcf]))

processes = run_vcfstats(vcf_files, args.varsim_jar.name, args.out_dir, args.log_dir)

# Run VarSim 
varsim_stdout = open(os.path.join(args.log_dir, "som_varsim.out"), "w")
varsim_stderr = open(os.path.join(args.log_dir, "som_varsim.log"), "w")

vcf_arg_list = ["--vcfs"] + vcf_files

# need to fix the store true ones
filter_arg_list = ["--filter"] if args.filter else []
disable_sim_arg_list = ["--disable_sim"] if args.disable_sim else []
force_five_base_encoding_arg_list = ["--force_five_base_encoding"] if args.force_five_base_encoding else []
keep_temp_arg_list = ["--keep_temp"] if args.keep_temp else []
profile_1_arg_list = ["--profile_1", args.profile_1] if args.profile_1 is not None else []
profile_2_arg_list = ["--profile_2", args.profile_2] if args.profile_2 is not None else []
varsim_command = ["python", os.path.realpath(args.varsim_py.name),
                  "--out_dir", str(os.path.realpath(args.out_dir)),
                  "--work_dir", str(os.path.realpath(args.work_dir)),
                  "--log_dir", str(os.path.realpath(args.log_dir)),
                  "--reference", str(os.path.realpath(args.reference.name)),
                  "--seed", str(args.seed),
                  "--sex", str(args.sex),
                  "--id", str(args.id),
                  "--simulator", str(args.simulator),
                  "--simulator_executable", str(args.simulator_executable.name),
                  "--varsim_jar", str(os.path.realpath(args.varsim_jar.name)),
                  "--read_length", str(args.read_length),
                  "--nlanes", str(args.nlanes),
                  "--total_coverage", str(args.total_coverage),
                  "--mean_fragment_size", str(args.mean_fragment_size),
                  "--sd_fragment_size", str(args.sd_fragment_size),
                  "--dwgsim_start_e", str(args.dwgsim_start_e),
                  "--dwgsim_end_e", str(args.dwgsim_end_e),
                  "--dwgsim_options", str(args.dwgsim_options),
                  "--art_options", str(args.art_options),
                  "--disable_rand_vcf",
                  "--disable_rand_dgv"] + vcf_arg_list + filter_arg_list + disable_sim_arg_list \
                 + force_five_base_encoding_arg_list + keep_temp_arg_list + profile_1_arg_list + profile_2_arg_list
p_varsim = subprocess.Popen(varsim_command, stdout=varsim_stdout, stderr=varsim_stderr)
logger.info("Executing command " + " ".join(varsim_command) + " with pid " + str(p_varsim.pid))
processes.append(p_varsim)

processes = monitor_processes(processes)

# Split the tumor truth VCF into normal variants and somatic variants
tumor_vcf = os.path.realpath(os.path.join(args.out_dir, "%s.truth.vcf" % args.id))
normal_vcf = os.path.join(args.out_dir, "%s_norm.vcf" % args.id)
somatic_vcf = os.path.join(args.out_dir, "%s_somatic.vcf" % args.id)
logger.info("Splitting the truth VCF into normal and somatic VCFs")
with open(tumor_vcf, "r") as tumor_truth_fd, \
    open(normal_vcf, "w") as normal_vcf_fd, \
    open(somatic_vcf, "w") as somatic_vcf_fd:
    for line in tumor_truth_fd:
        if line.startswith("#"):
            somatic_vcf_fd.write(line)
            normal_vcf_fd.write(line)
            continue
        if line.find("SOMATIC") >= 0:
            somatic_vcf_fd.write(line)
        else:
            normal_vcf_fd.write(line)

monitor_processes(run_vcfstats([normal_vcf, somatic_vcf], args.varsim_jar.name, args.out_dir, args.log_dir))

logger.info("Done! (%g hours)" % ((time.time() - t_s) / 3600.0))
