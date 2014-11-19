#!/usr/bin/python

# this is a temporary attempt to get a workable somatic workflow

import argparse
import os
import sys
import subprocess
import logging
import shutil
import time
import signal
import multiprocessing as mp
from multiprocessing import Process

def get_contigs_list(reference):
  with open("%s.fai" % (reference)) as fai_file:
    contigs = [line.strip().split()[0] for line in fai_file.readlines()]
  return contigs

my_dir = os.path.dirname(os.path.realpath(__file__))

default_vcf2diploid = os.path.join(my_dir, "target/build/vcf2diploid.jar")
default_randvcf2vcf = os.path.join(my_dir, "target/build/randvcf2vcf.jar")
default_liftover    = os.path.join(my_dir, "target/build/fastq_liftover.jar")
default_vcfstats    = os.path.join(my_dir, "target/build/vcfstats.jar")
default_varsim      = os.path.join(my_dir, "varsim.py")

require_randvcf2vcf = not os.path.isfile(default_randvcf2vcf)
require_vcf2diploid = not os.path.isfile(default_vcf2diploid)
require_liftover    = not os.path.isfile(default_liftover)
require_vcfstats    = not os.path.isfile(default_vcfstats)
require_varsim      = not os.path.isfile(default_varsim)

if not os.path.isfile(default_randvcf2vcf): default_randvcf2vcf = None
if not os.path.isfile(default_vcf2diploid): default_vcf2diploid = None
if not os.path.isfile(default_liftover):    default_liftover    = None
if not os.path.isfile(default_vcfstats):    default_vcfstats    = None
if not os.path.isfile(default_varsim):      require_varsim      = None

main_parser = argparse.ArgumentParser(description="VarSim: somatic workflow", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
main_parser.add_argument("--out_dir", metavar="Out directory", help="Output directory", required=False, default="somatic_out")
main_parser.add_argument("--work_dir", metavar="Work directory", help="Work directory", required=False, default="somatic_work")
main_parser.add_argument("--log_dir", metavar="Log directory", help="Directory to log to", required=False, default="somatic_log")
main_parser.add_argument("--reference", metavar="Reference", help="Reference file", required=True, type=file)
main_parser.add_argument("--seed", metavar="seed", help="Random number seed", type=int, default=0)
main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (male/female)", required=False, type=str, choices=["male", "female"], default="male")
main_parser.add_argument("--id", metavar="id", help="Sample ID", required=True)
main_parser.add_argument("--simulator", metavar="simulator", help="Read simulator", required=False, type=str, choices=["art", "dwgsim"], default="art")
main_parser.add_argument("--vcf2diploid_jar", metavar="vcf2diploid_jar", help="vcf2diploid jar", type=file, default=default_vcf2diploid, required=require_vcf2diploid)
main_parser.add_argument("--read_length",    metavar="read_length",    help="Length of reads", default=100, type=int)
main_parser.add_argument("--nlanes",         metavar="nlanes",         help="Number of lanes to generate", default=3, type=int)
main_parser.add_argument("--total_coverage", metavar="total_coverage", help="Total coverage", default=1.0, type=float)
main_parser.add_argument("--mean_fragment_size", metavar="mean_fragment_size", help="Mean fragment size", default=350, type=int)
main_parser.add_argument("--sd_fragment_size", metavar="sd_fragment_size", help="Standard deviation of fragment size", default=50, type=int)
main_parser.add_argument("--cosmic_vcf", metavar="cosmic_vcf", help="COSMIC database VCF", default=[], required=True)
main_parser.add_argument("--normal_vcf", metavar="normal_vcf", help="Normal VCF from previous VarSim run", default=[], required=True)
main_parser.add_argument("--liftover_jar", metavar="liftover_jar", help="LiftOver jar", type=file, default=default_liftover, required=require_liftover)
main_parser.add_argument("--force_five_base_encoding", action="store_true", help="Force bases to be ACTGN")
main_parser.add_argument("--filter", action="store_true", help="Only use PASS variants")
main_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
main_parser.add_argument("--vcfstats_jar", metavar="JAR", help="VCFStats jar", type=file, default=default_vcfstats, required=require_vcfstats)
main_parser.add_argument("--varsim_py", metavar="PYTHON", help="VarSim python script", type=file, default=default_varsim, required=require_varsim)

pipeline_control_group = main_parser.add_argument_group("Pipeline control options. Disable parts of the pipeline.")
pipeline_control_group.add_argument("--disable_rand_vcf", action="store_true", help="Disable RandVCF2VCF somatic")
pipeline_control_group.add_argument("--disable_vcf2diploid", action="store_true", help="Disable vcf2diploid")
pipeline_control_group.add_argument("--disable_sim", action="store_true", help="Disable read simulation")

#RandVCF2VCF seed num_SNP num_INS num_DEL num_MNP num_COMPLEX percent_novel min_length_lim max_length_lim reference_file file.vcf
rand_vcf_group = main_parser.add_argument_group("RandVCF2VCF somatic options")
rand_vcf_group.add_argument("--som_num_snp", metavar="num_snp", help="Number of somatic SNPs", default=9000, type=int)
rand_vcf_group.add_argument("--som_num_ins", metavar="num_ins", help="Number of somatic insertions", default=1000, type=int);
rand_vcf_group.add_argument("--som_num_del", metavar="num_del", help="Number of somatic deletions", default=1000, type=int);
rand_vcf_group.add_argument("--som_num_mnp", metavar="num_mnp", help="Number of somatic MNPs", default=100, type=int)
rand_vcf_group.add_argument("--som_num_complex", metavar="num_complex", help="Number of somatic complex variants", default=100, type=int)
#rand_vcf_group.add_argument("--som_percent_novel", metavar="percent_novel", help="Percent novel", default=0, type=float)
rand_vcf_group.add_argument("--som_min_length_lim", metavar="min_length_lim", help="Min length lim", default=0, type=int)
rand_vcf_group.add_argument("--som_max_length_lim", metavar="max_length_lim", help="Max length lim", default=49, type=int)
#rand_vcf_group.add_argument("--som_vcf", metavar="in_vcf", help="Input somatic variant database VCF", type=file, required=False)
rand_vcf_group.add_argument("--som_prop_het", metavar="vc_prop_het", help="Proportion of somatic heterozygous variants", default=1.0, type=float)
rand_vcf_group.add_argument("--rand_vcf_jar", metavar="rand_vcf_jar", help="RandVCF2VCF jar", type=file, default=default_randvcf2vcf, required=require_randvcf2vcf)


dwgsim_group = main_parser.add_argument_group("DWGSIM options")
dwgsim_group.add_argument("--dwgsim_start_e", metavar="first_base_error_rate", help="Error rate on the first base", default=0.0001, type=float)
dwgsim_group.add_argument("--dwgsim_end_e", metavar="last_base_error_rate", help="Error rate on the last base", default=0.0015, type=float)
dwgsim_group.add_argument("--dwgsim_options", help="DWGSIM command-line options", default="", required=False)
dwgsim_group.add_argument("--dwgsim", metavar="dwgsim_executable", help="DWGSIM executable", type=file, required=False, default=None)


art_group = main_parser.add_argument_group("ART options")
art_group.add_argument("--profile_1", metavar="profile_file1", help="Profile for first end", default=None, type=file)
art_group.add_argument("--profile_2", metavar="profile_file2", help="Profile for second end", default=None, type=file)
art_group.add_argument("--art_options", help="ART command-line options", default="", required=False)
art_group.add_argument("--art", metavar="art_executable", help="ART executable", type=file, required=False, default=None)

args = main_parser.parse_args()

for d in [args.log_dir, args.out_dir, args.work_dir]:
  if not os.path.exists(d):
    os.makedirs(d)

FORMAT = '%(levelname)s %(asctime)-15s %(message)s'
logging.basicConfig(filename=os.path.join(args.log_dir, "varsim.log"), filemode="w", level=logging.DEBUG, format=FORMAT)
logger = logging.getLogger(__name__)

def run_shell_command(cmd, cmd_stdout, cmd_stderr, cmd_dir="."):
  subproc = subprocess.Popen(cmd, stdout=cmd_stdout, stderr=cmd_stderr, cwd=cmd_dir, shell=True, preexec_fn=os.setsid)
  retcode = subproc.wait()
  sys.exit(retcode)

def monitor_multiprocesses(processes, logger):
  for p in processes:
    p.join()
    if p.exitcode != 0:
      logger.error("Process with pid %d failed with exit code %d" % (p.pid, pid.exitcode))
    else:
      logger.info("Process with pid %d finished successfully" % (p.pid))

def monitor_processes(processes, logger):
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
              logger.write ("Could not kill the process\n")
      sys.exit(1)
    processes = processes_running
  return []

def check_executable(fpath):
  if not os.path.isfile(fpath):
    sys.stderr.write("ERROR: File %s does not exist\n" % (fpath))
    sys.exit(1)
  if not os.access(fpath, os.X_OK):
    sys.stderr.write("ERROR: File %s is not executable\n" % (fpath))
    sys.exit(1)

if not args.disable_sim:
  if args.simulator == "dwgsim":
    if args.dwgsim is None:
      sys.stderr.write("ERROR: Please specify the DWGSIM binary with --dwgsim option\n")
      sys.exit(1)
    check_executable(args.dwgsim.name)
  if args.simulator == "art":
    if args.art is None:
      sys.stderr.write("ERROR: Please specify the ART binary with --art option\n")
      sys.exit(1)
    check_executable(args.art.name)

processes = []

t_s = time.time()

vcf_files = [os.path.realpath(args.normal_vcf)]

if not args.disable_rand_vcf:
  rand_vcf_stdout = open(os.path.join(args.out_dir, "random.cosmic.vcf"), "w")
  rand_vcf_stderr = open(os.path.join(args.log_dir, "random.cosmic.err"), "w")
  vcf_files.insert(0,os.path.realpath(rand_vcf_stdout.name))

  rand_vcf_command = ["java", "-jar", os.path.realpath(args.rand_vcf_jar.name), "-seed", str(args.seed),
                 "-num_snp", str(args.som_num_snp), 
                 "-num_ins", str(args.som_num_ins), 
                 "-num_del", str(args.som_num_del),
                 "-num_mnp", str(args.som_num_mnp), 
                 "-num_complex", str(args.som_num_complex), 
                 #"-novel", str(args.som_percent_novel),
                 "-min_len", str(args.som_min_length_lim), 
                 "-max_len", str(args.som_max_length_lim), 
                 "-ref", os.path.realpath(args.reference.name),
                 "-prop_het", str(args.som_prop_het), 
                 "-vcf", os.path.realpath(args.cosmic_vcf)]

  p_rand_vcf = subprocess.Popen(rand_vcf_command, stdout=rand_vcf_stdout, stderr=rand_vcf_stderr)
  logger.info("Executing command " + " ".join(rand_vcf_command) + " with pid " + str(p_rand_vcf.pid))
  processes.append(p_rand_vcf)

processes = monitor_processes(processes, logger)

processes = []
for in_vcf in vcf_files:
  out_prefix = os.path.basename(in_vcf)
  vcfstats_stdout = open(os.path.join(args.out_dir, "%s.stats" % (out_prefix)), "w")
  vcfstats_stderr = open(os.path.join(args.log_dir, "%s.vcfstats.err" % (out_prefix)), "w")
  vcfstats_command = ["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(args.vcfstats_jar.name), "-vcf", in_vcf]
  p_vcfstats = subprocess.Popen(vcfstats_command, stdout=vcfstats_stdout, stderr=vcfstats_stderr)
  logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
  processes.append(p_vcfstats)


# Run VarSim 
varsim_stdout = open(os.path.join(args.log_dir, "som_varsim.out"), "w")
varsim_stderr = open(os.path.join(args.log_dir, "som_varsim.log"), "w")

vcf_arg_list = ["--vcfs"] + vcf_files

# need to fix the store true ones
filter_arg_list = ["--filter"] if args.filter else []
disable_sim_arg_list = ["--disable_sim"] if args.disable_sim else []
force_five_base_encoding_arg_list = ["--force_five_base_encoding"] if args.force_five_base_encoding else []
keep_temp_arg_list = ["--keep_temp"] if args.keep_temp else []
dwgsim_arg_list = ["--dwgsim", args.dwgsim] if args.dwgsim is not None else []
profile_1_arg_list = ["--profile_1", args.profile_1] if args.profile_1 is not None else []
profile_2_arg_list = ["--profile_2", args.profile_2] if args.profile_2 is not None else []
art_arg_list = ["--art", args.art.name] if args.art is not None else []
varsim_command = ["python", os.path.realpath(args.varsim_py.name), 
				  "--out_dir", str(os.path.realpath(args.out_dir)),
                  "--work_dir", str(os.path.realpath(args.work_dir)), 
                  "--log_dir", str(os.path.realpath(args.log_dir)), 
                  "--reference", str(os.path.realpath(args.reference.name)),
                  "--seed", str(args.seed), 
                  "--sex", str(args.sex), 
                  "--id", str(args.id),
                  "--simulator", str(args.simulator), 
                  "--vcf2diploid_jar", str(os.path.realpath(args.vcf2diploid_jar.name)), 
                  "--read_length", str(args.read_length),
                  "--nlanes", str(args.nlanes), 
                  "--total_coverage", str(args.total_coverage), 
                  "--mean_fragment_size", str(args.mean_fragment_size), 
                  "--sd_fragment_size", str(args.sd_fragment_size), 
                  "--liftover_jar", str(os.path.realpath(args.liftover_jar.name)), 
                  "--vcfstats_jar", str(os.path.realpath(args.vcfstats_jar.name)), 
                  "--dwgsim_start_e", str(args.dwgsim_start_e), 
                  "--dwgsim_end_e", str(args.dwgsim_end_e), 
                  "--dwgsim_options", str(args.dwgsim_options), 
                  "--art_options", str(args.art_options), 
                  "--disable_rand_vcf",
                  "--disable_rand_dgv" ] + vcf_arg_list + filter_arg_list + disable_sim_arg_list + force_five_base_encoding_arg_list + keep_temp_arg_list + dwgsim_arg_list + profile_1_arg_list + profile_2_arg_list + art_arg_list
p_varsim = subprocess.Popen(varsim_command, stdout=varsim_stdout, stderr=varsim_stderr)
logger.info("Executing command " + " ".join(varsim_command) + " with pid " + str(p_varsim.pid))
processes.append(p_varsim)


processes = monitor_processes(processes, logger)

# grep out the cosmic variants
# This is a bit dodgy
# grep -v "COS" art_cosmic/out/sv.truth.vcf > out/sv_norm.vcf &
# grep  "COS" art_cosmic/out/sv.truth.vcf > out/sv_cosmic.vcf &

grep_norm_stdout = open(os.path.join(args.out_dir, str(args.id) + "_norm.vcf"), "w")
grep_norm_stderr = open(os.path.join(args.log_dir, str(args.id) + "_norm.err"), "w")

grep_norm_command = ["grep", "-v", "COS", os.path.realpath(str(args.out_dir) + "/" + str(args.id) +  ".truth.vcf")]

p_grep_norm = subprocess.Popen(grep_norm_command, stdout=grep_norm_stdout, stderr=grep_norm_stderr)
logger.info("Executing command " + " ".join(grep_norm_command) + " with pid " + str(p_grep_norm.pid))
processes.append(p_grep_norm)

grep_cos_stdout = open(os.path.join(args.out_dir, str(args.id) + "_somatic.vcf"), "w")
grep_cos_stderr = open(os.path.join(args.log_dir, str(args.id) + "_somatic.err"), "w")

grep_cos_command = ["grep", "COS", os.path.realpath(str(args.out_dir) + "/" + str(args.id) +  ".truth.vcf")]

p_grep_cos = subprocess.Popen(grep_cos_command, stdout=grep_cos_stdout, stderr=grep_cos_stderr)
logger.info("Executing command " + " ".join(grep_cos_command) + " with pid " + str(p_grep_cos.pid))
processes.append(p_grep_cos)


processes = monitor_processes(processes, logger)

vcfstats_stdout = open(os.path.join(args.out_dir, "%s.norm.vcf.stats" %  (args.id)), "w")
vcfstats_stderr = open(os.path.join(args.log_dir, "%s.norm.vcf.vcfstats.err" %  (args.id)), "w")
p_vcfstats = subprocess.Popen(["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(args.vcfstats_jar.name), "-vcf", os.path.join(args.out_dir, str(args.id) + "_norm.vcf")], stdout=vcfstats_stdout, stderr=vcfstats_stderr)
logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
processes.append(p_vcfstats)

vcfstats_stdout = open(os.path.join(args.out_dir, "%s.somatic.vcf.stats" %  (args.id)), "w")
vcfstats_stderr = open(os.path.join(args.log_dir, "%s.somatic.vcf.vcfstats.err" %  (args.id)), "w")
p_vcfstats = subprocess.Popen(["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(args.vcfstats_jar.name), "-vcf", os.path.join(args.out_dir, str(args.id) + "_somatic.vcf")], stdout=vcfstats_stdout, stderr=vcfstats_stderr)
logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
processes.append(p_vcfstats)

monitor_processes(processes, logger) 

logger.info("Done! (%g hours)" % ((time.time() - t_s)/3600.0))
