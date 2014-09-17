#!/usr/bin/python

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
default_randdgv2vcf = os.path.join(my_dir, "target/build/randdgv2vcf.jar")
default_liftover    = os.path.join(my_dir, "target/build/fastq_liftover.jar")
default_vcfstats    = os.path.join(my_dir, "target/build/vcfstats.jar")

require_randvcf2vcf = not os.path.isfile(default_randvcf2vcf)
require_randdgv2vcf = not os.path.isfile(default_randdgv2vcf)
require_vcf2diploid = not os.path.isfile(default_vcf2diploid)
require_liftover    = not os.path.isfile(default_liftover)
require_vcfstats    = not os.path.isfile(default_vcfstats)

if not os.path.isfile(default_randvcf2vcf): default_randvcf2vcf = None
if not os.path.isfile(default_randdgv2vcf): default_randdgv2vcf = None
if not os.path.isfile(default_vcf2diploid): default_vcf2diploid = None
if not os.path.isfile(default_liftover):    default_liftover    = None
if not os.path.isfile(default_vcfstats):    default_vcfstats    = None

main_parser = argparse.ArgumentParser(description="VarSim: An accurate reads simulator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
main_parser.add_argument("--out_dir", metavar="Out directory", help="Output directory", required=False, default="out")
main_parser.add_argument("--work_dir", metavar="Work directory", help="Work directory", required=False, default="work")
main_parser.add_argument("--log_dir", metavar="Log directory", help="Directory to log to", required=False, default="log")
main_parser.add_argument("--reference", metavar="Reference", help="Reference file", required=True, type=file)
main_parser.add_argument("--seed", metavar="seed", help="Random number seed", type=int, default=0)
main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (male/female)", required=False, type=str, choices=["male", "female"], default="male")
main_parser.add_argument("--id", metavar="id", help="Sample ID", required=True)
main_parser.add_argument("--simulator", metavar="simulator", help="Read simulator", required=False, type=str, choices=["art", "dwgsim"], default="art")
main_parser.add_argument("--vcf2diploid_jar", metavar="vcf2diploid_jar", help="vcf2diploid jar", type=file, default=default_vcf2diploid, required=require_vcf2diploid)
main_parser.add_argument("--read_length",    metavar="read_length",    help="Length of reads", default=100, type=int)
main_parser.add_argument("--nlanes",         metavar="nlanes",         help="Number of lanes to generate", default=1, type=int)
main_parser.add_argument("--total_coverage", metavar="total_coverage", help="Total coverage", default=1.0, type=float)
main_parser.add_argument("--mean_fragment_size", metavar="mean_fragment_size", help="Mean fragment size", default=350, type=int)
main_parser.add_argument("--sd_fragment_size", metavar="sd_fragment_size", help="Standard deviation of fragment size", default=50, type=int)
main_parser.add_argument("--vcfs", metavar="vcfs", help="VCF list", nargs="+", default=[])
main_parser.add_argument("--liftover_jar", metavar="liftover_jar", help="LiftOver jar", type=file, default=default_liftover, required=require_liftover)
main_parser.add_argument("--force_five_base_encoding", action="store_true", help="Force bases to be ACTGN")
main_parser.add_argument("--filter", action="store_true", help="Only use PASS variants")
main_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
main_parser.add_argument("--vcfstats_jar", metavar="JAR", help="VCFStats jar", type=file, default=default_vcfstats, required=require_vcfstats)

pipeline_control_group = main_parser.add_argument_group("Pipeline control options. Disable parts of the pipeline.")
pipeline_control_group.add_argument("--disable_rand_vcf", action="store_true", help="Disable RandVCF2VCF")
pipeline_control_group.add_argument("--disable_rand_dgv", action="store_true", help="Disable RandDGV2VCF")
pipeline_control_group.add_argument("--disable_vcf2diploid", action="store_true", help="Disable vcf2diploid")
pipeline_control_group.add_argument("--disable_sim", action="store_true", help="Disable read simulation")

#RandVCF2VCF seed num_SNP num_INS num_DEL num_MNP num_COMPLEX percent_novel min_length_lim max_length_lim reference_file file.vcf
rand_vcf_group = main_parser.add_argument_group("RandVCF2VCF options")
rand_vcf_group.add_argument("--vc_num_snp", metavar="num_snp", help="Number of SNPs", default=0, type=int)
rand_vcf_group.add_argument("--vc_num_ins", metavar="num_ins", help="Number of insertions", default=0, type=int);
rand_vcf_group.add_argument("--vc_num_del", metavar="num_del", help="Number of deletions", default=0, type=int);
rand_vcf_group.add_argument("--vc_num_mnp", metavar="num_mnp", help="Number of MNPs", default=0, type=int)
rand_vcf_group.add_argument("--vc_num_complex", metavar="num_complex", help="Number of complex variants", default=0, type=int)
rand_vcf_group.add_argument("--vc_percent_novel", metavar="percent_novel", help="Percent novel", default=0, type=float)
rand_vcf_group.add_argument("--vc_min_length_lim", metavar="min_length_lim", help="Min length lim", default=0, type=int)
rand_vcf_group.add_argument("--vc_max_length_lim", metavar="max_length_lim", help="Max length lim", default=50, type=int)
rand_vcf_group.add_argument("--vc_in_vcf", metavar="in_vcf", help="Input VCF", type=file, required=True)
rand_vcf_group.add_argument("--vc_prop_het", metavar="vc_prop_het", help="Proportion of heterozygous vars", default=0.6, type=float)
rand_vcf_group.add_argument("--rand_vcf_jar", metavar="rand_vcf_jar", help="RandVCF2VCF jar", type=file, default=default_randvcf2vcf, required=require_randvcf2vcf)

#RandDGV2VCF seed num_INS num_DEL num_DUP num_INV percent_novel min_length_lim max_length_lim reference_file insert_seq.txt dgv_file.txt
rand_dgv_group = main_parser.add_argument_group("RandDGV2VCF options")
rand_dgv_group.add_argument("--sv_num_ins", metavar="num_ins", help="Number of insertions", default=20, type=int);
rand_dgv_group.add_argument("--sv_num_del", metavar="num_del", help="Number of deletions", default=20, type=int);
rand_dgv_group.add_argument("--sv_num_dup", metavar="num_dup", help="Number of duplications", default=20, type=int);
rand_dgv_group.add_argument("--sv_num_inv", metavar="num_inv", help="Number of inversions", default=20, type=int);
rand_dgv_group.add_argument("--sv_percent_novel", metavar="percent_novel", help="Percent novel", default=0, type=float)
rand_dgv_group.add_argument("--sv_min_length_lim", metavar="min_length_lim", help="Min length lim", default=50, type=int)
rand_dgv_group.add_argument("--sv_max_length_lim", metavar="max_length_lim", help="Max length lim", default=1000000, type=int)
rand_dgv_group.add_argument("--sv_insert_seq", metavar="insert_seq", help="Insert seq", type=file, required=True)
rand_dgv_group.add_argument("--sv_dgv", metavar="dgv", help="DGV file", type=file, required=True)
rand_dgv_group.add_argument("--rand_dgv_jar", metavar="rand_dgv_jar", help="RandDGV2VCF jar", type=file, default=default_randdgv2vcf, required=require_randdgv2vcf)

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

args.vcfs = [os.path.realpath(vcf) for vcf in args.vcfs]

if not args.disable_rand_vcf:
  rand_vcf_stdout = open(os.path.join(args.out_dir, "random.vc.vcf"), "w")
  rand_vcf_stderr = open(os.path.join(args.log_dir, "RandVCF2VCF.err"), "w")
  args.vcfs.append(os.path.realpath(rand_vcf_stdout.name))

  rand_vcf_command = ["java", "-jar", os.path.realpath(args.rand_vcf_jar.name), "-seed", str(args.seed),
                 "-num_snp", str(args.vc_num_snp), "-num_ins", str(args.vc_num_ins), "-num_del", str(args.vc_num_del),
                 "-num_mnp", str(args.vc_num_mnp), "-num_complex", str(args.vc_num_complex), "-novel", str(args.vc_percent_novel),
                 "-min_len", str(args.vc_min_length_lim), "-max_len", str(args.vc_max_length_lim), "-ref", os.path.realpath(args.reference.name),
                 "-prop_het", str(args.vc_prop_het), "-vcf", os.path.realpath(args.vc_in_vcf.name)]

  p_rand_vcf = subprocess.Popen(rand_vcf_command, stdout=rand_vcf_stdout, stderr=rand_vcf_stderr)
  logger.info("Executing command " + " ".join(rand_vcf_command) + " with pid " + str(p_rand_vcf.pid))
  processes.append(p_rand_vcf)

if not args.disable_rand_dgv:
  rand_dgv_stdout = open(os.path.join(args.out_dir, "random.sv.vcf"), "w")
  rand_dgv_stderr = open(os.path.join(args.log_dir, "RandDGV2VCF.err"), "w")
  args.vcfs.append(os.path.realpath(rand_dgv_stdout.name))

  rand_dgv_command = ["java", "-Xms10g", "-Xmx10g", "-jar", os.path.realpath(args.rand_dgv_jar.name), "-seed", str(args.seed),
                    "-num_ins", str(args.sv_num_ins), "-num_del", str(args.sv_num_del), "-num_dup", str(args.sv_num_dup), "-num_inv", str(args.sv_num_inv),
                    "-novel", str(args.sv_percent_novel), "-min_len", str(args.sv_min_length_lim), "-max_len", str(args.sv_max_length_lim), "-ref", os.path.realpath(args.reference.name),
                    "-ins", os.path.realpath(args.sv_insert_seq.name), "-dgv", os.path.realpath(args.sv_dgv.name)]

  p_rand_dgv = subprocess.Popen(rand_dgv_command, stdout=rand_dgv_stdout, stderr=rand_dgv_stderr)
  logger.info("Executing command " + " ".join(rand_dgv_command) + " with pid " + str(p_rand_dgv.pid))
  processes.append(p_rand_dgv)

processes = monitor_processes(processes, logger)

merged_reference = os.path.join(args.out_dir, "%s.fa" % (args.id))
merged_truth_vcf = os.path.join(args.out_dir, "%s.truth.vcf" % (args.id))
merged_chain = os.path.join(args.out_dir, "merged.chain")

processes = []
for in_vcf in args.vcfs:
  out_prefix = os.path.basename(in_vcf)
  vcfstats_stdout = open(os.path.join(args.out_dir, "%s.stats" % (out_prefix)), "w")
  vcfstats_stderr = open(os.path.join(args.log_dir, "%s.vcfstats.err" % (out_prefix)), "w")
  vcfstats_command = ["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(args.vcfstats_jar.name), "-vcf", in_vcf]
  p_vcfstats = subprocess.Popen(vcfstats_command, stdout=vcfstats_stdout, stderr=vcfstats_stderr)
  logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
  processes.append(p_vcfstats)

if not args.disable_vcf2diploid:
  args.vcfs.reverse()
  vcf2diploid_stdout = open(os.path.join(args.out_dir, "vcf2diploid.out"), "w")
  vcf2diploid_stderr = open(os.path.join(args.log_dir, "vcf2diploid.err"), "w")
  vcf_arg_list = sum([["-vcf", v] for v in args.vcfs], [])
  filter_arg_list = ["-pass"] if args.filter else []
  vcf2diploid_command = ["java", "-jar", os.path.realpath(args.vcf2diploid_jar.name), "-t", args.sex, "-id", args.id, "-chr", os.path.realpath(args.reference.name)] + filter_arg_list + vcf_arg_list

  p_vcf2diploid = subprocess.Popen(vcf2diploid_command, stdout=vcf2diploid_stdout, stderr=vcf2diploid_stderr, cwd=args.out_dir)
  logger.info("Executing command " + " ".join(vcf2diploid_command) + " with pid " + str(p_vcf2diploid.pid))
  processes.append(p_vcf2diploid)

  processes = monitor_processes(processes, logger)

  # Now concatenate the .fa from vcf2diploid
  contigs = get_contigs_list(args.reference.name)
  with open(merged_reference, "w") as merged_fa:
    for contig in contigs:
      #continue
      for strand in ["maternal", "paternal"]:
        fasta = os.path.join(args.out_dir, "%s_%s_%s.fa" % (contig, args.id, strand))
        if not os.path.isfile(fasta): continue
        with open(fasta) as chromosome_fasta:
          shutil.copyfileobj(chromosome_fasta, merged_fa)
    #break
  # contatenate the vcfs
  with open(merged_truth_vcf, "w") as merged_vcf:
    for contig in contigs:
      chr_vcf = os.path.join(args.out_dir, "%s_%s.vcf" % (contig, args.id))
      if not os.path.isfile(chr_vcf): continue
      with open(chr_vcf) as chr_vcf_file:
        shutil.copyfileobj(chr_vcf_file, merged_vcf)

  # Merge the chain files
  with open(merged_chain, "w") as merged_chain_file:
    for strand in ["maternal", "paternal"]:
      with open(os.path.join(args.out_dir, "%s.chain" % (strand))) as strand_chain:
        shutil.copyfileobj(strand_chain, merged_chain_file)

  vcfstats_stdout = open(os.path.join(args.out_dir, "%s.truth.vcf.stats" %  (args.id)), "w")
  vcfstats_stderr = open(os.path.join(args.log_dir, "%s.truth.vcf.vcfstats.err" %  (args.id)), "w")
  p_vcfstats = subprocess.Popen(["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(args.vcfstats_jar.name), merged_truth_vcf], stdout=vcfstats_stdout, stderr=vcfstats_stderr)
  logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
  monitor_processes([p_vcfstats], logger) 

if processes:
  processes = monitor_processes(processes, logger)

merged_map = os.path.join(args.out_dir, "%s.map" % (args.id))

# Now generate the reads using dwgsim
tmp_files = []
if not args.disable_sim:
  fifos = []
  fastqs = []
  sim_ts = time.time()
  coverage_per_lane = args.total_coverage * 0.5 / args.nlanes
  processes = []

  fifo_src_dst = []
  if args.simulator == "dwgsim":
    for i in xrange(args.nlanes):
      for end in [1, 2]:
        fifo_src_dst.append(("simulated.lane%d.bwa.read%d.fastq" % (i, end), "simulated.lane%d.read%d.fq.gz" % (i, end)))
  if args.simulator == "art":
    for i in xrange(args.nlanes):
      for end in [1, 2]:
        for suffix in ["fq", "aln"]:
          fifo_src_dst.append(("simulated.lane%d.read%d.%s" % (i, end, suffix), "simulated.lane%d.read%d.%s.gz" % (i, end, suffix)))

  for fifo_name, dst in fifo_src_dst:
    fifos.append(os.path.join(args.out_dir, fifo_name))
    if os.path.exists(fifos[-1]): os.remove(fifos[-1])
    os.mkfifo(fifos[-1])

    gzip_stderr = open(os.path.join(args.log_dir, "gzip.%s" % (fifo_name)), "w")
    gzip_command = "cat %s | gzip -2 > %s" % (fifos[-1], os.path.join(args.out_dir, dst))
    gzip_p = Process(target=run_shell_command, args=(gzip_command, None, gzip_stderr))
    gzip_p.start()
    processes.append(gzip_p)
    logger.info("Executing command %s" % (gzip_command) + " with pid " + str(gzip_p.pid))
    tmp_files.append(os.path.join(args.out_dir, dst))

  if args.simulator == "dwgsim":
    for i in xrange(args.nlanes):
      dwgsim_command = [os.path.realpath(args.dwgsim.name), "-e", "%s,%s" % (args.dwgsim_start_e, args.dwgsim_end_e),
                      "-E", "%s,%s" % (args.dwgsim_start_e, args.dwgsim_end_e), args.dwgsim_options,
                      "-d", str(args.mean_fragment_size), "-s", str(args.sd_fragment_size), "-C", str(coverage_per_lane), "-1", str(args.read_length), "-2", str(args.read_length),
                      "-z", str(i), merged_reference, os.path.join(args.out_dir, "simulated.lane%d" % (i))]
      dwgsim_command = " ".join(dwgsim_command)

      dwgsim_stdout = open(os.path.join(args.log_dir, "dwgsim.lane%d.out" % (i)), "w")
      dwgsim_stderr = open(os.path.join(args.log_dir, "dwgsim.lane%d.err" % (i)), "w")
      dwgsim_p = Process(target=run_shell_command, args=(dwgsim_command, dwgsim_stdout, dwgsim_stderr))
      dwgsim_p.start()
      processes.append(dwgsim_p)
      logger.info("Executing command " + dwgsim_command + " with pid " + str(dwgsim_p.pid))

  if args.simulator == "art":
    profile_opts = []
    if args.profile_1 is not None and args.profile_2 is not None:
      profile_opts = ["-1", args.profile_1.name, "-2", args.profile_2.name]

    for i in xrange(args.nlanes):
      art_command = [args.art.name] + profile_opts + ["-i", merged_reference, "-p",
                     "-l", str(args.read_length), "-f", str(coverage_per_lane), "-m", str(args.mean_fragment_size),
                     "-s", str(args.sd_fragment_size), "-rs", str(i), args.art_options, "-o", os.path.join(args.out_dir, "simulated.lane%d.read" % (i))]
      art_command = " ".join(art_command)
      art_stdout = open(os.path.join(args.log_dir, "art.lane%d.out" % (i)), "w")
      art_stderr = open(os.path.join(args.log_dir, "art.lane%d.err" % (i)), "w")
      art_p = Process(target=run_shell_command, args=(art_command, art_stdout, art_stderr))
      art_p.start()
      processes.append(art_p)
      logger.info("Executing command " + art_command + " with pid " + str(art_p.pid))

  monitor_multiprocesses(processes, logger)
  processes = []

  logger.info("Read generation took %g seconds" % (time.time() - sim_ts))

  sim_t_liftover = time.time()

  # Now start lifting over the gzipped files
  for i in xrange(args.nlanes):
    liftover_stdout = open(os.path.join(args.log_dir, "lane%d.out" % (i)), "w")
    liftover_stderr = open(os.path.join(args.log_dir, "liftover%d.log" % (i)), "w")
    fastq_liftover_command = "java -server -Xms4g -Xmx4g -jar %s -map %s -id %d -fastq <(gunzip -c %s/simulated.lane%d.read1.fq.gz) -fastq <(gunzip -c %s/simulated.lane%d.read2.fq.gz) -out >(gzip -1 > %s/lane%d.read1.fq.gz) -out >(gzip -1 > %s/lane%d.read2.fq.gz)" % (os.path.realpath(args.liftover_jar.name), merged_map, i, args.out_dir, i, args.out_dir, i, args.out_dir, i, args.out_dir, i)
    if args.force_five_base_encoding: fastq_liftover_command += " -force_five_base_encoding "
    if args.simulator == "art": fastq_liftover_command += " -type art -aln <(gunzip -c %s/simulated.lane%d.read1.aln.gz) -aln <(gunzip -c %s/simulated.lane%d.read2.aln.gz)" % (args.out_dir, i, args.out_dir, i)
    fastq_liftover_command = "bash -c \"%s\"" % (fastq_liftover_command)
    liftover_p = Process(target=run_shell_command, args=(fastq_liftover_command, liftover_stdout, liftover_stderr))
    liftover_p.start()
    processes.append(liftover_p)
    fastqs.append(os.path.join(args.out_dir, "lane%d.read%d.fq.gz" % (i, end)))
    logger.info("Executing command " + fastq_liftover_command + " with pid " + str(liftover_p.pid))

  monitor_multiprocesses(processes, logger)

  logger.info("Liftover took %g seconds" % (time.time() - sim_t_liftover))

  sim_te = max(sim_ts + 1, time.time())
  bytes_written = sum([os.path.getsize(fastq) for fastq in fastqs])
  logger.info("Took %g seconds, %ld Mbytes written, %g MB/s" % (sim_te - sim_ts, bytes_written/1024.0/1024.0, bytes_written/1024.0/1024.0/(sim_te - sim_ts)))

  for fifo in fifos:
    os.remove(fifo)

if not args.keep_temp:
  for f in tmp_files:
    os.remove(f)
logger.info("Done! (%g hours)" % ((time.time() - t_s)/3600.0))
