#!/usr/bin/python

import argparse
import os
import sys
import subprocess
import logging
import shutil
import time
import signal
from multiprocessing import Process


def get_contigs_list(reference):
    with open("%s.fai" % (reference)) as fai_file:
        contigs = [line.strip().split()[0] for line in fai_file.readlines()]
    return contigs


my_dir = os.path.dirname(os.path.realpath(__file__))

default_varsim_jar = os.path.join(my_dir, "VarSim.jar")

require_varsim_jar = not os.path.isfile(default_varsim_jar)

if not os.path.isfile(default_varsim_jar): require_varsim_jar = None

main_parser = argparse.ArgumentParser(description="VarSim: A high-fidelity simulation validation framework",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
main_parser.add_argument("--out_dir", metavar="DIR", help="Output directory for the simulated genome, reads and variants", required=False, default="out")
main_parser.add_argument("--work_dir", metavar="DIR", help="Work directory, currently not used", required=False, default="work")
main_parser.add_argument("--log_dir", metavar="DIR", help="Log files of all steps are kept here", required=False,
                         default="log")
main_parser.add_argument("--reference", metavar="FASTA", help="Reference genome that variants will be inserted into", required=True, type=file)
main_parser.add_argument("--seed", metavar="seed", help="Random number seed for reproducibility", type=int, default=0)
main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (MALE/FEMALE)", required=False, type=str,
                         choices=["MALE", "FEMALE"], default="MALE")
main_parser.add_argument("--id", metavar="ID", help="Sample ID to be put in output VCF file", required=True)
main_parser.add_argument("--simulator", metavar="SIMULATOR", help="Read simulator to use", required=False, type=str,
                         choices=["art", "dwgsim", "pbsim"], default="art")
main_parser.add_argument("--simulator_executable", metavar="PATH", help="Path to the executable of the read simulator chosen"
                         , required=True, type=file)
main_parser.add_argument("--varsim_jar", metavar="PATH", help="Path to VarSim.jar", type=file, default=default_varsim_jar,
                         required=require_varsim_jar)
main_parser.add_argument("--read_length", metavar="LENGTH", help="Length of read to simulate", default=100, type=int)
main_parser.add_argument("--nlanes", metavar="INTEGER", help="Number of lanes to generate, coverage will be divided evenly over the lanes. Simulation is parallized over lanes. Each lane will have its own pair of files", default=1, type=int)
main_parser.add_argument("--total_coverage", metavar="FLOAT", help="Total coverage to simulate", default=1.0, type=float)
main_parser.add_argument("--mean_fragment_size", metavar="INT", help="Mean fragment size to simulate", default=350,
                         type=int)
main_parser.add_argument("--sd_fragment_size", metavar="INT", help="Standard deviation of fragment size to simulate",
                         default=50, type=int)
main_parser.add_argument("--vcfs", metavar="VCF", help="Addtional list of VCFs to insert into genome, priority is lowest ... highest", nargs="+", default=[])
main_parser.add_argument("--force_five_base_encoding", action="store_true", help="Force output bases to be only ACTGN")
main_parser.add_argument("--filter", action="store_true", help="Only use PASS variants for simulation")
main_parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files after simulation")

pipeline_control_group = main_parser.add_argument_group("Pipeline control options. Disable parts of the pipeline.")
pipeline_control_group.add_argument("--disable_rand_vcf", action="store_true", help="Disable sampling from the provided small variant VCF")
pipeline_control_group.add_argument("--disable_rand_dgv", action="store_true", help="Disable sampline from the provided DGV file")
pipeline_control_group.add_argument("--disable_vcf2diploid", action="store_true", help="Disable diploid genome simulation")
pipeline_control_group.add_argument("--disable_sim", action="store_true", help="Disable read simulation")

# RandVCF2VCF seed num_SNP num_INS num_DEL num_MNP num_COMPLEX percent_novel min_length_lim max_length_lim reference_file file.vcf
rand_vcf_group = main_parser.add_argument_group("Small variant simulation options")
rand_vcf_group.add_argument("--vc_num_snp", metavar="INTEGER", help="Number of SNPs to sample from small variant VCF", default=0, type=int)
rand_vcf_group.add_argument("--vc_num_ins", metavar="INTEGER", help="Number of insertions to sample from small variant VCF", default=0, type=int)
rand_vcf_group.add_argument("--vc_num_del", metavar="INTEGER", help="Number of deletions to sample from small variant VCF", default=0, type=int)
rand_vcf_group.add_argument("--vc_num_mnp", metavar="INTEGER", help="Number of MNPs to sample from small variant VCF", default=0, type=int)
rand_vcf_group.add_argument("--vc_num_complex", metavar="INTEGER", help="Number of complex variants to sample from small variant VCF", default=0,
                            type=int)
rand_vcf_group.add_argument("--vc_percent_novel", metavar="FLOAT", help="Percent variants sampled from small variant VCF that will be moved to novel positions", default=0, type=float)
rand_vcf_group.add_argument("--vc_min_length_lim", metavar="INTEGER", help="Min length of small variant to accept [inclusive]", default=0, type=int)
rand_vcf_group.add_argument("--vc_max_length_lim", metavar="INTEGER", help="Max length of small variant to accept [inclusive]", default=99,
                            type=int)
rand_vcf_group.add_argument("--vc_in_vcf", metavar="VCF", help="Input small variant VCF, usually dbSNP", type=file, required=False)
rand_vcf_group.add_argument("--vc_prop_het", metavar="FLOAT", help="Proportion of heterozygous small variants", default=0.6,
                            type=float)

# RandDGV2VCF seed num_INS num_DEL num_DUP num_INV percent_novel min_length_lim max_length_lim reference_file insert_seq.txt dgv_file.txt
rand_dgv_group = main_parser.add_argument_group("Structural variant simulation options")
rand_dgv_group.add_argument("--sv_num_ins", metavar="INTEGER", help="Number of insertions to sample from DGV", default=20, type=int)
rand_dgv_group.add_argument("--sv_num_del", metavar="INTEGER", help="Number of deletions to sample from DGV", default=20, type=int)
rand_dgv_group.add_argument("--sv_num_dup", metavar="INTEGER", help="Number of duplications to sample from DGV", default=20, type=int)
rand_dgv_group.add_argument("--sv_num_inv", metavar="INTEGER", help="Number of inversions to sample from DGV", default=20, type=int)
rand_dgv_group.add_argument("--sv_percent_novel", metavar="FLOAT", help="Percent variants sampled from DGV that will be moved to novel positions", default=0, type=float)
rand_dgv_group.add_argument("--sv_min_length_lim", metavar="min_length_lim", help="Min length of structural variant to accept [inclusive]", default=100,
                            type=int)
rand_dgv_group.add_argument("--sv_max_length_lim", metavar="max_length_lim", help="Max length of structural variant to accept [inclusive]", default=1000000,
                            type=int)
rand_dgv_group.add_argument("--sv_insert_seq", metavar="FILE", help="Path to file containing concatenation of real insertion sequences", type=file, required=False)
rand_dgv_group.add_argument("--sv_dgv", metavar="DGV_FILE", help="DGV file containing structural variants", type=file, required=False)

dwgsim_group = main_parser.add_argument_group("DWGSIM options")
dwgsim_group.add_argument("--dwgsim_start_e", metavar="first_base_error_rate", help="Error rate on the first base",
                          default=0.0001, type=float)
dwgsim_group.add_argument("--dwgsim_end_e", metavar="last_base_error_rate", help="Error rate on the last base",
                          default=0.0015, type=float)
dwgsim_group.add_argument("--dwgsim_options", help="DWGSIM command-line options", default="", required=False)

art_group = main_parser.add_argument_group("ART options")
art_group.add_argument("--profile_1", metavar="profile_file1", help="ART error profile for first end", default=None, type=file)
art_group.add_argument("--profile_2", metavar="profile_file2", help="ART error profile for second end", default=None, type=file)
art_group.add_argument("--art_options", help="ART command-line options", default="", required=False)

pbsim_group = main_parser.add_argument_group("PBSIM options")
pbsim_group.add_argument("--model_qc", metavar="model_qc", help="PBSIM QC model", default=None, type=str)
pbsim_group.add_argument("--accuracy_mean", metavar="accuracy_mean", help="PBSIM mean accuracy", default=None, type=float)

args = main_parser.parse_args()

# make the directories we need
for d in [args.log_dir, args.out_dir, args.work_dir]:
    if not os.path.exists(d):
        os.makedirs(d)

# Setup logging
FORMAT = '%(levelname)s %(asctime)-15s %(message)s'
logging.basicConfig(filename=os.path.join(args.log_dir, "varsim.log"), filemode="w", level=logging.DEBUG, format=FORMAT)
logger = logging.getLogger(__name__)


# ####### Some functions here ##########
def run_shell_command(cmd, cmd_stdout, cmd_stderr, cmd_dir="."):
    subproc = subprocess.Popen(cmd, stdout=cmd_stdout, stderr=cmd_stderr, cwd=cmd_dir, shell=True, preexec_fn=os.setsid)
    retcode = subproc.wait()
    sys.exit(retcode)


def monitor_multiprocesses(processes, logger):
    for p in processes:
        p.join()
        if p.exitcode != 0:
            logger.error("Process with pid %d failed with exit code %d" % (p.pid, p.exitcode))  # Marghoob: pid?
        else:
            logger.info("Process with pid %d finished successfully" % p.pid)


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
                    logger.error("Process %s failed. Will kill the remaining processes." % p.pid)
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
                            logger.write("Could not kill the process\n")
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

# ####### END Some functions here ##########



# Make sure we can actually execute the executable
if not args.disable_sim:
    check_executable(args.simulator_executable.name)

processes = []

t_s = time.time()

args.vcfs = [os.path.realpath(vcf) for vcf in args.vcfs]

if not args.disable_rand_vcf:
    rand_vcf_stdout = open(os.path.join(args.out_dir, "random.vc.vcf"), "w")
    rand_vcf_stderr = open(os.path.join(args.log_dir, "RandVCF2VCF.err"), "w")
    args.vcfs.append(os.path.realpath(rand_vcf_stdout.name))

    rand_vcf_command = ["java", "-jar", os.path.realpath(args.varsim_jar.name), "randvcf2vcf", "-seed", str(args.seed),
                        "-t", args.sex,
                        "-num_snp", str(args.vc_num_snp), "-num_ins", str(args.vc_num_ins), "-num_del",
                        str(args.vc_num_del),
                        "-num_mnp", str(args.vc_num_mnp), "-num_complex", str(args.vc_num_complex), "-novel",
                        str(args.vc_percent_novel),
                        "-min_len", str(args.vc_min_length_lim), "-max_len", str(args.vc_max_length_lim), "-ref",
                        os.path.realpath(args.reference.name),
                        "-prop_het", str(args.vc_prop_het), "-vcf", os.path.realpath(args.vc_in_vcf.name)]

    p_rand_vcf = subprocess.Popen(rand_vcf_command, stdout=rand_vcf_stdout, stderr=rand_vcf_stderr)
    logger.info("Executing command " + " ".join(rand_vcf_command) + " with pid " + str(p_rand_vcf.pid))
    processes.append(p_rand_vcf)

if not args.disable_rand_dgv:
    rand_dgv_stdout = open(os.path.join(args.out_dir, "random.sv.vcf"), "w")
    rand_dgv_stderr = open(os.path.join(args.log_dir, "RandDGV2VCF.err"), "w")
    args.vcfs.append(os.path.realpath(rand_dgv_stdout.name))

    rand_dgv_command = ["java", "-Xms10g", "-Xmx10g", "-jar", os.path.realpath(args.varsim_jar.name), "randdgv2vcf",
                        "-t", args.sex,
                        "-seed", str(args.seed),
                        "-num_ins", str(args.sv_num_ins), "-num_del", str(args.sv_num_del), "-num_dup",
                        str(args.sv_num_dup), "-num_inv", str(args.sv_num_inv),
                        "-novel", str(args.sv_percent_novel), "-min_len", str(args.sv_min_length_lim), "-max_len",
                        str(args.sv_max_length_lim), "-ref", os.path.realpath(args.reference.name),
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
    vcfstats_command = ["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(args.varsim_jar.name), "vcfstats", "-vcf",
                        in_vcf]
    p_vcfstats = subprocess.Popen(vcfstats_command, stdout=vcfstats_stdout, stderr=vcfstats_stderr)
    logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
    processes.append(p_vcfstats)

if not args.disable_vcf2diploid:
    args.vcfs.reverse()
    vcf2diploid_stdout = open(os.path.join(args.out_dir, "vcf2diploid.out"), "w")
    vcf2diploid_stderr = open(os.path.join(args.log_dir, "vcf2diploid.err"), "w")
    vcf_arg_list = sum([["-vcf", v] for v in args.vcfs], [])
    filter_arg_list = ["-pass"] if args.filter else []
    vcf2diploid_command = ["java", "-jar", os.path.realpath(args.varsim_jar.name), "vcf2diploid",
                           "-t", args.sex,
                           "-id", args.id,
                           "-chr", os.path.realpath(args.reference.name)] + filter_arg_list + vcf_arg_list

    p_vcf2diploid = subprocess.Popen(vcf2diploid_command, stdout=vcf2diploid_stdout, stderr=vcf2diploid_stderr,
                                     cwd=args.out_dir)
    logger.info("Executing command " + " ".join(vcf2diploid_command) + " with pid " + str(p_vcf2diploid.pid))
    processes.append(p_vcf2diploid)

    processes = monitor_processes(processes, logger)

    # Now concatenate the .fa from vcf2diploid
    contigs = get_contigs_list(args.reference.name)
    with open(merged_reference, "w") as merged_fa:
        for contig in contigs:
            for strand in ["maternal", "paternal"]:
                fasta = os.path.join(args.out_dir, "%s_%s_%s.fa" % (contig, args.id, strand))
                if not os.path.isfile(fasta):
                    continue
                with open(fasta) as chromosome_fasta:
                    shutil.copyfileobj(chromosome_fasta, merged_fa)

    # contatenate the vcfs
    with open(merged_truth_vcf, "w") as merged_vcf:
        first_file = True
        for contig in contigs:
            chr_vcf = os.path.join(args.out_dir, "%s_%s.vcf" % (contig, args.id))
            if not os.path.isfile(chr_vcf):
                continue
            with open(chr_vcf) as chr_vcf_file:
                for line in chr_vcf_file:
                    line = line.strip()
                    if len(line) == 0:
                        continue
                    if not first_file:
                        if line[0] == "#":
                            continue
                    merged_vcf.write(line)
                    merged_vcf.write('\n')

                first_file = False


    # Merge the chain files
    with open(merged_chain, "w") as merged_chain_file:
        for strand in ["maternal", "paternal"]:
            with open(os.path.join(args.out_dir, "%s.chain" % (strand))) as strand_chain:
                shutil.copyfileobj(strand_chain, merged_chain_file)

    vcfstats_stdout = open(os.path.join(args.out_dir, "%s.truth.vcf.stats" % (args.id)), "w")
    vcfstats_stderr = open(os.path.join(args.log_dir, "%s.truth.vcf.vcfstats.err" % (args.id)), "w")
    p_vcfstats = subprocess.Popen(
        ["java", "-Xmx1g", "-Xms1g", "-jar", os.path.realpath(args.varsim_jar.name), "vcfstats", "-vcf",
         merged_truth_vcf], stdout=vcfstats_stdout, stderr=vcfstats_stderr)
    logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
    monitor_processes([p_vcfstats], logger)

if processes:
    processes = monitor_processes(processes, logger)

merged_map = os.path.join(args.out_dir, "%s.map" % (args.id))

# Now generate the reads using art/pbsim/dwgsim
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
                fifo_src_dst.append(
                    ("simulated.lane%d.read%d.fastq" % (i, end),
                     "simulated.lane%d.read%d.fq.gz" % (i, end)))
    elif args.simulator == "art":
        for i in xrange(args.nlanes):
            for end in [1, 2]:
                for suffix in ["fq", "aln"]:
                    fifo_src_dst.append(("simulated.lane%d.read%d.%s" % (i, end, suffix),
                                         "simulated.lane%d.read%d.%s.gz" % (i, end, suffix)))
    elif args.simulator == "pbsim":
        for i in xrange(args.nlanes):
            for end in [1, 2]: # the '2' read files are empty, and for compatibility only
                for suffix in ["fq", "maf"]:
                    fifo_src_dst.append(("simulated.lane%d.read%d.%s" % (i, end, suffix),
                                         "simulated.lane%d.read%d.%s.gz" % (i, end, suffix)))
    else:
        raise NotImplementedError("simulation method "+args.simulator+" not implemented");

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
            dwgsim_command = [os.path.realpath(args.simulator_executable.name), "-e",
                              "%s,%s" % (args.dwgsim_start_e, args.dwgsim_end_e),
                              "-E", "%s,%s" % (args.dwgsim_start_e, args.dwgsim_end_e), args.dwgsim_options,
                              "-d", str(args.mean_fragment_size), "-s", str(args.sd_fragment_size), "-C",
                              str(coverage_per_lane), "-1", str(args.read_length), "-2", str(args.read_length),
                              "-z", str(i), merged_reference, os.path.join(args.out_dir, "simulated.lane%d" % (i))]
            dwgsim_command = " ".join(dwgsim_command)

            dwgsim_stdout = open(os.path.join(args.log_dir, "dwgsim.lane%d.out" % (i)), "w")
            dwgsim_stderr = open(os.path.join(args.log_dir, "dwgsim.lane%d.err" % (i)), "w")
            dwgsim_p = Process(target=run_shell_command, args=(dwgsim_command, dwgsim_stdout, dwgsim_stderr))
            dwgsim_p.start()
            processes.append(dwgsim_p)
            logger.info("Executing command " + dwgsim_command + " with pid " + str(dwgsim_p.pid))
    elif args.simulator == "art":
        profile_opts = []
        if args.profile_1 is not None and args.profile_2 is not None:
            profile_opts = ["-1", args.profile_1.name, "-2", args.profile_2.name]

        for i in xrange(args.nlanes):
            art_command = [args.simulator_executable.name] + profile_opts + ["-i", merged_reference, "-p",
                                                                             "-l", str(args.read_length), "-f",
                                                                             str(coverage_per_lane),
                                                                             "-m", str(args.mean_fragment_size),
                                                                             "-s", str(args.sd_fragment_size), "-rs",
                                                                             str(i),
                                                                             args.art_options, "-o",
                                                                             os.path.join(args.out_dir,
                                                                                          "simulated.lane%d.read" % (
                                                                                              i))]
            art_command = " ".join(art_command)
            art_stdout = open(os.path.join(args.log_dir, "art.lane%d.out" % (i)), "w")
            art_stderr = open(os.path.join(args.log_dir, "art.lane%d.err" % (i)), "w")
            art_p = Process(target=run_shell_command, args=(art_command, art_stdout, art_stderr))
            art_p.start()
            processes.append(art_p)
            logger.info("Executing command " + art_command + " with pid " + str(art_p.pid))
    elif args.simulator == "pbsim":
        nRef = 0;
        with open(merged_reference,'r') as fa:
            nRef = sum( 1 for line in fa if len(line)>0 and line[0] == '>' )
        assert nRef > 0 and nRef < 10000

        for i in xrange(args.nlanes):
            tmp_prefix = os.path.join(args.out_dir,"simulated.lane%d"%(i));
            tmp_fastq_list = " ".join( "%s_%s.fastq"%(tmp_prefix,"0"*(4-len(str(idx)))+str(idx)) for idx in range(1,nRef+1) )
            tmp_maf_list = " ".join( "%s_%s.maf"%(tmp_prefix,"0"*(4-len(str(idx)))+str(idx)) for idx in range(1,nRef+1) )
            tmp_ref_list = " ".join( "%s_%s.ref"%(tmp_prefix,"0"*(4-len(str(idx)))+str(idx)) for idx in range(1,nRef+1) )
            pbsim_command = [os.path.realpath(args.simulator_executable.name),
                             "--data-type", "CLR",
                             "--difference-ratio", "1:12:2",
                             "--depth", str(coverage_per_lane),
                             "--model_qc", args.model_qc,
                             "--seed", str(2089*(i+1)),
                             "--length-max", "50000",
                             "--length-mean", str(args.mean_fragment_size),
                             "--length-sd", str(args.sd_fragment_size),
                             "--accuracy-mean", str(args.accuracy_mean),
                             merged_reference,
                             "--prefix", tmp_prefix,
                             "&& ( cat", tmp_fastq_list, "> %s.read1.fq"%(tmp_prefix), ")", #this cat is i/o bound, need optimization of piping
                             "&& ( cat", tmp_maf_list, "> %s.read1.maf"%(tmp_prefix), ")" #this cat is i/o bound, need optimization of piping
                             "&& ( head -q -n1 ", tmp_ref_list, "> %s.ref"%(tmp_prefix), ")" #make reference header list
                             "&& ( echo > %s.read2.fq"%(tmp_prefix), ")", #dummy file
                             "&& ( echo > %s.read2.maf"%(tmp_prefix), ")" #dummy file
                             ]

            pbsim_command = " ".join(pbsim_command)

            pbsim_stdout = open(os.path.join(args.log_dir, "pbsim.lane%d.out" % (i)), "w")
            pbsim_stderr = open(os.path.join(args.log_dir, "pbsim.lane%d.err" % (i)), "w")
            pbsim_p = Process(target=run_shell_command, args=(pbsim_command, pbsim_stdout, pbsim_stderr))
            pbsim_p.start()
            processes.append(pbsim_p)
            logger.info("Executing command " + pbsim_command + " with pid " + str(pbsim_p.pid))
    else:
        raise NotImplementedError("simulation method "+args.simulator+" not implemented");



    monitor_multiprocesses(processes, logger)
    processes = []

    logger.info("Read generation took %g seconds" % (time.time() - sim_ts))

    sim_t_liftover = time.time()

    # Now start lifting over the gzipped files
    for i in xrange(args.nlanes):
        liftover_stdout = open(os.path.join(args.log_dir, "lane%d.out" % (i)), "w")
        liftover_stderr = open(os.path.join(args.log_dir, "liftover%d.log" % (i)), "w")
        fastq_liftover_command = "java -server -Xms4g -Xmx4g -jar %s fastq_liftover -map %s -id %d " \
                                 "-fastq <(gunzip -c %s/simulated.lane%d.read1.fq.gz) " \
                                 "-fastq <(gunzip -c %s/simulated.lane%d.read2.fq.gz) " \
                                 "-out >(gzip -1 > %s/lane%d.read1.fq.gz) " \
                                 "-out >(gzip -1 > %s/lane%d.read2.fq.gz)" % (
                                     os.path.realpath(args.varsim_jar.name), merged_map, i, args.out_dir, i,
                                     args.out_dir, i, args.out_dir, i,
                                     args.out_dir, i)
        if args.force_five_base_encoding:
            fastq_liftover_command += " -force_five_base_encoding "
        if args.simulator == "art":
            fastq_liftover_command += " -type art " \
                                      "-aln <(gunzip -c %s/simulated.lane%d.read1.aln.gz) " \
                                      "-aln <(gunzip -c %s/simulated.lane%d.read2.aln.gz)" % (
                                          args.out_dir, i, args.out_dir, i)
        elif args.simulator == "pbsim":
            fastq_liftover_command += " -type pbsim " \
                                      "-maf <(gunzip -c %s/simulated.lane%d.read1.maf.gz) " \
                                      "-ref %s/simulated.lane%d.ref "% (args.out_dir, i, args.out_dir, i)
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
    logger.info("Took %g seconds, %ld Mbytes written, %g MB/s" % (
        sim_te - sim_ts, bytes_written / 1024.0 / 1024.0, bytes_written / 1024.0 / 1024.0 / (sim_te - sim_ts)))

    for fifo in fifos:
        os.remove(fifo)

if not args.keep_temp:
    for f in tmp_files:
        os.remove(f)
logger.info("Done! (%g hours)" % ((time.time() - t_s) / 3600.0))
