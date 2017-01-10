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
from multiprocessing import Process
from liftover_restricted_vcf_map import lift_vcfs, lift_maps

MY_DIR = os.path.dirname(os.path.realpath(__file__))
VARSIMJAR = os.path.realpath(os.path.join(MY_DIR, "VarSim.jar"))
DEFAULT_VARSIMJAR = os.path.join(MY_DIR, "VarSim.jar")
REQUIRE_VARSIMJAR = not os.path.isfile(DEFAULT_VARSIMJAR)
if REQUIRE_VARSIMJAR: DEFAULT_VARSIMJAR = None

def get_loglevel(string):
    if string == "info":
        return logging.INFO
    if string == "warn":
        return logging.WARN
    if string == "debug":
        return logging.DEBUG
    return logging.INFO


def convertCN(filenames, operation):
    """
    convert '2/1'-like copy number to a single number(e.g. 2)
    0 will be considered same as 1
    by default the max number will be kept
    the change is in place
    """
    if operation != "two2one" and operation != "one2two":
        raise ValueError("Only two2one or one2two allowed")
    two2one = operation == "two2one"
    delimiter = re.compile('[/|]')
    for name in filenames:
        with open(name, 'r') as file_fd:
            output = tempfile.NamedTemporaryFile(mode = 'r+w', delete = False)
            for l in file_fd:
                l = l.rstrip()
                fields = l.split("\t")
                if l.startswith("#") or 'CN' not in fields[8]:
                    if l.startswith('##FORMAT=<ID=CN'):
                        if two2one:
                            l = l.replace("Type=String","Type=Integer")
                        else:
                            l = l.replace("Type=Integer", "Type=String")
                    output.write(l + "\n")
                else:
                    info = fields[8].split(':')
                    cnIndex = info.index('CN')
                    gtIndex = info.index('GT')
                    #change CN field in all samples
                    for sampleIndex in range(9,len(fields)):
                        sampleInfo = fields[sampleIndex].split(':')
                        if two2one:
                            cn = delimiter.split(sampleInfo[cnIndex])
			    #here cn is list of strings
			    sampleInfo[cnIndex] = str(max(map(int, cn)))
                        elif len(delimiter.split(sampleInfo[cnIndex])) == 1:
                            #only split when there is only one number
                            gt = delimiter.split(sampleInfo[gtIndex])
                            cn = sampleInfo[cnIndex]
                            for i in range(len(gt)):
                                gt[i] = '1' if gt[i] == '0' else cn
                            if sampleInfo[gtIndex].find('/') >= 0:
                                sampleInfo[cnIndex] = '/'.join(gt)
                            else:
                                sampleInfo[cnIndex] = '|'.join(gt)
                        fields[sampleIndex] = ":".join(sampleInfo)
                    output.write("\t".join(fields) + "\n")
            output.close()
            shutil.copyfile(output.name, name)
            os.remove(output.name)
    return


def get_contigs_list(reference):
    with open("%s.fai" % (reference)) as fai_file:
        contigs = [line.strip().split()[0] for line in fai_file.readlines()]
    return contigs

# Check java version to make sure it is Java 8
def check_java():
    jv = filter(lambda x: x.startswith("java version"), subprocess.check_output("java -version", stderr=subprocess.STDOUT, shell=True).split("\n"))[0].split()[2].replace("\"", "")
    if LooseVersion(jv) < LooseVersion("1.8"):
        logger.error("VarSim requires Java 1.8 to be on the path.")
        raise EnvironmentError("VarSim requires Java 1.8 to be on the path")


def run_shell_command(cmd, cmd_stdout, cmd_stderr, cmd_dir="."):
    subproc = subprocess.Popen(cmd, stdout=cmd_stdout, stderr=cmd_stderr, cwd=cmd_dir, shell=True, preexec_fn=os.setsid)
    retcode = subproc.wait()
    sys.exit(retcode)


def makedirs(dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)


def monitor_multiprocesses(processes, logger):
    for p in processes:
        p.join()
        if p.exitcode != 0:
            logger.error("Process with pid %d failed with exit code %d" % (p.pid, p.exitcode))  # Marghoob: pid?
        else:
            logger.info("Process with pid %d finished successfully" % p.pid)


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


def concatenate_files(files, merged, header_str="", simple_cat=True, remove_original=False):
    logger = logging.getLogger(concatenate_files.__name__)
    logger.info("Concatenating " + " ".join(files) + " as " + merged)
    with open(merged, "w") as merged_fd:
        for index, f in enumerate(files):
            with open(f) as fd:
                if simple_cat:
                    shutil.copyfileobj(fd, merged_fd)
                else:
                    for line in fd:
                        if line.strip() and (not index or not header_str or not line.startswith(header_str)):
                            merged_fd.write(line)
            if remove_original:
                logger.info("Removing " + f)
                os.remove(f)


def check_executable(fpath):
    logger = logging.getLogger(check_executable.__name__)
    if not os.path.isfile(fpath):
        logger.error("ERROR: File %s does not exist\n" % (fpath))
        sys.exit(os.EX_NOINPUT)
    if not os.access(fpath, os.X_OK):
        logger.error("ERROR: File %s is not executable\n" % (fpath))
        sys.exit(os.EX_NOINPUT)


def fill_missing_sequences(vcf, seq_file, reference, work_dir, log_dir):
    logger = logging.getLogger(fill_missing_sequences.__name__)

    out_vcf = os.path.join(work_dir, os.path.basename(vcf))
    if out_vcf.endswith(".gz"):
        out_vcf = out_vcf[:-3]
    out_log = os.path.join(log_dir, "%s_fill_missing.log" % (os.path.basename(vcf)))

    command = ["java", "-Xmx10g", "-Xms10g", "-jar", VARSIMJAR, "randsequencevcf", "-in_vcf", vcf, "-seq", seq_file, "-out_vcf", out_vcf, "-ref", reference]
    with open(out_log, "w") as log_fd:
        logger.info("Running command " + " ".join(command))
        subprocess.check_call(" ".join(command), shell=True, stderr=log_fd)
    return out_vcf
        

def run_vcfstats(vcfs, out_dir, log_dir):
    logger = logging.getLogger(run_vcfstats.__name__)
    processes = []
    for in_vcf in vcfs:
        out_prefix = os.path.basename(in_vcf)
        vcfstats_stdout = open(os.path.join(out_dir, "%s.stats" % (out_prefix)), "w")
        vcfstats_stderr = open(os.path.join(log_dir, "%s.vcfstats.err" % (out_prefix)), "w")
        vcfstats_command = ["java", "-Xmx1g", "-Xms1g", "-jar", VARSIMJAR, "vcfstats", "-vcf",
                        in_vcf]
        p_vcfstats = subprocess.Popen(vcfstats_command, stdout=vcfstats_stdout, stderr=vcfstats_stderr)
        logger.info("Executing command " + " ".join(vcfstats_command) + " with pid " + str(p_vcfstats.pid))
        processes.append(p_vcfstats)
    return processes


def run_randvcf(sampling_vcf, out_vcf_fd, log_file_fd, seed, sex, num_snp, num_ins, num_del, num_mnp, num_complex, percent_novel, min_length, max_length, reference, prop_het):
    logger = logging.getLogger(run_randvcf.__name__)

    rand_vcf_command = ["java", "-jar", VARSIMJAR, "randvcf2vcf", "-seed", str(seed),
                        "-t", sex,
                        "-num_snp", str(num_snp), "-num_ins", str(num_ins), "-num_del",
                        str(num_del),
                        "-num_mnp", str(num_mnp), "-num_complex", str(num_complex), "-novel",
                        str(percent_novel),
                        "-min_len", str(min_length), "-max_len", str(max_length), "-ref",
                        os.path.realpath(reference),
                        "-prop_het", str(prop_het), "-vcf", sampling_vcf]

    p_rand_vcf = subprocess.Popen(rand_vcf_command, stdout=out_vcf_fd, stderr=log_file_fd)
    logger.info("Executing command " + " ".join(rand_vcf_command) + " with pid " + str(p_rand_vcf.pid))
    return p_rand_vcf


def get_version():
    return subprocess.check_output("java -jar {} -version".format(VARSIMJAR), shell=True).strip()


if __name__ == "__main__":
    check_java()

    main_parser = argparse.ArgumentParser(description="VarSim: A high-fidelity simulation validation framework",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("--out_dir", metavar="DIR",
                             help="Output directory for the simulated genome, reads and variants", required=False,
                             default="out")
    main_parser.add_argument("--work_dir", metavar="DIR", help="Work directory, currently not used", required=False,
                             default="work")
    main_parser.add_argument("--log_dir", metavar="DIR", help="Log files of all steps are kept here", required=False,
                             default="log")
    main_parser.add_argument("--reference", metavar="FASTA", help="Reference genome that variants will be inserted into",
                             required=True, type=file)
    main_parser.add_argument("--seed", metavar="seed", help="Random number seed for reproducibility", type=int, default=0)
    main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (MALE/FEMALE)", required=False, type=str,
                             choices=["MALE", "FEMALE"], default="MALE")
    main_parser.add_argument("--id", metavar="ID", help="Sample ID to be put in output VCF file", required=True)
    main_parser.add_argument("--simulator", metavar="SIMULATOR", help="Read simulator to use", required=False, type=str,
                             choices=["art", "dwgsim", "longislnd"], default="art")
    main_parser.add_argument("--simulator_executable", metavar="PATH",
                             help="Path to the executable of the read simulator chosen"
                             , required=True, type=file)
    main_parser.add_argument("--varsim_jar", metavar="PATH", help="Path to VarSim.jar (deprecated)", type=file,
                             default=DEFAULT_VARSIMJAR,
                             required=False)
    main_parser.add_argument("--read_length", metavar="LENGTH", help="Length of read to simulate", default=100, type=int)
    main_parser.add_argument("--nlanes", metavar="INTEGER",
                             help="Number of lanes to generate, coverage will be divided evenly over the lanes. Simulation is parallized over lanes. Each lane will have its own pair of files",
                             default=1, type=int)
    main_parser.add_argument("--total_coverage", metavar="FLOAT", help="Total coverage to simulate", default=1.0,
                             type=float)
    main_parser.add_argument("--mean_fragment_size", metavar="INT", help="Mean fragment size to simulate", default=350,
                             type=int)
    main_parser.add_argument("--sd_fragment_size", metavar="INT", help="Standard deviation of fragment size to simulate",
                             default=50, type=int)
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
    rand_vcf_group.add_argument("--vc_in_vcf", metavar="VCF", help="Input small variant VCF, usually dbSNP", type=file,
                                required=False)
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
                                help="Path to file containing concatenation of real insertion sequences", type=file,
                                required=False)
    rand_dgv_group.add_argument("--sv_dgv", metavar="DGV_FILE", help="DGV file containing structural variants", type=file,
                                required=False)

    dwgsim_group = main_parser.add_argument_group("DWGSIM options")
    dwgsim_group.add_argument("--dwgsim_start_e", metavar="first_base_error_rate", help="Error rate on the first base",
                              default=0.0001, type=float)
    dwgsim_group.add_argument("--dwgsim_end_e", metavar="last_base_error_rate", help="Error rate on the last base",
                              default=0.0015, type=float)
    dwgsim_group.add_argument("--dwgsim_options", help="DWGSIM command-line options", default="", required=False)

    art_group = main_parser.add_argument_group("ART options")
    art_group.add_argument("--profile_1", metavar="profile_file1", help="ART error profile for first end", default=None,
                           type=file)
    art_group.add_argument("--profile_2", metavar="profile_file2", help="ART error profile for second end", default=None,
                           type=file)
    art_group.add_argument("--art_options", help="ART command-line options", default="", required=False)

    pbsim_group = main_parser.add_argument_group("PBSIM options")
    pbsim_group.add_argument("--model_qc", metavar="model_qc", help="PBSIM QC model", default=None, type=str)

    longislnd_group = main_parser.add_argument_group("LongISLND options")
    longislnd_group.add_argument("--longislnd_options", help="LongISLND options", default="")

    args = main_parser.parse_args()

    # make the directories we need
    makedirs([args.log_dir, args.out_dir])

    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = get_loglevel(args.loglevel)
    if not args.log_to_stderr:
        logging.basicConfig(filename=os.path.join(args.log_dir, "varsim.log"), filemode="w", level=loglevel, format=FORMAT)
    else:
        logging.basicConfig(level=loglevel, format=FORMAT)
    logger = logging.getLogger(__name__)

    # Make sure we can actually execute the executable
    if not args.disable_sim:
        check_executable(args.simulator_executable.name)

    processes = []

    t_s = time.time()

    args.vcfs = map(os.path.realpath, args.vcfs)
    in_vcfs = []
    for i, vcf in enumerate(args.vcfs):
        tool_work_dir = os.path.join(args.out_dir, "filled_in", str(i))
        makedirs([tool_work_dir])
        in_vcfs.append(fill_missing_sequences(vcf, os.path.realpath(args.sv_insert_seq.name), args.reference.name, tool_work_dir, tool_work_dir))
    args.vcfs = map(os.path.realpath, in_vcfs)

    open_fds = []
    if not args.disable_rand_vcf:
        rand_vcf_out_fd = open(os.path.join(args.out_dir, "random.vc.vcf"), "w")
        rand_vcf_log_fd = open(os.path.join(args.log_dir, "RandVCF2VCF.err"), "w")
        args.vcfs.append(os.path.realpath(rand_vcf_out_fd.name))
        processes.append(run_randvcf(os.path.realpath(args.vc_in_vcf.name), rand_vcf_out_fd, rand_vcf_log_fd,
                    args.seed, args.sex, args.vc_num_snp, args.vc_num_ins, args.vc_num_del, args.vc_num_mnp,
                    args.vc_num_complex, args.vc_percent_novel, args.vc_min_length_lim, args.vc_max_length_lim,
                    args.reference.name, args.vc_prop_het))
        open_fds += [rand_vcf_out_fd, rand_vcf_log_fd]

    if not args.disable_rand_dgv:
        rand_dgv_stdout = open(os.path.join(args.out_dir, "random.sv.vcf"), "w")
        rand_dgv_stderr = open(os.path.join(args.log_dir, "RandDGV2VCF.err"), "w")
        args.vcfs.append(os.path.realpath(rand_dgv_stdout.name))

        rand_dgv_command = ["java", "-Xms10g", "-Xmx10g", "-jar", VARSIMJAR, "randdgv2vcf",
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
        open_fds += [rand_dgv_stdout, rand_dgv_stderr]

    processes = monitor_processes(processes)
    for open_fd in open_fds:
        open_fd.close()

    merged_reference = os.path.join(args.out_dir, "%s.fa" % (args.id))
    merged_truth_vcf = os.path.join(args.out_dir, "%s.truth.vcf" % (args.id))
    merged_map = os.path.join(args.out_dir, "%s.map" % (args.id))

    processes = run_vcfstats(args.vcfs, args.out_dir, args.log_dir)

    if not args.disable_vcf2diploid:
        args.vcfs.reverse()
        vcf2diploid_stdout = open(os.path.join(args.out_dir, "vcf2diploid.out"), "w")
        vcf2diploid_stderr = open(os.path.join(args.log_dir, "vcf2diploid.err"), "w")
        vcf_arg_list = sum([["-vcf", v] for v in args.vcfs], [])
        filter_arg_list = ["-pass"] if args.filter else []
        vcf2diploid_command = ["java", "-jar", VARSIMJAR, "vcf2diploid",
                               "-t", args.sex,
                               "-id", args.id,
                               "-chr", os.path.realpath(args.reference.name)] + filter_arg_list + vcf_arg_list

        p_vcf2diploid = subprocess.Popen(vcf2diploid_command, stdout=vcf2diploid_stdout, stderr=vcf2diploid_stderr,
                                         cwd=args.out_dir)
        logger.info("Executing command " + " ".join(vcf2diploid_command) + " with pid " + str(p_vcf2diploid.pid))
        processes.append(p_vcf2diploid)

        processes = monitor_processes(processes)

        # Now concatenate the .fa from vcf2diploid
        contigs = get_contigs_list(args.reference.name)
        contig_fastas = map(lambda (x, y): os.path.join(args.out_dir, "%s_%s_%s.fa" % (x, args.id, y)), itertools.product(contigs, ["maternal", "paternal"]))
        fastas_to_cat = filter(os.path.isfile, contig_fastas)
        concatenate_files(fastas_to_cat, merged_reference, remove_original=True)

        if os.path.getsize(merged_reference) == 0:
            logger.error("Merged FASTA is empty. Something bad happened. Exiting")
            raise RuntimeError("Empty FASTA generated by vcf2diploid")

        # contatenate the vcfs
        vcfs_to_cat = filter(os.path.isfile, map(lambda x: os.path.join(args.out_dir, "%s_%s.vcf" % (x, args.id)), contigs))
        concatenate_files(vcfs_to_cat, merged_truth_vcf, header_str="#", simple_cat=False, remove_original=True)

        monitor_processes(run_vcfstats([merged_truth_vcf], args.out_dir, args.log_dir))

        if args.lift_ref:
            lifted_dir = os.path.join(args.out_dir, "lifted")
            makedirs([lifted_dir])
            #quick fix for issue of CN
            convertCN([merged_truth_vcf], "two2one")
            merged_truth_vcf = lift_vcfs([merged_truth_vcf], os.path.join(lifted_dir, "truth.vcf"), None)
            #quick fix for issue of CN
            convertCN([merged_truth_vcf], "one2two")
            merged_map = lift_maps([merged_map], os.path.join(lifted_dir, "truth.map"))

    if processes:
        processes = monitor_processes(processes)

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
            for i, end in itertools.product(xrange(args.nlanes), [1, 2]):
                fifo_src_dst.append(
                    ("simulated.lane%d.read%d.fastq" % (i, end),
                     "simulated.lane%d.read%d.fq.gz" % (i, end)))
        elif args.simulator == "art":
            for i, end, suffix in itertools.product(xrange(args.nlanes), [1, 2], ["fq", "aln"]):
                fifo_src_dst.append(("simulated.lane%d.read%d.%s" % (i, end, suffix),
                                     "simulated.lane%d.read%d.%s.gz" % (i, end, suffix)))
        elif args.simulator == "pbsim":
            for i, end, suffix in itertools.product(xrange(args.nlanes), [1, 2], ["fq", "maf"]): # the '2' read files are empty and for compatibility only
                fifo_src_dst.append(("simulated.lane%d.read%d.%s" % (i, end, suffix),
                                     "simulated.lane%d.read%d.%s.gz" % (i, end, suffix)))
        elif args.simulator == "longislnd":
            pass
        else:
            raise NotImplementedError("simulation method " + args.simulator + " not implemented");

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
            with open(merged_reference, 'r') as fa:
                nRef = sum(1 for line in fa if len(line) > 0 and line[0] == '>')
            assert 0 < nRef < 10000

            for i in xrange(args.nlanes):
                tmp_prefix = os.path.join(args.out_dir, "simulated.lane%d" % (i));
                tmp_fastq_list = " ".join(
                    "%s_%s.fastq" % (tmp_prefix, "0" * (4 - len(str(idx))) + str(idx)) for idx in range(1, nRef + 1))
                tmp_maf_list = " ".join(
                    "%s_%s.maf" % (tmp_prefix, "0" * (4 - len(str(idx))) + str(idx)) for idx in range(1, nRef + 1))
                tmp_ref_list = " ".join(
                    "%s_%s.ref" % (tmp_prefix, "0" * (4 - len(str(idx))) + str(idx)) for idx in range(1, nRef + 1))
                pbsim_command = [os.path.realpath(args.simulator_executable.name),
                                 "--data-type", "CLR",
                                 "--depth", str(coverage_per_lane),
                                 "--model_qc", args.model_qc,
                                 "--seed", str(2089 * (i + 1)),
                                 merged_reference,
                                 "--prefix", tmp_prefix,
                                 "&& ( cat", tmp_fastq_list, "> %s.read1.fq" % (tmp_prefix), ")",
                                 # this cat is i/o bound, need optimization of piping
                                 "&& ( cat", tmp_maf_list, "> %s.read1.maf" % (tmp_prefix),
                                 ")"  # this cat is i/o bound, need optimization of piping
                                 "&& ( head -q -n1 ", tmp_ref_list, "> %s.ref" % (tmp_prefix),
                                 ")"  # make reference header list
                                 "&& ( echo > %s.read2.fq" % (tmp_prefix), ")",  # dummy file
                                 "&& ( echo > %s.read2.maf" % (tmp_prefix), ")"  # dummy file
                                 ]

                pbsim_command = " ".join(pbsim_command)

                pbsim_stdout = open(os.path.join(args.log_dir, "pbsim.lane%d.out" % (i)), "w")
                pbsim_stderr = open(os.path.join(args.log_dir, "pbsim.lane%d.err" % (i)), "w")
                pbsim_p = Process(target=run_shell_command, args=(pbsim_command, pbsim_stdout, pbsim_stderr))
                pbsim_p.start()
                processes.append(pbsim_p)
                logger.info("Executing command " + pbsim_command + " with pid " + str(pbsim_p.pid))
        elif args.simulator == "longislnd":
            longislnd_command = [args.simulator_executable.name, args.longislnd_options, "--coverage", str(args.total_coverage), "--out", os.path.join(args.out_dir, "longislnd_sim"), "--fasta", merged_reference]
            longislnd_command = " ".join(longislnd_command)
            longislnd_stdout = open(os.path.join(args.log_dir, "longislnd.out"), "w")
            longislnd_stderr = open(os.path.join(args.log_dir, "longislnd.err"), "w")
            longislnd_p = Process(target=run_shell_command, args=(longislnd_command, longislnd_stdout, longislnd_stderr))
            longislnd_p.start()
            processes.append(longislnd_p)
            logger.info("Executing command " + longislnd_command + " with pid " + str(longislnd_p.pid))
        else:
            raise NotImplementedError("simulation method " + args.simulator + " not implemented");

        monitor_multiprocesses(processes, logger)
        processes = []

        logger.info("Read generation took %g seconds" % (time.time() - sim_ts))

        sim_t_liftover = time.time()

        # Now start lifting over the gzipped files
        if args.simulator != "longislnd":
            for i in xrange(args.nlanes):
                liftover_stdout = open(os.path.join(args.log_dir, "lane%d.out" % (i)), "w")
                liftover_stderr = open(os.path.join(args.log_dir, "liftover%d.log" % (i)), "w")
                fastq_liftover_command = "java -server -Xms4g -Xmx4g -jar %s fastq_liftover -map %s -id %d " \
                                         "-fastq <(gunzip -c %s/simulated.lane%d.read1.fq.gz) " \
                                         "-fastq <(gunzip -c %s/simulated.lane%d.read2.fq.gz) " \
                                         "-out >(gzip -1 > %s/lane%d.read1.fq.gz) " \
                                         "-out >(gzip -1 > %s/lane%d.read2.fq.gz)" % (
                                             VARSIMJAR, merged_map, i, args.out_dir, i,
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
                                              "-ref %s/simulated.lane%d.ref " % (args.out_dir, i, args.out_dir, i)
                fastq_liftover_command = "bash -c \"%s\"" % (fastq_liftover_command)
                liftover_p = Process(target=run_shell_command, args=(fastq_liftover_command, liftover_stdout, liftover_stderr))
                liftover_p.start()
                processes.append(liftover_p)
                fastqs.append(os.path.join(args.out_dir, "lane%d.read%d.fq.gz" % (i, end)))
                logger.info("Executing command " + fastq_liftover_command + " with pid " + str(liftover_p.pid))
        else:
            # liftover the read map files
            read_map_files = list(glob.glob(os.path.join(args.out_dir, "longislnd_sim", "*.bed")))
            merged_raw_readmap = os.path.join(args.out_dir, "longislnd_sim", "merged_readmap.bed")
            concatenate_files(read_map_files, merged_raw_readmap)
            read_maps = "-longislnd %s" % merged_raw_readmap 
            read_map_liftover_command = "java -server -jar %s longislnd_liftover " % VARSIMJAR + read_maps + " -map %s " % merged_map + " -out %s" % (os.path.join(args.out_dir, args.id + ".truth.map"))
            read_map_liftover_stderr = open(os.path.join(args.log_dir, "longislnd_liftover.err"), "w")
            read_map_liftover_p = Process(target=run_shell_command, args=(read_map_liftover_command, None, read_map_liftover_stderr))
            read_map_liftover_p.start()
            processes.append(read_map_liftover_p)
            logger.info("Executing command " + read_map_liftover_command + " with pid " + str(read_map_liftover_p.pid))

        monitor_multiprocesses(processes, logger)

        logger.info("Liftover took %g seconds" % (time.time() - sim_t_liftover))

        sim_te = max(sim_ts + 1, time.time())
        bytes_written = sum([os.path.getsize(fastq) for fastq in fastqs])
        logger.info("Took %g seconds, %ld Mbytes written, %g MB/s" % (
            sim_te - sim_ts, bytes_written / 1024.0 / 1024.0, bytes_written / 1024.0 / 1024.0 / (sim_te - sim_ts)))

        for fifo in fifos:
            os.remove(fifo)

    if not args.keep_temp:
        logger.info("Cleaning up intermediate files")
        for f in tmp_files:
            os.remove(f)
    logger.info("Done! (%g hours)" % ((time.time() - t_s) / 3600.0))
