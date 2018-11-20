#!/bin/sh
"exec" "`dirname $0`/opt/miniconda2/bin/python" "$0" "$@"

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
import pysam
from liftover_restricted_vcf_map import lift_vcfs, lift_maps
from generate_small_test_ref import gen_restricted_ref_and_vcfs 
from utils import makedirs, run_shell_command, versatile_open, get_loglevel, check_java, MY_DIR, VARSIMJAR, get_version
import utils

REQUIRE_VARSIMJAR = not os.path.isfile(VARSIMJAR)
if REQUIRE_VARSIMJAR: VARSIMJAR = None

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
        with versatile_open(name, 'r') as file_fd:
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
                            os.killpg(p.pid, signal.SIGKILL)
                        except OSError, ex:
                            logger.error("Could not kill the process " + str(p.pid))
            raise Exception('Aborting... Please check log for details.')
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


def fill_missing_sequences(vcf, id, seq_file, reference, work_dir, log_dir):
    logger = logging.getLogger(fill_missing_sequences.__name__)

    out_vcf = os.path.join(work_dir, os.path.basename(vcf))
    if out_vcf.endswith(".gz"):
        out_vcf = out_vcf[:-3]
    out_log = os.path.join(log_dir, "%s_fill_missing.log" % (os.path.basename(vcf)))

    command = ["java", utils.JAVA_XMX, "-jar", VARSIMJAR, "randsequencevcf", "-id", id, "-in_vcf", vcf, "-seq", seq_file, "-out_vcf", out_vcf, "-ref", reference]
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
        vcfstats_command = ["java", utils.JAVA_XMX, "-jar", VARSIMJAR, "vcfstats", "-vcf",
                        in_vcf]
        logger.info("Executing command " + " ".join(vcfstats_command))
        p_vcfstats = subprocess.Popen(vcfstats_command, stdout=vcfstats_stdout, stderr=vcfstats_stderr)
        logger.info(" with pid " + str(p_vcfstats.pid))
        processes.append(p_vcfstats)
    return processes


class RandVCFOptions:
    def __init__(self, num_snp, num_ins, num_del, num_mnp, num_complex, percent_novel, min_length, max_length, prop_het, num_dup = 0, num_inv = 0):
        self.num_snp = num_snp
        self.num_ins = num_ins
        self.num_del = num_del
        self.num_mnp = num_mnp
        self.num_dup = num_dup
        self.num_inv = num_inv
        self.num_complex = num_complex
        self.percent_novel = percent_novel
        self.min_length = min_length
        self.max_length = max_length
        self.prop_het = prop_het

class RandDGVOptions:
    def __init__(self, num_ins, num_del, num_dup, num_inv, percent_novel, min_length, max_length, prop_het, output_all = " "):
        self.num_ins = num_ins
        self.num_del = num_del
        self.num_dup = num_dup
        self.num_inv = num_inv
        self.percent_novel = percent_novel
        self.min_length = min_length
        self.max_length = max_length
        self.prop_het = prop_het
        self.output_all = output_all

def randdgv_options2randvcf_options(randdgv_options):
    '''
    automatically set up shared fields between RandVCFOptions and RandDGVOptions
    :param randdgv_options:
    :return: RandVCFOptions instance
    '''
    return RandVCFOptions(
        num_snp= 0,
        num_ins = randdgv_options.num_ins,
        num_del = randdgv_options.num_del,
        num_mnp = 0,
        num_complex = 0,
        percent_novel= randdgv_options.percent_novel,
        min_length = randdgv_options.min_length,
        max_length = randdgv_options.max_length,
        prop_het=randdgv_options.prop_het,
        num_dup = randdgv_options.num_dup,
        num_inv = randdgv_options.num_inv
    )

def run_randvcf(sampling_vcf, out_vcf_fd, log_file_fd, seed, sex, randvcf_options, reference):
    logger = logging.getLogger(run_randvcf.__name__)

    rand_vcf_command = ["java", utils.JAVA_XMX, "-jar", VARSIMJAR, "randvcf2vcf",
                        "-seed", str(seed),
                        "-t", sex,
                        "-num_snp", str(randvcf_options.num_snp),
                        "-num_ins", str(randvcf_options.num_ins),
                        "-num_del", str(randvcf_options.num_del),
                        "-num_mnp", str(randvcf_options.num_mnp),
                        "-num_complex", str(randvcf_options.num_complex),
                        "-num_dup", str(randvcf_options.num_dup),
                        "-num_inv", str(randvcf_options.num_inv),
                        "-novel", str(randvcf_options.percent_novel),
                        "-min_len", str(randvcf_options.min_length),
                        "-max_len", str(randvcf_options.max_length),
                        "-prop_het", str(randvcf_options.prop_het),
                        "-ref", os.path.realpath(reference),
                        "-vcf", sampling_vcf]

    logger.info("Executing command " + " ".join(rand_vcf_command))
    p_rand_vcf = subprocess.Popen(rand_vcf_command, stdout=out_vcf_fd, stderr=log_file_fd)
    logger.info(" with pid " + str(p_rand_vcf.pid))
    return p_rand_vcf


def run_randdgv(dgv_file, out_vcf_fd, log_file_fd, seed, sex, options, reference, insert_seq_file):
    logger = logging.getLogger(run_randvcf.__name__)

    rand_dgv_command = ["java", utils.JAVA_XMX, "-jar", VARSIMJAR, "randdgv2vcf",
                        "-t", sex,
                        "-seed", str(seed),
                        options.output_all,
                        "-num_ins", str(options.num_ins),
                        "-num_del", str(options.num_del),
                        "-num_dup", str(options.num_dup),
                        "-num_inv", str(options.num_inv),
                        "-novel", str(options.percent_novel),
                        "-min_len", str(options.min_length),
                        "-max_len", str(options.max_length),
                        "-prop_het", str(options.prop_het),
                        "-ref", os.path.realpath(reference),
                        "-ins", os.path.realpath(insert_seq_file),
                        "-dgv", os.path.realpath(dgv_file)]

    logger.info("Executing command " + " ".join(rand_dgv_command))
    p_rand_dgv = subprocess.Popen(rand_dgv_command, stdout=out_vcf_fd, stderr=log_file_fd)
    logger.info(" with pid " + str(p_rand_dgv.pid))

    return p_rand_dgv

def varsim_main(reference,
                simulator, # use None to disable simulation
                simulator_exe,
                total_coverage,
                variant_vcfs=[],
                sampling_vcf=None,
                dgv_file=None,
                randvcf_options=None, # use None to disable RandVCF
                randdgv_options=None, # use None to disable RandDGV 
                nlanes=1,
                simulator_options="",
                sample_id="VarSim_Sample",
                log_dir="log",
                out_dir="out",
                sv_insert_seq=None,
                seed=0,
                sex="MALE",
                remove_filtered=False,
                keep_temp=False,
                force_five_base_encoding=False,
                lift_ref=False,
                disable_vcf2diploid=False):
    check_java()

    # make the directories we need
    makedirs([log_dir, out_dir])

    logger = logging.getLogger(varsim_main.__name__)

    # Make sure we can actually execute the executable
    if simulator:
        if simulator not in ["dwgsim", "art", "longislnd"]:
            raise NotImplementedError("Simulation method {} not implemented".format(simulator))
        check_executable(simulator_exe)

    processes = []

    t_s = time.time()

    variant_vcfs = map(os.path.realpath, variant_vcfs)

    if sv_insert_seq:
        in_vcfs = []
        for i, vcf in enumerate(variant_vcfs):
            tool_work_dir = os.path.join(out_dir, "filled_in", str(i))
            makedirs([tool_work_dir])
            in_vcfs.append(fill_missing_sequences(vcf, sample_id, os.path.realpath(sv_insert_seq), reference, tool_work_dir, tool_work_dir))
        variant_vcfs = map(os.path.realpath, in_vcfs)
    else:
        logger.warn("Not filling in SV sequences since no insert sequence file provided")

    open_fds = []
    if randvcf_options:
        if not sampling_vcf:
            logger.error("Need to provide the VCF for random sampling")
            raise ValueError("Sampling VCF missing")

        rand_vcf_out_fd = open(os.path.join(out_dir, "random.vc.vcf"), "w")
        rand_vcf_log_fd = open(os.path.join(log_dir, "RandVCF2VCF.err"), "w")
        variant_vcfs.append(os.path.realpath(rand_vcf_out_fd.name))
        processes.append(run_randvcf(os.path.realpath(sampling_vcf), rand_vcf_out_fd, rand_vcf_log_fd, seed, sex, randvcf_options, reference))
        open_fds += [rand_vcf_out_fd, rand_vcf_log_fd]

    if randdgv_options:
        if not sv_insert_seq:
            raise ValueError("Need SV sequence file to fill in SV sequences")

        if not dgv_file:
            logger.error("Need to provide the DGV file for random sampling")
            raise ValueError("DGV file missing")

        rand_dgv_stdout = open(os.path.join(out_dir, "random.sv.vcf"), "w")
        rand_dgv_stderr = open(os.path.join(log_dir, "RandDGV2VCF.err"), "w")
        variant_vcfs.append(os.path.realpath(rand_dgv_stdout.name))
        processes.append(run_randdgv(dgv_file, rand_dgv_stdout, rand_dgv_stderr, seed, sex, randdgv_options, reference, sv_insert_seq))
        open_fds += [rand_dgv_stdout, rand_dgv_stderr]

    processes = monitor_processes(processes)
    for open_fd in open_fds:
        open_fd.close()

    merged_reference = os.path.join(out_dir, "%s.fa" % (sample_id))
    merged_truth_vcf = os.path.join(out_dir, "%s.truth.vcf" % (sample_id))
    merged_map = os.path.join(out_dir, "%s.map" % (sample_id))

    processes = run_vcfstats(variant_vcfs, out_dir, log_dir)

    if not disable_vcf2diploid:
        vcf2diploid_stdout = open(os.path.join(out_dir, "vcf2diploid.out"), "w")
        vcf2diploid_stderr = open(os.path.join(log_dir, "vcf2diploid.err"), "w")
        vcf_arg_list = sum([["-vcf", v] for v in variant_vcfs], [])
        filter_arg_list = ["-pass"] if remove_filtered else []
        vcf2diploid_command = ["java", utils.JAVA_XMX, "-jar", VARSIMJAR, "vcf2diploid",
                               "-t", sex,
                               "-id", sample_id,
                               "-chr", os.path.realpath(reference)] + filter_arg_list + vcf_arg_list

        logger.info("Executing command " + " ".join(vcf2diploid_command))
        p_vcf2diploid = subprocess.Popen(vcf2diploid_command, stdout=vcf2diploid_stdout, stderr=vcf2diploid_stderr,
                                         cwd=out_dir)
        logger.info(" with pid " + str(p_vcf2diploid.pid))
        processes.append(p_vcf2diploid)

        processes = monitor_processes(processes)

        # Now concatenate the .fa from vcf2diploid
        contigs = get_contigs_list(reference)
        contig_fastas = map(lambda (x, y): os.path.join(out_dir, "%s_%s_%s.fa" % (x, sample_id, y)), itertools.product(contigs, ["maternal", "paternal"]))
        fastas_to_cat = filter(os.path.isfile, contig_fastas)
        concatenate_files(fastas_to_cat, merged_reference, remove_original=True)

        if os.path.getsize(merged_reference) == 0:
            logger.error("Merged FASTA is empty. Something bad happened. Exiting")
            raise RuntimeError("Empty FASTA generated by vcf2diploid")

        # contatenate the vcfs
        vcfs_to_cat = filter(os.path.isfile, map(lambda x: os.path.join(out_dir, "%s_%s.vcf" % (x, sample_id)), contigs))
        concatenate_files(vcfs_to_cat, merged_truth_vcf, header_str="#", simple_cat=False, remove_original=True)

        monitor_processes(run_vcfstats([merged_truth_vcf], out_dir, log_dir))

        if lift_ref:
            lifted_dir = os.path.join(out_dir, "lifted")
            makedirs([lifted_dir])
            #quick fix for issue of CN
            convertCN([merged_truth_vcf], "two2one")
            merged_truth_vcf = lift_vcfs([merged_truth_vcf], os.path.join(lifted_dir, "truth.vcf"), None, tabix_index=False)
            #quick fix for issue of CN
            convertCN([merged_truth_vcf], "one2two")
            pysam.tabix_index(merged_truth_vcf, force=True, preset='vcf')
            merged_map = lift_maps([merged_map], os.path.join(lifted_dir, "truth.map"))

    if processes:
        processes = monitor_processes(processes)

    # Now generate the reads using art/pbsim/dwgsim
    tmp_files = []
    if simulator:
        fifos = []
        fastqs = []
        sim_ts = time.time()
        coverage_per_lane = total_coverage * 0.5 / nlanes
        processes = []

        fifo_src_dst = []
        if simulator == "dwgsim":
            for i, end in itertools.product(xrange(nlanes), [1, 2]):
                fifo_src_dst.append(
                    ("simulated.lane%d.read%d.fastq" % (i, end),
                     "simulated.lane%d.read%d.fq.gz" % (i, end)))
        elif simulator == "art":
            for i, end, suffix in itertools.product(xrange(nlanes), [1, 2], ["fq", "aln"]):
                fifo_src_dst.append(("simulated.lane%d.read%d.%s" % (i, end, suffix),
                                     "simulated.lane%d.read%d.%s.gz" % (i, end, suffix)))
        else: # simulator == "longislnd":
            pass

        for fifo_name, dst in fifo_src_dst:
            fifos.append(os.path.join(out_dir, fifo_name))
            if os.path.exists(fifos[-1]): os.remove(fifos[-1])
            os.mkfifo(fifos[-1])

            gzip_stderr = open(os.path.join(log_dir, "gzip.%s" % (fifo_name)), "w")
            gzip_command = "cat %s | gzip -2 > %s" % (fifos[-1], os.path.join(out_dir, dst))
            logger.info("Executing command %s" % (gzip_command) )
            gzip_p = subprocess.Popen(gzip_command, stdout = None, stderr = gzip_stderr, shell = True)
            logger.info( " with pid " + str(gzip_p.pid))
            processes.append(gzip_p)
            tmp_files.append(os.path.join(out_dir, dst))

        simulator_commands_files = []
        if simulator == "dwgsim":
            for i in xrange(nlanes):
                simulator_command = "{} {} -C {} -z {} {} {}".format(os.path.realpath(simulator_exe), simulator_options, coverage_per_lane, seed + i, merged_reference, os.path.join(out_dir, "simulated.lane%d" % (i)))
                simulator_commands_files.append((simulator_command, os.path.join(log_dir, "dwgsim.lane%d.out" % (i)), os.path.join(log_dir, "dwgsim.lane%d.err" % (i))))
        elif simulator == "art":
            for i in xrange(nlanes):
                simulator_command = "{} {} -i {} -f {} -rs {} -o {}".format(simulator_exe, simulator_options, merged_reference, coverage_per_lane, seed + i, os.path.join(out_dir, "simulated.lane%d.read" % (i)))
                simulator_commands_files.append((simulator_command, os.path.join(log_dir, "art.lane%d.out" % (i)), os.path.join(log_dir, "art.lane%d.err" % (i))))
        else: # simulator == "longislnd":
            simulator_command = "{} {} --coverage {} --out {} --fasta {}".format(simulator_exe, simulator_options, total_coverage * 0.5, os.path.join(out_dir, "longislnd_sim"), merged_reference)
            simulator_commands_files.append((simulator_command, os.path.join(log_dir, "longislnd.out"), os.path.join(log_dir, "longislnd.err")))

        simulator_fds = []
        for command, stdout, stderr in simulator_commands_files:
            stdout_fd = open(stdout, "w")
            stderr_fd = open(stderr, "w")
            process = subprocess.Popen(command, stdout=stdout_fd, stderr=stderr_fd, shell=True, close_fds=True)
            logger.info("Executing command {} with pid {}".format(command, process.pid))
            processes.append(process)
            simulator_fds += [stdout_fd, stderr_fd]

        monitor_processes(processes)

        for fd in simulator_fds:
            fd.close()

        processes = []

        logger.info("Read generation took %g seconds" % (time.time() - sim_ts))

        sim_t_liftover = time.time()

        # Now start lifting over the gzipped files
        if simulator != "longislnd":
            for i in xrange(nlanes):
                liftover_stdout = open(os.path.join(log_dir, "lane%d.out" % (i)), "w")
                liftover_stderr = open(os.path.join(log_dir, "liftover%d.log" % (i)), "w")
                fastq_liftover_command = "java -server %s -jar %s fastq_liftover -map %s -id %d " \
                                         "-fastq <(gunzip -c %s/simulated.lane%d.read1.fq.gz) " \
                                         "-fastq <(gunzip -c %s/simulated.lane%d.read2.fq.gz) " \
                                         "-out >(gzip -1 > %s/lane%d.read1.fq.gz) " \
                                         "-out >(gzip -1 > %s/lane%d.read2.fq.gz)" % (
                                             utils.JAVA_XMX, 
                                             VARSIMJAR, merged_map, i, out_dir, i,
                                             out_dir, i, out_dir, i,
                                             out_dir, i)
                if force_five_base_encoding:
                    fastq_liftover_command += " -force_five_base_encoding "
                if simulator == "art":
                    fastq_liftover_command += " -type art " \
                                              "-aln <(gunzip -c %s/simulated.lane%d.read1.aln.gz) " \
                                              "-aln <(gunzip -c %s/simulated.lane%d.read2.aln.gz)" % (
                                                  out_dir, i, out_dir, i)
                elif simulator == "pbsim":
                    fastq_liftover_command += " -type pbsim " \
                                              "-maf <(gunzip -c %s/simulated.lane%d.read1.maf.gz) " \
                                              "-ref %s/simulated.lane%d.ref " % (out_dir, i, out_dir, i)
                fastq_liftover_command = "bash -c \"%s\"" % (fastq_liftover_command)
                logger.info("Executing command " + fastq_liftover_command)
		liftover_p = subprocess.Popen(fastq_liftover_command, stdout = liftover_stdout, stderr = liftover_stderr, shell = True)
                logger.info(" with pid " + str(liftover_p.pid))
                processes.append(liftover_p)
                fastqs.append(os.path.join(out_dir, "lane%d.read%d.fq.gz" % (i, end)))
        else:
            # liftover the read map files
            read_map_files = list(glob.glob(os.path.join(out_dir, "longislnd_sim", "*.bed")))
            merged_raw_readmap = os.path.join(out_dir, "longislnd_sim", "merged_readmap.bed")
            concatenate_files(read_map_files, merged_raw_readmap)
            read_maps = "-longislnd %s" % merged_raw_readmap 
            read_map_liftover_command = "java %s -server -jar %s longislnd_liftover " % (utils.JAVA_XMX, VARSIMJAR) + read_maps + " -map %s " % merged_map + " -out %s" % (os.path.join(out_dir, sample_id + ".truth.map"))
            read_map_liftover_stderr = open(os.path.join(log_dir, "longislnd_liftover.err"), "w")
            logger.info("Executing command " + read_map_liftover_command )
            read_map_liftover_p = subprocess.Popen(read_map_liftover_command, stdout = None, stderr = read_map_liftover_stderr, shell = True)
            processes.append(read_map_liftover_p)
            logger.info(" with pid " + str(read_map_liftover_p.pid))

        monitor_processes(processes)

        logger.info("Liftover took %g seconds" % (time.time() - sim_t_liftover))

        sim_te = max(sim_ts + 1, time.time())
        bytes_written = sum([os.path.getsize(fastq) for fastq in fastqs])
        logger.info("Took %g seconds, %ld Mbytes written, %g MB/s" % (
            sim_te - sim_ts, bytes_written / 1024.0 / 1024.0, bytes_written / 1024.0 / 1024.0 / (sim_te - sim_ts)))

        for fifo in fifos:
            os.remove(fifo)

    if not keep_temp:
        logger.info("Cleaning up intermediate files")
        for f in tmp_files:
            os.remove(f)
    logger.info("Done! (%g hours)" % ((time.time() - t_s) / 3600.0))


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
                             required=True)
    main_parser.add_argument("--seed", metavar="seed", help="Random number seed for reproducibility", type=int, default=0)
    main_parser.add_argument("--sex", metavar="Sex", help="Sex of the person (MALE/FEMALE)", required=False, type=str,
                             choices=["MALE", "FEMALE"], default="MALE")
    main_parser.add_argument("--id", metavar="ID", help="Sample ID to be put in output VCF file", required=True)
    main_parser.add_argument("--simulator", metavar="SIMULATOR", help="Read simulator to use", required=False, type=str,
                             choices=["art", "dwgsim", "longislnd"], default="art")
    main_parser.add_argument("--simulator_executable", metavar="PATH",
                             help="Path to the executable of the read simulator chosen"
                             , required=True)
    main_parser.add_argument("--varsim_jar", metavar="PATH", help="Path to VarSim.jar (deprecated)",
                             default=VARSIMJAR,
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
    main_parser.add_argument("--java_max_mem", metavar="XMX", help="max java memory", default="10g", type = str)
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
    rand_vcf_group.add_argument("--vc_in_vcf", metavar="VCF", help="Input small variant VCF, usually dbSNP",
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
                                help="Path to file containing concatenation of real insertion sequences",
                                required=False)
    rand_dgv_group.add_argument("--sv_dgv", metavar="DGV_FILE", help="DGV file containing structural variants",
                                required=False)
    rand_dgv_group.add_argument("--sv_prop_het", metavar="FLOAT", help="Proportion of heterozygous structural variants",
                                default=0.6,
                                type=float)

    dwgsim_group = main_parser.add_argument_group("DWGSIM options")
    dwgsim_group.add_argument("--dwgsim_start_e", metavar="first_base_error_rate", help="Error rate on the first base",
                              default=0.0001, type=float)
    dwgsim_group.add_argument("--dwgsim_end_e", metavar="last_base_error_rate", help="Error rate on the last base",
                              default=0.0015, type=float)
    dwgsim_group.add_argument("--dwgsim_options", help="DWGSIM command-line options", default="", required=False)

    art_group = main_parser.add_argument_group("ART options")
    art_group.add_argument("--profile_1", metavar="profile_file1", help="ART error profile for first end", default="")
    art_group.add_argument("--profile_2", metavar="profile_file2", help="ART error profile for second end", default="")
    art_group.add_argument("--art_options", help="ART command-line options", default="")

    pbsim_group = main_parser.add_argument_group("PBSIM options")
    pbsim_group.add_argument("--model_qc", metavar="model_qc", help="PBSIM QC model", default=None, type=str)

    longislnd_group = main_parser.add_argument_group("LongISLND options")
    longislnd_group.add_argument("--longislnd_options", help="LongISLND options", default="")

    args = main_parser.parse_args()

    utils.JAVA_XMX = utils.JAVA_XMX + args.java_max_mem
    makedirs([args.log_dir, args.out_dir])

    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = get_loglevel(args.loglevel)
    if not args.log_to_stderr:
        logging.basicConfig(filename=os.path.join(args.log_dir, "varsim.log"), filemode="w", level=loglevel, format=FORMAT)
    else:
        logging.basicConfig(level=loglevel, format=FORMAT)

    simulator = None if args.disable_sim else args.simulator
    simulator_opts = ""
    if args.simulator == "dwgsim":
        simulator_opts = "-e {1},{2} -E {1},{2} -d {3} -s {4} -1 {5} -2 {5} {6}".format(args.dwgsim_start_e, args.dwgsim_end_e, args.mean_fragment_size, args.sd_fragment_size, args.read_length, args.dwgsim_options)
    elif args.simulator == "art":
        profile_opts = "-1 {} -2 {}".format(args.profile_1, args.profile_2) if (args.profile_1 and args.profile_2) else ""
        simulator_opts = "-p -l {} -m {} -s {} {} {}".format(args.read_length, args.mean_fragment_size, args.sd_fragment_size, profile_opts, args.art_options)
    elif args.simulator == "longislnd":
        simulator_opts = args.longislnd_options
    elif args.simulator == "pbsim":
        raise NotImplementedError("pbsim is no longer supported")

    randvcf_options = None if args.disable_rand_vcf else RandVCFOptions(args.vc_num_snp, args.vc_num_ins, args.vc_num_del, args.vc_num_mnp, args.vc_num_complex, args.vc_percent_novel, args.vc_min_length_lim, args.vc_max_length_lim, args.vc_prop_het)
    randdgv_options = None if args.disable_rand_dgv else RandDGVOptions(args.sv_num_ins, args.sv_num_del, args.sv_num_dup, args.sv_num_inv, args.sv_percent_novel, args.sv_min_length_lim, args.sv_max_length_lim, args.sv_prop_het)

    logger = logging.getLogger()
    logger.info(str(args))
    varsim_main(args.reference,
                simulator,
                args.simulator_executable,
                args.total_coverage,
                variant_vcfs=args.vcfs,
                sampling_vcf=args.vc_in_vcf,
                dgv_file=args.sv_dgv,
                randvcf_options=randvcf_options,
                randdgv_options=randdgv_options,
                nlanes=args.nlanes,
                simulator_options=simulator_opts,
                sample_id=args.id,
                log_dir=args.log_dir,
                out_dir=args.out_dir,
                sv_insert_seq=args.sv_insert_seq,
                seed=args.seed,
                sex=args.sex,
                remove_filtered=args.filter,
                keep_temp=args.keep_temp,
                force_five_base_encoding=args.force_five_base_encoding,
                lift_ref=args.lift_ref,
                disable_vcf2diploid=args.disable_vcf2diploid)
