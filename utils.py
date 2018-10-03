import logging
import subprocess
import os
import sys
import pysam
from distutils.version import LooseVersion
# Check java version to make sure it is Java 8
MY_DIR = os.path.dirname(os.path.realpath(__file__))
VARSIMJAR = os.path.realpath(os.path.join(MY_DIR, "VarSim.jar"))
RTGJAR = os.path.realpath(os.path.join(MY_DIR, "RTG.jar"))
SORT_VCF = os.path.realpath(os.path.join(MY_DIR, "src","sort_vcf.sh"))

COMBINE_KEEP_ALL_DUPLICATE = 1
COMBINE_KEEP_FIRST_DUPLICATE = 2
COMBINE_KEEP_NO_DUPLICATE = 3

def count_variants(vcf):
    '''
    count number of variants
    :param vcf:
    :return:
    '''
    count = 0
    with versatile_open(vcf, 'r') as fh:
        for l in fh:
            if l.rstrip() and (not l.startswith('#')):
                count += 1
    return count

def check_java():
    logger = logging.getLogger(check_java.__name__)
    try:
        jv = subprocess.check_output("java -version", stderr=subprocess.STDOUT, shell=True)
        if "openjdk" in jv or "OpenJDK" in jv:
            raise EnvironmentError("Please replace OpenJDK with Oracle JDK")
        jv = filter(lambda x: x.startswith("java version"), jv.split("\n"))[0].split()[2].replace("\"", "")
        if LooseVersion(jv) < LooseVersion("1.8"):
            logger.error("VarSim requires Java 1.8 to be on the path.")
            raise EnvironmentError("VarSim requires Java 1.8 to be on the path")
    except subprocess.CalledProcessError:
        raise EnvironmentError("No java (>=1.8) found")

def get_version():
    return subprocess.check_output("java -jar {} -version".format(VARSIMJAR), shell=True).strip()

def run_shell_command(cmd, cmd_stdout, cmd_stderr, cmd_dir="."):
    '''
    run command (list of str or str), redirect stdout, stderr to user-specified file handles
    :param cmd:
    :param cmd_stdout:
    :param cmd_stderr:
    :param cmd_dir:
    :return:
    '''
    logger = logging.getLogger(run_shell_command.__name__)
    if type(cmd) == list:
        cmd = ' '.join(cmd)
    logger.info('running ' + cmd)
    subproc = subprocess.Popen(cmd, stdout=cmd_stdout, stderr=cmd_stderr, cwd=cmd_dir, shell=True, preexec_fn=os.setsid)
    logger.info('PID ' + str(subproc.pid))
    retcode = subproc.wait()
    if retcode != 0:
        raise Exception('{0} failed'.format(cmd))
    return(retcode)

def makedirs(dirs):
    if type(dirs) == list:
        for d in dirs:
            if not os.path.exists(d):
                os.makedirs(d)
    else:
        if not os.path.exists(dirs):
            os.makedirs(dirs)

def versatile_open(filename, mode):
    '''
    open regular file, gzipped files
    :param filename: filename string
    :param mode: mode string
    :return: file handle
    '''
    if filename.endswith('.gz'):
        import gzip
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def check_simulator_opts(simulator, simulator_opts):
    '''
    check opts string for a particular simulator
    :param simulator:
    :param opts:
    :return:
    '''
    #TODO: we need to check if the options are really required
    required_opts = {
        'dwgsim' : ['-e', '-E', '-d', '-s', '-1', '-2'],
        'art' : ['-p','-l','-m','-s'],
        'longislnd' : []
    }
    if simulator:
        if simulator in required_opts:
            for i in required_opts[simulator]:
                if not (i in simulator_opts):
                    raise ValueError('{0} is missing for {1}'.format(i, simulator))
        else:
            raise NotImplementedError("simulator {0} is not supported".format(simulator))
    else:
        return

def get_loglevel(string):
    '''
    take literal loglevel
    return loglevel defined in logging module
    :param string:
    :return:
    '''
    if string == "info":
        return logging.INFO
    if string == "warn":
        return logging.WARN
    if string == "debug":
        return logging.DEBUG
    return logging.INFO

def combine_vcf(combined_vcf, vcfs, duplicate_handling_mode = COMBINE_KEEP_ALL_DUPLICATE, gzip = True):
    '''
    combine multiple VCFs, sort, optionally remove duplicate
    :param combined_vcf:
    :param vcfs:
    :param rm_duplicate: if true, remove duplicate variants (by chr+pos+ref+alt)
    :return:
    '''
    logger = logging.getLogger(combine_vcf.__name__)
    logger.info("Merging {0}".format(" ".join(map(str, vcfs))))
    if not vcfs or len(vcfs) < 2:
        raise ValueError('at least 2 VCFs required')

    sort_command = [SORT_VCF]
    sort_command.extend(vcfs)
    gz_vcf = "{}.gz".format(combined_vcf)
    with open(combined_vcf, "w") as sorted_out:
        run_shell_command(sort_command, cmd_stdout=sorted_out, cmd_stderr=sys.stderr)
    if duplicate_handling_mode == COMBINE_KEEP_FIRST_DUPLICATE or duplicate_handling_mode == COMBINE_KEEP_NO_DUPLICATE:
        previous_line = None
        current_count = 0
        uniq_vcf = combined_vcf + '.uniq'
        with open(combined_vcf, "r") as input, open(uniq_vcf, 'w') as output:
            for l in input:
                if l.startswith('#'):
                    output.write(l)
                elif previous_line:
                    #assume no empty field
                    chr0, pos0, id0, ref0, alt0 = previous_line.rstrip().split()[0:5]
                    chr1, pos1, id1, ref1, alt1 = l.rstrip().split()[0:5]
                    if (chr0, pos0, ref0, alt0) == (chr1, pos1, ref1, alt1):
                        #duplicate
                        current_count += 1
                    else:
                        if duplicate_handling_mode == COMBINE_KEEP_FIRST_DUPLICATE or\
                                (duplicate_handling_mode == COMBINE_KEEP_NO_DUPLICATE and current_count == 1):
                            output.write(previous_line)
                        elif current_count > 1:
                            logger.debug('{0} duplicated {1} times, are discarded'.format(previous_line, current_count))
                        previous_line = l
                        current_count = 1
                else:
                    previous_line = l
                    current_count = 1
            #process last variant record
            if previous_line:
                if duplicate_handling_mode == COMBINE_KEEP_FIRST_DUPLICATE or \
                    (duplicate_handling_mode == COMBINE_KEEP_NO_DUPLICATE and current_count == 1):
                    output.write(previous_line)
        os.rename(uniq_vcf, combined_vcf)
    if gzip:
        pysam.tabix_index(combined_vcf, force=True, preset='vcf')
        return gz_vcf
    else:
        return combined_vcf


def sort_and_compress(vcf):
    '''
    sort and compress vcf and return compressed filename
    :param vcf:
    :return:
    '''
    gz_vcf = "{}.gz".format(vcf)
    sorted_vcf = "{}.sorted".format(vcf)

    sort_command = [SORT_VCF, vcf]
    with open(sorted_vcf, "w") as sorted_out:
        run_shell_command(sort_command, cmd_stdout=sorted_out, cmd_stderr=sys.stderr)
    os.rename(sorted_vcf, vcf)
    pysam.tabix_index(vcf, force=True, preset='vcf')
    return gz_vcf
