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

def check_java():
    logger = logging.getLogger(check_java.__name__)
    jv = filter(lambda x: x.startswith("java version"), subprocess.check_output("java -version", stderr=subprocess.STDOUT, shell=True).split("\n"))[0].split()[2].replace("\"", "")
    if LooseVersion(jv) < LooseVersion("1.8"):
        logger.error("VarSim requires Java 1.8 to be on the path.")
        raise EnvironmentError("VarSim requires Java 1.8 to be on the path")

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

def run_bgzip(vcf):
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
