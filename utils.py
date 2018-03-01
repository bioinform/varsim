import logging
import subprocess
import os
from distutils.version import LooseVersion
# Check java version to make sure it is Java 8
MY_DIR = os.path.dirname(os.path.realpath(__file__))
VARSIMJAR = os.path.realpath(os.path.join(MY_DIR, "VarSim.jar"))

def check_java():
    logger = logging.getLogger(check_java.__name__)
    jv = filter(lambda x: x.startswith("java version"), subprocess.check_output("java -version", stderr=subprocess.STDOUT, shell=True).split("\n"))[0].split()[2].replace("\"", "")
    if LooseVersion(jv) < LooseVersion("1.8"):
        logger.error("VarSim requires Java 1.8 to be on the path.")
        raise EnvironmentError("VarSim requires Java 1.8 to be on the path")

def get_version():
    return subprocess.check_output("java -jar {} -version".format(VARSIMJAR), shell=True).strip()


def run_shell_command(cmd, cmd_stdout, cmd_stderr, cmd_dir="."):
    subproc = subprocess.Popen(cmd, stdout=cmd_stdout, stderr=cmd_stderr, cwd=cmd_dir, shell=True, preexec_fn=os.setsid)
    retcode = subproc.wait()
    sys.exit(retcode)


def makedirs(dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

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
