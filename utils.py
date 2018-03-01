#!/usr/bin/env python
import os
import subprocess
import sys

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
    required_opts = {
        dwgsim : ['-e', '-E', '-d', '-s', '-1', '-2'],
        art : ['-1','-2','-p','-l','-m','-s'],
        longislnd_required_opts : ['--fasta']
    }
    if simulator and simulator_opts:
        if simulator in required_opts:
            for i in required_opts[simulator]:
                if not (i in simulator_opts):
                    raise ValueError('{0} is missing for {1}'.format(i, simulator))
        else:
            raise NotImplementedError("simulator {0} is not supported".format(simulator))
    else:
        return
