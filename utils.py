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



