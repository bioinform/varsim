#!/bin/sh
"exec" "`dirname $0`/opt/miniconda2/bin/python" "$0" "$@"

import argparse
import logging
import os
import sys
logger = None

def process_args(args):
    """main function
    Arguments:
        args: dictionary containing arguments
    Returns:
        None
    Raises:
        None
    """
    outfile = args.prefix + ".map"
    outfile_paternal = args.prefix + ".paternal.map"
    outfile_maternal = args.prefix + ".maternal.map"
    with open(outfile, 'w') as fh, open(outfile_paternal, 'w') as fh_paternal, open(outfile_maternal, 'w') as fh_maternal:
        for l in args.map:
            if not l.rstrip():
                continue
            fields = l.split("\t")
            # exchange source and destination chromosomes
            fields[1], fields[3] = fields[3], fields[1]
            # exchnage source and destination positions
            fields[2], fields[4] = fields[4], fields[2]
            # set INS -> DEL and DEL -> INS
            if fields[6] == 'DEL':
                fields[6] = 'INS'
            elif fields[6] == 'INS':
                fields[6] = 'DEL'
            flipped_line = "\t".join(fields)
            fh.write(flipped_line)
            if args.split_haplotype:
                if fields[3].endswith('_paternal'):
                    fh_paternal.write(flipped_line)
                if fields[3].endswith('_maternal'):
                    fh_maternal.write(flipped_line)
    logger.info("{} done.".format(outfile))
    if args.split_haplotype:
        logger.info("{} {} done.".format(outfile_paternal, outfile_maternal))
    return
if __name__ == "__main__":
    INFO_FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    DEBUG_FORMAT = '%(levelname)s %(asctime)-15s %(name)-15s %(funcName)-20s %(message)s'

    parser = argparse.ArgumentParser(description="Flip map file (for internal use)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("prefix", type = str, help = 'output prefix')
    parser.add_argument("map", type = argparse.FileType('r'), help = 'VarSim map file')
    parser.add_argument("--split_haplotype", action = 'store_true', help = 'Split destination (after flipping) *_paternal and *_maternal (if any) chromosomes into separate files')
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')
    parser.add_argument("--verbose","-v", action = 'count',  help='increase verbosity')
    args = parser.parse_args()

    if (args.verbose is not None) and (args.verbose >= 1):
        logging.basicConfig(level=logging.DEBUG, format = DEBUG_FORMAT)
    else:
        logging.basicConfig(level=logging.INFO, format= INFO_FORMAT)
    logger = logging.getLogger(__name__)

    logger.info("running {}".format(" ".join(sys.argv)))
    process_args(args)
    logger.info("Done.")
