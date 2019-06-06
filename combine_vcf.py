#!/bin/sh
"exec" "`dirname $0`/opt/miniconda2/bin/python" "$0" "$@"

import argparse
import os
import sys
import logging
import utils
import tempfile, shutil
LOGGER = None

def process(args):
    '''
    main
    :param args:
    :return:
    '''
    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = utils.get_loglevel(args.loglevel)
    logging.basicConfig(level=loglevel, format=FORMAT)

    global LOGGER
    LOGGER = logging.getLogger(__name__)
    LOGGER.info('running {}'.format(' '.join(sys.argv)))
    LOGGER.info('sort-and-index vcfs and then combine')

    dup_mode = None
    if args.mode == 'first_duplicate':
        dup_mode = utils.COMBINE_KEEP_FIRST_DUPLICATE
    elif args.mode == 'all_duplicate':
        dup_mode = utils.COMBINE_KEEP_ALL_DUPLICATE
    elif args.mode == 'no_duplicate':
        dup_mode = utils.COMBINE_KEEP_NO_DUPLICATE
    else:
        raise ValueError

    """
    scenarios:
    vcf
    vcf.gz
    vcf.gz + vcf.gz.tbi
    """
    input_vcfs = args.vcfs
    for i in range(len(args.vcfs)):
        current_vcf = args.vcfs[i]
        if current_vcf.endswith(".gz") and os.path.isfile(current_vcf + ".tbi"):
            input_vcfs[i] = current_vcf
        elif current_vcf.endswith(".gz"):
            LOGGER.info('indexing {}'.format(current_vcf))
            utils.index_vcf_gz(current_cvf)
            input_vcfs[i] = current_vcf
        else:
            LOGGER.info('sort and index {}'.format(current_vcf))
            input_vcfs[i] = utils.sort_and_compress(current_vcf, mode = 2)
    output_vcf = args.output_prefix + '.vcf'
    output_vcf = utils.combine_vcf(output_vcf, input_vcfs,
            duplicate_handling_mode = dup_mode)
    LOGGER.info('{} done'.format(output_vcf))
    return

if __name__ == "__main__":
    main_parser = argparse.ArgumentParser(description="VarSim: A high-fidelity simulation validation framework",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("--output_prefix", metavar="PREFIX", help="output prefix", required = True)
    main_parser.add_argument("--vcfs", metavar="VCF", help="variant calls to be evaluated", nargs="+", default=[], required = True)
    main_parser.add_argument("--mode", metavar="MODE", help="mode for duplicates", choices = ['first_duplicate','all_duplicate','no_duplicate'], default = 'first_duplicate')
    main_parser.add_argument("--loglevel", help="Set logging level", choices=["debug", "warn", "info"], default="info")

    args = main_parser.parse_args()
    process(args)
