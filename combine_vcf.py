#!/usr/bin/env python

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
    augmented_t = utils.combine_vcf(augmented_t, [varsim_tp, varsim_fn], duplicate_handling_mode=utils.COMBINE_KEEP_FIRST_DUPLICATE)

    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = utils.get_loglevel(args.loglevel)
    if args.log_to_file:
        logging.basicConfig(filename=args.log_to_file, filemode="w", level=loglevel, format=FORMAT)
    else:
        logging.basicConfig(level=loglevel, format=FORMAT)

    global LOGGER
    LOGGER = logging.getLogger(__name__)
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

    with tempfile.TemporaryDirectory() as tmpdirname:
        LOGGER.info('created temporary directory', tmpdirname)
        input_vcfs = args.vcfs
        for i in range(len(args.vcfs)):
            temp_vcf = os.path.join(tmpdirname, "{}.vcf".format(i))
            shutil.copyfile(args.vcfs[i], temp_vcf)
            input_vcfs[i] = utils.sort_and_compress(temp_vcf)
        output_vcf = args.output_prefix + '.vcf'
        output_vcf = utils.combine_vcf(output_vcf, input_vcfs,
                duplicate_handling_mode = dup_mode)
    LOGGER.info('{} done'.format(output_vcf))
    return

if __name__ == "__main__":
    main_parser = argparse.ArgumentParser(description="VarSim: A high-fidelity simulation validation framework",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("--ouput_prefix", metavar="PREFIX", help="output prefix", required = True)
    main_parser.add_argument("--vcfs", metavar="VCF", help="variant calls to be evaluated", nargs="+", default=[], required = True)
    main_parser.add_argument("--mode", metavar="MODE", help="mode for duplicates", choices = ['first_duplicate','all_duplicate','no_duplicate'], default = 'first_duplicate')
    main_parser.add_argument("--loglevel", help="Set logging level", choices=["debug", "warn", "info"], default="info")

    args = main_parser.parse_args()
    process(args)
