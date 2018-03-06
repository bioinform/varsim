#!/usr/bin/env python

import argparse
import os
import sys
import json
import logging
import shutil
import utils
LOGGER = None

def merge_results(outdir, varsim_tp, varsim_fn, vcfeval_tp,
                  varsim_fp, vcfeval_tp_predict):
    '''
    generate augmented TP, FN, FP
    :param varsim_tp:
    :param varsim_fn:
    :param vcfeval_tp:
    :param varsim_fp:
    :param vcfeval_tp_predict:
    :return:
    '''
    #some implementation philosiphy (subject to change)
    #retain any variant recognized by VarSim (they might not be recognized by vcfeval, e.g. <DUP>)
    #assume variants are uniquely identifiable by chr+pos+ref+alt
    #varsim_tp + vcfeval_tp = augmented_tp
    #varsim_tp + varsim_fn = T
    #T - augmented_tp = augmented_fn
    #varsim_fp - vcfeval_tp_predict = augmented_fp
    augmented_tp = os.path.join(outdir, "augmented_tp.vcf")
    augmented_t = os.path.join(outdir, "augmented_t.vcf")
    augmented_fn = os.path.join(outdir, "augmented_fn.vcf")
    augmented_fp = os.path.join(outdir, "augmented_fp.vcf")
    augmented_tp = utils.combine_vcf(augmented_tp, [varsim_tp, vcfeval_tp], duplicate_handling_mode=utils.COMBINE_KEEP_FIRST_DUPLICATE)
    augmented_t = utils.combine_vcf(augmented_t, [varsim_tp, varsim_fn], duplicate_handling_mode=utils.COMBINE_KEEP_FIRST_DUPLICATE)

    #assumption: augmented_tp is subset of augmented_t
    augmented_fn = utils.combine_vcf(augmented_fn, [augmented_t, augmented_tp], duplicate_handling_mode=utils.COMBINE_KEEP_NO_DUPLICATE)
    #assumption: vcfeval_tp_predict is subset of varsim_fp
    augmented_fp = utils.combine_vcf(augmented_fp, [varsim_fp, vcfeval_tp_predict], duplicate_handling_mode=utils.COMBINE_KEEP_NO_DUPLICATE)

    return augmented_tp, augmented_fn, augmented_fp


class VCFComparator(object):
    def __init__(self, prefix, true_vcf, reference, regions, sample, vcfs, exclude_filtered, match_geno, log_to_file, opts):
        self.prefix = prefix
        self.true_vcf = true_vcf
        self.reference = reference
        self.sample = sample
        self.vcfs = vcfs
        self.exclude_filtered = exclude_filtered
        self.match_geno = match_geno
        self.log_to_file = log_to_file
        self.regions = regions
        self.opts = opts #additional options
        self.tp,self.tp_predict,self.fp,self.fn = None, None, None, None

    def run(self):
        '''
        generate TP, FN, FP
        :return:
        '''
        pass

    def get_tp(self):
        '''
        :return: TP (based on truth) file
        '''
        if not self.tp:
            self.run()
        return self.tp

    def get_tp_predict(self):
        '''
        :return: TP (based on prediction) file
        '''
        if not self.tp_predict:
            self.run()
        return self.tp_predict

    def get_fp(self):
        '''
        :return: FP file
        '''
        if not self.fp:
            self.run()
        return self.fp

    def get_fn(self):
        '''
        :return: FN file
        '''
        if not self.fn:
            self.run()
        return self.fn

class VarSimVCFComparator(VCFComparator):
    def get_tp_predict(self):
        '''
        varsim does not generate TP based off of predictions
        :return:
        '''
        return None

    def run(self):
        '''

        :return:
        '''
        cmd = ['java', '-jar', utils.VARSIMJAR, 'vcfcompare',
           '-prefix', self.prefix, '-true_vcf',
           self.true_vcf,
           '-reference', self.reference,
           ]
        if self.exclude_filtered:
            cmd.append('-exclude_filtered')
        if self.match_geno:
            cmd.append('-match_geno')
        if self.sample:
            cmd.append('-sample')
            cmd.append(self.sample)
        if self.regions:
            cmd.append('-bed')
            cmd.append(self.regions)
        if self.opts:
            cmd.append(self.opts)
        cmd.extend(self.vcfs)

        if self.log_to_file:
            with utils.versatile_open(self.log_to_file, 'a') as logout:
                utils.run_shell_command(cmd, sys.stdout, logout)
        else:
            utils.run_shell_command(cmd, sys.stdout, sys.stderr)
        tp = self.prefix + '_TP.vcf'
        fn = self.prefix + '_FN.vcf'
        fp = self.prefix + '_FP.vcf'
        for i in (tp, fn, fp):
            if not os.path.exists(i):
                raise Exception('{0} was not generated by VarSim vcfcompare. Please check and rerun.'.format(i))
        self.tp, self.fn, self.fp = tp, fn, fp

class RTGVCFComparator(VCFComparator):
    def run(self):
        '''

        :return:
        '''
        #command example
        #rtg-tools-3.8.4-bdba5ea_install/rtg vcfeval --baseline truth.vcf.gz \
        #--calls compare1.vcf.gz -o vcfeval_split_snp -t ref.sdf --output-mode=annotate --sample xx --squash-ploidy --regions ?? \
        cmd = ['java', '-jar', utils.RTGJAR, 'vcfeval',
               '-o', self.prefix, '--baseline',
               self.true_vcf,
               '-t', self.reference,
               ]
        if not self.exclude_filtered:
            cmd.append('--all-records')
        if not self.match_geno:
            cmd.append('--squash-ploidy')
        if self.sample:
            cmd.append('--sample')
            cmd.append(self.sample)
        if self.regions:
            cmd.append('--bed-regions')
            cmd.append(self.regions)
        if self.opts:
            cmd.append(self.opts)
        if len(self.vcfs) != 1:
            raise ValueError('vcfeval only takes 1 prediction VCF and 1 truth VCF: {0}'.format(self.vcfs))
        cmd.append('--calls')
        cmd.append(self.vcfs[0])

        if self.log_to_file:
            with utils.versatile_open(self.log_to_file, 'a') as logout:
                utils.run_shell_command(cmd, sys.stdout, logout)
        else:
            utils.run_shell_command(cmd, sys.stdout, sys.stderr)
        tp = os.path.join(self.prefix, 'tp-baseline.vcf.gz')
        tp_predict = os.path.join(self.prefix, 'tp.vcf.gz')
        fn = os.path.join(self.prefix, 'fn.vcf.gz')
        fp = os.path.join(self.prefix, 'fp.vcf.gz')
        for i in (tp, tp_predict, fn, fp):
            if not os.path.exists(i):
                raise Exception('{0} was not generated by vcfeval. Please check and rerun.'.format(i))
        self.tp, self.tp_predict, self.fn, self.fp = tp, tp_predict, fn, fp

def generate_sdf(reference, log):
    '''
    take reference and generate SDF
    :param reference:
    :return:
    '''
    sdf = reference + '.sdf'
    if os.path.exists(sdf):
        LOGGER.info('{0} exists, doing nothing'.format(sdf))
        LOGGER.info('to rerun SDF generation, please remove or rename {0}'.format(sdf))
        return sdf
    cmd = ['java','-jar',utils.RTGJAR,'format',
           '-o', sdf, reference]
    if log:
        with utils.versatile_open(log, 'a') as logout:
            utils.run_shell_command(cmd, logout, logout)
    else:
        utils.run_shell_command(cmd, sys.stdout, sys.stderr)
    return sdf

def process(args):
    '''
    main
    :param args:
    :return:
    '''

    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = utils.get_loglevel(args.loglevel)
    if args.log_to_file:
        logging.basicConfig(filename=args.log_to_file, filemode="w", level=loglevel, format=FORMAT)
    else:
        logging.basicConfig(level=loglevel, format=FORMAT)

    if len(args.vcfs) > 1:
        raise NotImplementedError('right now only support one prediction VCF. Quick workaround: src/sort_vcf.sh vcf1 vcf2 > merged.vcf')

    global LOGGER
    LOGGER = logging.getLogger(__name__)
    LOGGER.info('working hard ...')

    args.out_dir = os.path.abspath(args.out_dir)
    args.reference = os.path.abspath(args.reference)
    utils.makedirs([args.out_dir])

    varsim_prefix = os.path.join(args.out_dir, 'varsim_compare_results')
    varsim_comparator = VarSimVCFComparator(prefix=varsim_prefix, true_vcf = args.true_vcf, reference = args.reference,
                                            regions = args.regions,
               sample = args.sample, vcfs = args.vcfs,
               exclude_filtered = args.exclude_filtered,
               match_geno = args.match_geno, log_to_file= args.log_to_file, opts = args.vcfcompare_options)
    varsim_tp, varsim_fn, varsim_fp = varsim_comparator.get_tp(), varsim_comparator.get_fn(), varsim_comparator.get_fp()
    varsim_tp = utils.run_bgzip(varsim_tp)
    varsim_fn = utils.run_bgzip(varsim_fn)
    varsim_fp = utils.run_bgzip(varsim_fp)
    #run vcfeval
    sdf = args.sdf
    if not sdf:
        LOGGER.info("user did not supply SDF-formatted reference, trying to generate one...")
        sdf = generate_sdf(args.reference, args.log_to_file)

    '''for vcfeval
    sample column must be present, and not empty
    if single-sample vcf, vcfeval doesn't check if samples match in truth and call
    in multi-sample vcf, sample name must be specified
    right now
    '''
    vcfeval_prefix = os.path.join(args.out_dir, 'vcfeval_compare_results')
    if os.path.exists(vcfeval_prefix):
        LOGGER.warn('{0} exists, removing ...'.format(vcfeval_prefix))
        shutil.rmtree(vcfeval_prefix)
    vcfeval_comparator = RTGVCFComparator(prefix=vcfeval_prefix, true_vcf = varsim_fn, reference = sdf,
                                          regions = args.regions,
                                            sample = args.sample, vcfs = [varsim_fp],
                                            exclude_filtered = args.exclude_filtered,
                                            match_geno = args.match_geno, log_to_file= args.log_to_file,
                                          opts = args.vcfeval_options)
    vcfeval_tp, vcfeval_tp_predict = vcfeval_comparator.get_tp(), vcfeval_comparator.get_tp_predict()
    augmented_tp, augmented_fn, augmented_fp = merge_results(
                      outdir = args.out_dir,
                      varsim_tp = varsim_tp, varsim_fn = varsim_fn,
                      vcfeval_tp = vcfeval_tp, varsim_fp = varsim_fp, vcfeval_tp_predict = vcfeval_tp_predict)
    LOGGER.info("Variant comparison done.\nTrue positive: {0}\nFalse negative: {1}\nFalse positive: {2}\n".
                format(augmented_tp, augmented_fn, augmented_fp))
    summarize_results(os.path.join(args.out_dir,"augmetned"),augmented_tp, augmented_fn, augmented_fp)


def print_stats(stats):
    print ("{0: <15}\t{1: <10}\t{2: <10}\t{3: <10}\t{4: <5}\t{5: <5}\t{6: <5}".format("VariantType","Recall","Precision","F1", "TP","T", "FP"))
    for vartype, value in stats.iteritems():
        try:
            recall = value['tp'] / float(value['t']) if float(value['t']) != 0 else float('NaN')
            precision = float(value['tp']) / (value['tp'] + value['fp']) if value['tp'] + value['fp'] != 0 else float('NaN')
            f1 = 'NA' if recall == float('NaN') or precision == float('NaN') or (recall + precision) == 0 else 2 * recall * precision / (recall + precision)
        except ValueError:
            sys.stderr.write("invalide values\n")
        print ("{0: <15}\t{1:.10f}\t{2:.10f}\t{3:.10f}\t{4:<5}\t{5:<5}\t{6: <5}".format(vartype, recall, precision, f1, value['tp'], value['t'], value['fp']))

def parse_jsons(jsonfile, stats, include_sv = False):
    var_types = stats.keys()
    metrics = stats[var_types[0]].keys()
    with utils.versatile_open(jsonfile, 'r') as fh:
        data = json.load(fh)
        for vt in var_types:
            if vt in data['num_true_correct']['data']:
                for mt in metrics:
                    try:
                        stats[vt][mt] += data['num_true_correct']['data'][vt]['sum_count'][mt]
                        if include_sv:
                            stats[vt][mt] -= data['num_true_correct']['data'][vt]['svSumCount'][mt]
                    except KeyError as err:
                        print ("error in {}. No {} field".format(jsonfile, err))
                        stats[vt][mt] += 0

def summarize_results(prefix, tp, fn, fp):
    '''
    count variants by type and tabulate
    :param augmented_tp:
    :param augmented_fn:
    :param augmented_fp:
    :return:
    '''
    cmd = ['java', '-jar', utils.VARSIMJAR, 'vcfcompareresultsparser',
           '-prefix', prefix, '-tp',tp,
           '-fn', fn, '-fp', fp,
           ]
    utils.run_shell_command(cmd, cmd_stdout=sys.stdout, cmd_stderr=sys.stderr)
    jsonfile = "{0}_report.json".format(prefix)
    var_types = ['SNP', 'Deletion', 'Insertion', 'Complex']
    metrics = ['tp', 'fp', 't', 'fn']
    stats = {k: {ii: 0 for ii in metrics} for k in var_types}
    parse_jsons(jsonfile, stats, include_sv=False)
    print_stats(stats)
    return


if __name__ == "__main__":
    utils.check_java()

    main_parser = argparse.ArgumentParser(description="VarSim: A high-fidelity simulation validation framework",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument("--reference", metavar="FASTA", help="reference filename", required=True, type=str)
    main_parser.add_argument("--sdf", metavar="SDF", help="SDF formatted reference folder", required=False, type=str, default='')
    main_parser.add_argument("--out_dir", metavar="OUTDIR", help="output folder", required=True, type=str)
    main_parser.add_argument("--vcfs", metavar="VCF", help="variant calls to be evaluated", nargs="+", default=[], required = True)
    main_parser.add_argument("--true_vcf", metavar="VCF", help="Input small variant sampling VCF, usually dbSNP", required = True)
    main_parser.add_argument("--regions", help="BED file to restrict analysis [Optional]", required = False, type=str)
    main_parser.add_argument("--sample", metavar = "SAMPLE", help="sample name", required = False, type=str)
    main_parser.add_argument("--exclude_filtered", action = 'store_true', help="only consider variants with PASS or . in FILTER column", required = False)
    main_parser.add_argument("--match_geno", action = 'store_true', help="compare genotype in addition to alleles", required = False)
    main_parser.add_argument('--version', action='version', version=utils.get_version())
    main_parser.add_argument("--log_to_file", metavar="LOGFILE", help="logfile. If not specified, log to stderr", required=False, type=str, default="")
    main_parser.add_argument("--loglevel", help="Set logging level", choices=["debug", "warn", "info"], default="info")
    main_parser.add_argument("--vcfcompare_options", metavar="OPT", help="additional options for VarSim vcfcompare", default="", type = str)
    main_parser.add_argument("--vcfeval_options", metavar="OPT", help="additional options for RTG vcfeval", default="", type = str)

    args = main_parser.parse_args()
    process(args)
