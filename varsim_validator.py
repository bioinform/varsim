#!/usr/bin/env python

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
from distutils.version import LooseVersion
from liftover_restricted_vcf_map import lift_vcfs, lift_maps
from generate_small_test_ref import gen_restricted_ref_and_vcfs 
from varsim import varsim_main, get_version, check_java, get_loglevel, makedirs, RandVCFOptions, RandDGVOptions
import json

MY_DIR = os.path.dirname(os.path.realpath(__file__))
VARSIMJAR = os.path.realpath(os.path.join(MY_DIR, "VarSim.jar"))

def aggregate_reports(sample_reports, samples, variant_type="all"):
    if variant_type == "all":
        summary_report = sample_reports[samples[0]]["num_true_correct"]["all_data"]["sum_count"]
    elif variant_type == "snp":
        summary_report = sample_reports[samples[0]]["num_true_correct"]["data"]["SNP"]["sum_count"]
        
    keys_for_aggregation = ["_TP", "_FP", "_TN", "_FN", "_T", "fp", "tn", "fn", "t"]
    for sample in samples[1:]:
        if variant_type == "all":
            sample_report = sample_reports[sample]["num_true_correct"]["all_data"]["sum_count"]
        elif variant_type == "snp":
            sample_report = sample_reports[sample]["num_true_correct"]["data"]["SNP"]["sum_count"]

        for key in keys_for_aggregation:
            summary_report[key] += sample_report[key]

    summary_report["fdr"] = float(summary_report["fp"]) / float(summary_report["fp"] + summary_report["tp"]) * 100
    summary_report["tpr"] = float(summary_report["tp"]) / float(summary_report["t"]) * 100
    summary_report["ppv"] = 100.0 - summary_report["fdr"]

    return summary_report
    

def varsim_multi_validation(reference, regions, samples, varsim_dir, variants_dir, out_dir, vcfcompare_options=""):
    logger = logging.getLogger(varsim_multi_validation.__name__)

    bed_options = "-bed {}".format(regions) if regions else ""
    vcfcompare_options += " " + bed_options

    sample_reports = {}
    for sample in samples:
        sample_dir = os.path.join(out_dir, sample)
        makedirs([sample_dir])

        sample_truth = os.path.join(varsim_dir, sample, "out", "lifted", "truth.vcf")
        sample_called = os.path.join(variants_dir, sample, "{}.vcf".format(sample))
        sample_called = os.path.join(variants_dir, sample, "fb_0.35.vcf")
        command = "java -jar {} vcfcompare -reference {} {} -true_vcf {} -prefix {} {}".format(VARSIMJAR, reference, vcfcompare_options, sample_truth, os.path.join(sample_dir, sample), sample_called)

        with open(os.path.join(sample_dir, "vcfcompare.out"), "w") as stdout, open(os.path.join(sample_dir, "vcfcompare.err"), "w") as stderr:
            subprocess.check_call(command, shell=True, stdout=stdout, stderr=stderr)

        report_json = os.path.join(sample_dir, "{}_report.json".format(sample))
        with open(report_json) as report_json_fd:
            sample_reports[sample] = json.load(report_json_fd)

        logger.info("Generated accuracy report for {}".format(sample))
        logger.info(json.dumps(sample_reports[sample]["num_true_correct"]["all_data"]["sum_count"], indent=2))

    # Now generate summary stats

    final_report = {}
    final_report["all"] = aggregate_reports(sample_reports, samples)
    final_report["snp"] = aggregate_reports(sample_reports, samples, "snp")
    print json.dumps(final_report, indent=2)

if __name__ == "__main__":
    check_java()

    parser = argparse.ArgumentParser(description="VarSim: A high-fidelity simulation validation framework",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--out_dir", metavar="DIR",
                             help="Output directory for the simulated genome, reads and variants", required=False,
                             default="out")
    parser.add_argument("--reference", metavar="FASTA", help="Reference genome that variants will be inserted into", required=True)
    parser.add_argument("--regions", help="Regions of interest for simulation. Skip for whole genome simulation")
    parser.add_argument("--samples", help="Samples to be simulated", required=True, nargs="+")
    parser.add_argument("--varsim", help="Root directory of multi-sample truth generation", required=True)
    parser.add_argument("--variants", help="Root directory of variant calls", required=True)
    parser.add_argument("--vcfcompare_options", help="Other VCFCompare options", default="")
    parser.add_argument('--version', action='version', version=get_version())
    parser.add_argument("--loglevel", help="Set logging level", choices=["debug", "warn", "info"], default="info")

    args = parser.parse_args()

    makedirs([args.out_dir])

    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = get_loglevel(args.loglevel)
    logging.basicConfig(level=loglevel, format=FORMAT)

    varsim_multi_validation(args.reference, args.regions, args.samples, args.varsim, args.variants, args.out_dir, args.vcfcompare_options)
