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
from copy import deepcopy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
VARSIMJAR = os.path.realpath(os.path.join(MY_DIR, "VarSim.jar"))

ALL_COUNTS = ["fp", "tp", "fn", "tn", "t"]

def sum_counts(count1, count2={}, keys_for_summation=["_TP", "_FP", "_FN", "_T", "fp", "fn", "tp", "t"]):    
    for key in keys_for_summation:
        if key not in count1:
            count1[key] = 0

    for key in keys_for_summation:
        if key in count2:
            count1[key] += count2[key]

    return count1
    

ENABLED_DICT = {"all": ["SNP", "Insertion", "Deletion"], "snp": ["SNP"], "ins": ["Insertion"], "del": ["Deletion"], "indel": ["Insertion", "Deletion"]}

def aggregate_reports(sample_reports, samples, variant_type="all"):
    logger = logging.getLogger(aggregate_reports.__name__)

    if not samples:
        samples = sample_reports.keys() 

    enabled = ENABLED_DICT[variant_type]

    summary_report = sum_counts({}, {}, ALL_COUNTS)
    for sample in samples:
        sample_s = sample_reports[sample]["num_true_correct"]["data"]
        for v in enabled:
            if v not in sample_s:
                logger.error("{} missing for {}".format(v, sample))
                continue
            summary_report = sum_counts(summary_report, sample_s[v]["sum_count"])
            summary_report = sum_counts(summary_report, sample_s[v]["sum_per_base_count"], keys_for_summation=["_TN", "tn"])

    summary_report["fdr"] = (float(summary_report["fp"]) / float(summary_report["fp"] + summary_report["tp"]) * 100) if (summary_report["fp"] + summary_report["tp"] > 0) else 0
    summary_report["tpr"] = (float(summary_report["tp"]) / float(summary_report["t"]) * 100) if summary_report["t"] > 0 else 0
    summary_report["ppv"] = 100.0 - summary_report["fdr"]
    summary_report["f1"] = (2.0 * float(summary_report["ppv"]  * summary_report["tpr"]) / float(summary_report["ppv"] + summary_report["tpr"])) if (summary_report["ppv"] + summary_report["tpr"] > 0) else 0
    summary_report["spc"] = float(summary_report["tn"]) / float(summary_report["tn"] + summary_report["fp"]) if (summary_report["tn"] + summary_report["fp"] > 0) else 0

    return summary_report
    

def varsim_multi_validation(reference, regions, samples, varsim_dir, variants_dir, out_dir, vcfcompare_options="", disable_vcfcompare=False):
    logger = logging.getLogger(varsim_multi_validation.__name__)

    bed_options = "-bed {}".format(regions) if regions else ""
    vcfcompare_options += " " + bed_options

    sample_reports = {}
   
    samples_found = {}
    for sample in samples:
        sample_dir = os.path.join(out_dir, sample)
        makedirs([sample_dir])

        sample_truth = os.path.join(varsim_dir, sample, "out", "lifted", "truth.vcf")
        sample_called = os.path.join(variants_dir, sample, "{}.vcf".format(sample))
        #sample_called = os.path.join(variants_dir, sample, "fb_0.35.vcf")

        if not os.path.isfile(sample_called):
            logger.error("{} missing".format(sample_called))
            continue

        samples_found[sample] = {"truth": sample_truth, "called": sample_called}

        if not disable_vcfcompare:
            command = "java -jar {} vcfcompare -reference {} {} -true_vcf {} -prefix {} {}".format(VARSIMJAR, reference, vcfcompare_options, sample_truth, os.path.join(sample_dir, sample), sample_called)

            with open(os.path.join(sample_dir, "vcfcompare.out"), "w") as stdout, open(os.path.join(sample_dir, "vcfcompare.err"), "w") as stderr:
                subprocess.check_call(command, shell=True, stdout=stdout, stderr=stderr)

        report_json = os.path.join(sample_dir, "{}_report.json".format(sample))
        with open(report_json) as report_json_fd:
            sample_reports[sample] = json.load(report_json_fd)

        logger.info("Generated accuracy report for {}".format(sample))
        #logger.info(json.dumps(sample_reports[sample]["num_true_correct"]["all_data"]["sum_count"], indent=2))

    # Now generate summary stats

    final_report = {}
    final_report["all"] = aggregate_reports(sample_reports, None)
    final_report["snp"] = aggregate_reports(sample_reports, None, "snp")
    final_report["ins"] = aggregate_reports(sample_reports, None, "ins")
    final_report["del"] = aggregate_reports(sample_reports, None, "del")
    final_report["indel"] = aggregate_reports(sample_reports, None, "indel")
    final_report["samples"] = {}
    for sample in sample_reports:
        sample_summary = {"report": {}, "truth": samples_found[sample]["truth"], "called": samples_found[sample]["called"]}
        for key in ["all", "snp", "ins", "del", "indel"]:
            sample_summary["report"][key] = aggregate_reports(sample_reports, [sample], key)
        final_report["samples"][sample] = sample_summary

    final_report["num_samples"] = len(final_report["samples"])
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
    parser.add_argument("--disable_vcfcompare", action="store_true", help="Do not run VCFcompare if already ran")
    parser.add_argument('--version', action='version', version=get_version())
    parser.add_argument("--loglevel", help="Set logging level", choices=["debug", "warn", "info"], default="info")

    args = parser.parse_args()

    makedirs([args.out_dir])

    # Setup logging
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    loglevel = get_loglevel(args.loglevel)
    logging.basicConfig(level=loglevel, format=FORMAT)

    varsim_multi_validation(args.reference, args.regions, args.samples, args.varsim, args.variants, args.out_dir, args.vcfcompare_options, args.disable_vcfcompare)
