#!/usr/bin/env python

from __future__ import division
import argparse
import os
import utils
import gzip
import math
import csv
from subprocess import Popen, PIPE
import pysam

def add_var_to_fn_fp(fn_fp_dict, call_tuple, wl_tuple, call_type, info, mapq):
    """Helper function to iterate_dirs. Add to the numbers of observed false calls in the dictionary fn_fp_dict.
    Detect if a fp is pure or not. Append AF and mapq information."""

    pure_fp = 0
    if call_tuple not in fn_fp_dict:
        fn_fp_dict[call_tuple] = [wl_tuple] +  [0, 0, 0, [], []]

    if call_type == "fn":
        fn_fp_dict[call_tuple][1] += 1

    elif call_type == "fp":
        fn_fp_dict[call_tuple][2] += 1

        if info and len(info)>2:
            if info[2] == "pure":
                pure_fp = 1
                fn_fp_dict[call_tuple][3] += 1

    #AF
    if info and info[0] != "N/A":
        fn_fp_dict[call_tuple][4].append(float(info[0]))

    #mapq
    if len(mapq) > 0:
        fn_fp_dict[call_tuple][5].extend(mapq)

    return fn_fp_dict, pure_fp


def count_number_of_var_in_file(fil, var_type=None):
    """Count the number of variants in a vcf"""

    num_lines = 0
    with utils.versatile_open(fil, "rt") as file_handle:
        for line in file_handle.readlines():
            if line[0] == "#" or not line.strip():
                continue

            line_split = line.strip().split()
            if var_type:
                if var_type == "snps":
                    if len(line_split[3]) == 1 and len(line_split[4]) == 1:
                        num_lines+=1 
                else:
                    if len(line_split[3]) != 1 or len(line_split[4]) != 1:
                        num_lines+=1
            else:
                num_lines+=1

    return num_lines


def calc_score(numerator, denominator):
    """Calculate scores and upper and lower boundaries for provided numerator and denominator"""

    if denominator > 0:
        score = float(numerator/denominator)
        sd = math.sqrt(score*(1-score)/denominator)
        delta = 3.0*sd

        if score < 1: # Use Wald formula
            CI_low = max(0, score-delta)
            CI_high = min(1, score+delta)
        else: # Use Wilson formula
            CI_low = 1.0/(1.0+3*3/denominator)
            CI_high = 1.0
    else:
        score = 'nan'
        CI_low = 'nan'
        CI_high = 'nan'
    return score, CI_low, CI_high


def calc_stats(totals, aligner, caller, sim_dir, whitelist):
    """Calculate stats such as sensitivity, precision, specificity and F1 for a set of results"""

    calc_F1 = lambda s, p : 2*s*p/(s+p)

    totals["all"] = {}
    for item in totals["snps"]:
        totals["all"][item] = totals["snps"][item] + totals["nonsnps"][item]

    out_list_all = []

    if whitelist:
        out_list_all.append(["var type", "T", "TP", "FN", "FP", "pure_FP", "TN", "sensitivity", "lower confidence interval", "upper confidence interval", "precision", "lower confidence interval", "upper confidence interval", "specificity", "lower confidence interval", "upper confidence interval", "F1", "lower confidence interval", "upper confidence interval"])
    else:
        out_list_all.append(["var type", "T", "TP", "FN", "FP", "pure_FP", "sensitivity", "lower confidence interval", "upper confidence interval", "precision", "lower confidence interval", "upper confidence interval", "F1", "lower confidence interval", "upper confidence interval"])

    for var_type in ["snps", "nonsnps", "all"]:
        out_list = [var_type]
        out_list.extend([totals[var_type]["t"], totals[var_type]["tp"], totals[var_type]["fn"], totals[var_type]["fp"], totals[var_type]["pure_fp"]])
        if whitelist:
            out_list.append(totals[var_type]["0/0 sims"]+totals[var_type]["tp"])

        #sensitivity
        sens = calc_score(totals[var_type]["t"]-totals[var_type]["fn"], totals[var_type]["t"])
        out_list.extend(sens)

        #precision
        prec = calc_score(totals[var_type]["t"]-totals[var_type]["fn"], totals[var_type]["t"]-totals[var_type]["fn"]+totals[var_type]["fp"])
        out_list.extend(prec)

        #specificity
        if whitelist:
            out_list.extend(calc_score(totals[var_type]["0/0 sims"]+totals[var_type]["tp"], totals[var_type]["0/0 sims"]+totals[var_type]["tp"]+totals[var_type]["fp"]))

        #F1
        out_list.extend(['nan' if sens[0]=='nan' or prec[0]=='nan' else calc_F1(sens[0], prec[0]), 'nan' if sens[1]=='nan' or prec[1]=='nan' else calc_F1(sens[1], prec[1]), 'nan' if sens[2]=='nan' or prec[2]=='nan' else calc_F1(sens[2], prec[2])])

        out_list_all.append(out_list)


    file_string = "stats_{}_{}.csv".format(aligner, caller)
    if whitelist:
        file_string = "wl_"+file_string   

    with open(file_string, mode='w') as outfile:
        outfile_writer = csv.writer(outfile, delimiter=',')
        for i in out_list_all:
            outfile_writer.writerow(i)


def get_mapq(var, bam):
    """return a list of the mapqs of all the reads in a bam that overlap a variant"""

    mapqs = []

    samfile = pysam.AlignmentFile(bam, "rb")

    for read in samfile.fetch(var[0], int(var[1])-1, int(var[1])):
        mapqs.append(read.mapping_quality)

    return mapqs


def iterate_dirs(fn_fp_dict, totals, files, dir_list, args, whitelist):
    """Iterate through all dirs for which stats and a report are required for. Build information and numbers for use by calc_stats and produce_report"""

    if whitelist:
        num_WL_var = {}
        num_WL_var["snps"] = count_number_of_var_in_file(whitelist, "snps")
        num_WL_var["nonsnps"] = count_number_of_var_in_file(whitelist, "nonsnps")

    for i in dir_list:
        for j in range(args.num_reps):

            pure_fp = [0, 0]

            for call_type in files:
                if whitelist:
                    file_string = "{0}.dup{1}/eval_wl_merged_{2}.dup{3}.{4}.{5}.filter/augmented_{6}.vcf.gz".format(i, j, i, j, args.read_aligner, args.variant_caller, files[call_type])
                else:
                    file_string = "{0}.dup{1}/eval_merged_{2}.dup{3}.{4}.{5}.filter/augmented_{6}.vcf.gz".format(i, j, i, j, args.read_aligner, args.variant_caller, files[call_type])

                with utils.versatile_open(os.path.join(args.sim_dir, file_string), "rt") as file_handle:

                    for line in file_handle.readlines():
                        line_split = line.strip().split()

                        if line[0] == "#":
                            continue

                        if len(line_split) > 6 and ';' in line_split[6]:
                            info = line_split[6].split(';')
                        else:
                            info = None

                        call_var = tuple([line_split[k] for k in [0, 1, 3, 4]])
                        if info and len(info) > 1 and info[1] != "N/A" and call_type == "fp":
                            wl_var = tuple(info[1].replace('/', '|').split('_'))
                        else:
                            wl_var = tuple([line_split[k].replace('/', '|') for k in [0, 1, 3, 4, 9]])

                        var_type = "snps" if len(wl_var[2]) == 1 and len(wl_var[3]) == 1 else "nonsnps"

                        totals[var_type][call_type] += 1

                        if call_type in ["tp", "t"]:
                            continue

                        #mapq
                        mapq = get_mapq(wl_var, os.path.join(args.sim_dir, "{}.dup{}/merged_{}.dup{}.{}.bam".format(i, j, i, j, args.read_aligner)))

                        fn_fp_dict, p_fp = add_var_to_fn_fp(fn_fp_dict, call_var, wl_var, call_type, info, mapq)
                        pure_fp[var_type == "nonsnps"] += p_fp

            for var_type in ["snps", "nonsnps"]:
                #Calculate number of Pure_FP and TN
                totals[var_type]["pure_fp"] += pure_fp[var_type == "nonsnps"]
                if whitelist:
                    totals[var_type]["0/0 sims"] += num_WL_var[var_type] - count_number_of_var_in_file(os.path.join(args.sim_dir,"{}.dup{}/wl_true_merged_{}.dup{}.{}.{}.filter.vcf".format(i, j, i, j, args.read_aligner, args.variant_caller)), var_type) - pure_fp[var_type == "nonsnps"]

    return fn_fp_dict, totals


def get_gene_cdot_pdot(wl_variant, whitelist):
    """Return the gene, c. and p. for the variant"""

    info = None
    gene = None
    c_dot = None
    p_dot = None

    info = variant_search(wl_variant, whitelist, 3, [7])
    if info:
        info_split = info[0].split(';')
        gene = info_split[0].split(':')[0]

        for field in info_split:
            if not c_dot:
                try:
                    c_dot = field.split('c.')[1]
                except:
                    pass
            if not p_dot:
                try:
                    p_dot = field.split('p.')[1]
                except:
                    pass
            if c_dot and p_dot:
                break

    p_dot = p_dot[:-1] if p_dot and p_dot[-1] == "|" else p_dot

    return gene, c_dot, p_dot


def get_neighbourhood(var, base, ref, vcfprimers, neighbourhood_size):
    """Return the flanking k-mers of the variant"""

    utils.write_vcf(["##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample", "{}\t{}\t.\t{}\t{}\t.\t.\t.\t.\t{}".format(*var)], "temp.vcf")
    utils.sort_and_compress("temp.vcf")

    # 5' and 3' sequences
    try:
        flanking = Popen([vcfprimers, "-f", ref, "-l", str(neighbourhood_size), "temp.vcf.gz"], stdout=PIPE, universal_newlines=True)
        s_out, s_err = flanking.communicate()
    except:
        print "Error running vcfprimers"

    three_prime=s_out.split()[-1]
    five_prime=s_out.split()[1]

    return five_prime+base+three_prime


def get_homopolymer_type_length(neighbourhood, neighbourhood_size):
    """Find the base a homopolymer consists of and it's length"""

    max_length = 0
    max_base = None
    max_position = -1
    length = 1
    for i in range(1, len(neighbourhood)):
        if neighbourhood[i].upper() == neighbourhood[i-1].upper():
            length+=1
        else:
            if (length > max_length) and (i >= neighbourhood_size):
                max_length = length
                max_base = neighbourhood[i-1].upper()
                max_position = i-1
            if i > len(neighbourhood)-neighbourhood_size:
                break
            length = 1

    return max_base, max_length, max_position-max_length+1


def get_error_rate(neighbourhood, homo_type, length, position, error_model_file, neighbourhood_size):
    """Calculate the HP error rate, if not possible calculate the kmer error rate"""

    error_rate = (None, None)

    if((position < 2) or (position+length > len(neighbourhood)-2)):
        error_rate=("Manual inspection required", "N/A")

    else:
        error_rate = calc_HP_error_rate(neighbourhood.upper(), homo_type, length, position, error_model_file)
    if not error_rate:
        error_rate = calc_kmer_error_rate(neighbourhood.upper(), neighbourhood_size, error_model_file)

    return error_rate


def search_kmer_error_rate(kmer, error_model_file):
    """Search for the kmer in the error model file and return the context specific error rate."""

    error_model = open(error_model_file, "r")
    for line in error_model:
        if "KMER_STATS" in line and kmer in line:
            matched_fields = line.split()
            break
    error_model.close()
    try:
        error_rate = 100.0-float(matched_fields[5])
    except:
        error_rate=None

    return error_rate


def calc_kmer_error_rate(neighbourhood, neighbourhood_size, error_model_file):
    """Calculate the KMER error rate"""

    error_model = open(error_model_file, "r")

    matched_fields = None
    for line in error_model:
        if "KMER_STATS" in line:
            matched_fields = line.split()
            break

    error_model.close()

    kmer_length = len(matched_fields[1])
    kmer = neighbourhood[neighbourhood_size-(kmer_length//2): neighbourhood_size-(kmer_length//2)+kmer_length]

    forward_error_rate = search_kmer_error_rate(kmer, error_model_file)
    reverse_error_rate = search_kmer_error_rate(reverse_complement(kmer), error_model_file)

    if forward_error_rate and reverse_error_rate:
        error_rate = (forward_error_rate, 'plus') if forward_error_rate >= reverse_error_rate else (reverse_error_rate, 'minus')
    else:
        error_rate = ("N/A", "N/A")
    return error_rate


def search_HP_error_rate(search_string, error_model_file, length):
    """Search for the homopolymer in the error model file and return the context specific error rate."""

    error_model = open(error_model_file, "r")

    matched_fields = None
    for line in error_model:
        if "HP_STATS" in line and search_string in line and "{}:".format(length) in line:
            matched_fields = line.split()
            break
    error_model.close()

    error_rate = None
    if matched_fields:
        for item in matched_fields:
            if item.startswith(str(length)+"-"):
                error_rate = 100.0-float(item[item.find("(")+1:item.find(")")])
                break
    return error_rate


def calc_HP_error_rate(neighbourhood, homo_type, length, position, error_model_file):
    """Calculate the HP error rate"""

    five_prime_pair = neighbourhood[position-2:position]
    three_prime_pair = neighbourhood[position+length: position+length+2]

    search_string = five_prime_pair + homo_type + three_prime_pair

    forward_error_rate = search_HP_error_rate(search_string, error_model_file, length)
    reverse_error_rate = search_HP_error_rate(reverse_complement(search_string), error_model_file, length)

    if forward_error_rate and reverse_error_rate:
        error_rate = (forward_error_rate, 'plus') if forward_error_rate >= reverse_error_rate else (reverse_error_rate, 'minus')
    else:
        error_rate = None

    return error_rate


def variant_search(variant, file_to_search, line_to_start, fields_to_return):
    """Find the variant in a vcf and return the desired fields"""

    info = None
    file_to_search_file = utils.versatile_open(file_to_search, "rt")

    for file_to_search_line in file_to_search_file.readlines()[line_to_start:]:
        file_to_search_line_split = file_to_search_line.strip().split()
        if (variant[0] == file_to_search_line_split[0] and variant[1] == file_to_search_line_split[1] and variant[2] == file_to_search_line_split[3] and variant[3] == file_to_search_line_split[4]):
            info = [file_to_search_line_split[i] for i in fields_to_return]
            break

    file_to_search_file.close()

    return info


def reverse_complement(seq):
    """Return reverse complement of provided sequence"""

    complement = {
        "A":"T",
        "T":"A",
        "G":"C",
        "C":"G",
        "a":"t",
        "t":"a",
        "g":"c",
        "c":"g",
        "N":"N",
        "n":"n"
    }

    reverse = ""
    for i in seq:
        reverse = complement[i] + reverse

    return reverse


def produce_report(fn_fp, args, whitelist):
    """Generate report on false calls made"""

    out_list_all = []
    if whitelist:
        out_list_all.append(["WL chr", "WL pos", "WL REF", "WL ALT", "WL genotype", "Call chr", "Call pos", "Call REF", "Call ALT", "Gene", "c.", "p.", "FN freq", "FP freq", "Pure FP freq", "Average AF", "Average MAPQ", "Ref neighbourhood", "Ref homopolymer", "Ref homopolymer length", "Ref context specific error rate", "Strand on which error rate is the highest", "Alt neighbourhood", "Alt homopolymer", "Alt homopolymer length", "Alt context specific error rate", "Strand on which error rate is the highest"])
    else:
        out_list_all.append(["Sim chr", "Sim pos", "Sim REF", "Sim ALT", "Sim genotype", "Call chr", "Call pos", "Call REF", "Call ALT", "FN freq", "FP freq", "Pure FP freq", "Average AF", "Average MAPQ", "Ref neighbourhood", "Ref homopolymer", "Ref homopolymer length", "Ref context specific error rate", "Strand on which error rate is the highest", "Alt neighbourhood", "Alt homopolymer", "Alt homopolymer length", "Alt context specific error rate", "Strand on which error rate is the highest"])

    for var in sorted(fn_fp.keys()):

        out_list = []
        wl_variant = fn_fp[var][0]
        out_list.extend(wl_variant)
        out_list.extend(var)

        #Get gene, c. p.
        if whitelist:
            gene, c_dot, p_dot = get_gene_cdot_pdot(wl_variant, whitelist)
            out_list.extend([gene, c_dot, p_dot])

        #output fn, fp, pure fp frequency
        out_list.extend(fn_fp[var][1:4])

        #Calc mapq, af
        for i in range(4, 6):
            if len(fn_fp[var][i]) > 0:
                out_list.append(sum(fn_fp[var][i])/len(fn_fp[var][i]))
            else:
                out_list.append("N/A")

        #REF neighbourhood
        ref_neighbourhood = get_neighbourhood(wl_variant, wl_variant[2], args.ref, args.vcfprimers, args.neighbourhood_size)
        out_list.append(ref_neighbourhood[args.neighbourhood_size-10:args.neighbourhood_size] + ' ' + wl_variant[2] + ' ' + ref_neighbourhood[len(ref_neighbourhood)-args.neighbourhood_size:len(ref_neighbourhood)-args.neighbourhood_size+10])

        #REF get homopolymer type and length
        ref_type, ref_length, ref_position = get_homopolymer_type_length(ref_neighbourhood, args.neighbourhood_size)
        out_list.extend([ref_type, ref_length])

        #REF get error rate
        ref_error, ref_strand = get_error_rate(ref_neighbourhood, ref_type, ref_length, ref_position, args.error_model, args.neighbourhood_size)
        out_list.append(ref_error)
        out_list.append(ref_strand)

        #ALT neighbourhood
        #if multiple alleles select the first
        alt = wl_variant[3].split(',')[0]
        alt_neighbourhood = get_neighbourhood(wl_variant, wl_variant[3], args.ref, args.vcfprimers, args.neighbourhood_size)
        out_list.append(alt_neighbourhood[args.neighbourhood_size-10:args.neighbourhood_size] + ' ' + wl_variant[3] + ' ' + alt_neighbourhood[len(alt_neighbourhood)-args.neighbourhood_size:len(alt_neighbourhood)-args.neighbourhood_size+10])

        #ALT get homopolymer type and length
        alt_type, alt_length, alt_position = get_homopolymer_type_length(alt_neighbourhood, args.neighbourhood_size)
        out_list.extend([alt_type, alt_length])

        #ALT get error rate
        alt_error, alt_strand = get_error_rate(alt_neighbourhood, alt_type, alt_length, alt_position, args.error_model, args.neighbourhood_size)
        out_list.append(alt_error)
        out_list.append(alt_strand)

        out_list_all.append(out_list)

    #Create csv summarizing the FN and FP calls

    file_string = "report_{}_{}.csv".format(args.read_aligner, args.variant_caller)
    if whitelist:
        file_string = "wl_"+file_string

    with open(file_string, 'w') as outfile:
        outfile_writer = csv.writer(outfile, delimiter=',')
        for i in out_list_all:
            outfile_writer.writerow(i)


def get_args():
    """Process command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Process command line aguments')
    )

    parser.add_argument(
        '--sim_dir', required=True, type=os.path.abspath,
        help='Simulation Directory'
    )

    parser.add_argument(
        '--dirs', required=True, type=str,
        help='Sample Directories'
    )

    parser.add_argument(
        '--vcfprimers', required=True, type=str,
        help='path to vcfprimers'
    )

    parser.add_argument(
        '--ref', required=True, type=str,
        help='reference sdf'
    )

    parser.add_argument(
        '--num_reps', required=True, type=int,
        help='Number of replicates'
    )

    parser.add_argument(
        '--error_model', required=True, type=os.path.abspath,
        help='Path to error model file'
    )

    parser.add_argument(
        '--read_aligner', required=True, type=str,
        help='Read aligner'
    )

    parser.add_argument(
        '--variant_caller', required=True, type=str,
        help='Variant caller'
    )

    parser.add_argument(
        '--neighbourhood_size', required=True, type=int,
        help='Size of neighbourhood on either side of a variant to search for a homopolymer. The length of the sequence search would be 2 times this value plus the length of the variant itself.'
    )

    parser.add_argument(
        '--whitelist', required=False, type=os.path.abspath,
        help='Whitelist in vcf format'
    )

    args = parser.parse_args()

    return args


def process(args):
    """Process arguments, run calc_stats and produce_report"""

    fn_fp_dict = {}
    totals = {"snps":{"fn": 0, "fp": 0, "pure_fp": 0, "tp": 0, "0/0 sims": 0, "t": 0}, "nonsnps":{"fn": 0, "fp": 0, "pure_fp": 0, "tp": 0, "0/0 sims": 0, "t": 0}}
    files = {"fn": "fn_annotated", "fp": "fp_annotated", "tp": "tp", "t": "t"}

    dir_list = args.dirs.split(':')

    whitelist = False if not args.whitelist else args.whitelist

    fn_fp_dict, totals = iterate_dirs(fn_fp_dict, totals, files, dir_list, args, whitelist)
    calc_stats(totals, args.read_aligner, args.variant_caller, args.sim_dir, whitelist)

    produce_report(fn_fp_dict, args, whitelist)

if __name__ == "__main__":
    process(get_args())
