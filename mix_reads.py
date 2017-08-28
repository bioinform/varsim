import os
import pysam
import random
import argparse


def mix_reads_bam(alignment_files, proportions, mixed, expected_basecount=0):
    # Open the readers for all
    handles = map(pysam.AlignmentFile, alignment_files)

    summed_proportions = [sum(proportions[:i+1]) for i in range(len(proportions))]

    total_basecount = 0
    out_handle = pysam.AlignmentFile(mixed, template=handles[0])
    while total_basecount < expected_basecount:
        random_value = random.random()
        for index, value in enumerate(summed_proportions):
            if random_value < value:
                read = handles[index].next()
                total_basecount += read.query_length
                out_handle.write(read)
                break

    for handle in handles:
        handle.close()
    out_handle.close()

    return mixed


def mix_reads(reads, proportions, mixed, seed=0, expected_basecount=0):
    # type: (iterable(str), iterable(float), str, int) -> object

    # Assume reads in the same format
    if not reads:
        return None

    normalized_proportions = map(lambda x: float(x)/float(sum(proportions)), proportions)

    random.seed(seed)

    if reads[0].endswith(".bam"):
        return mix_reads_bam(reads, normalized_proportions, mixed, expected_basecount=expected_basecount)

    return mixed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mix reads from FASTQs or BAMs in specified proportions", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--reads", help="Read files as BAMs or FASTQs", nargs="+", required=True, default=[])
    parser.add_argument("--proportions", help="Proportions to mix the reads in", nargs="+", required=True, default=[])
    parser.add_argument("--out", help="Output file", required=True)
    parser.add_argument("--seed", help="Random seed", default=0, type=int)
    parser.add_argument("--coverage", help="Total coverage", required=True, type=float)
    parser.add_argument("--reference", help="Reference FASTA", required=True)

    args = parser.parse_args()

    if len(args.reads) != len(args.proportions):
        raise "Length of reads not matching proportions"

    # Get reference base count
    with pysam.FastaFile()

    mix_reads(args.reads, args.proportions, args.out, )