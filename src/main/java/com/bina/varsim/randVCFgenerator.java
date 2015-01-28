package com.bina.varsim;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Random;


class Sample_params {
    public int added_num;
    public int seen_num;

    Sample_params() {
        added_num = 0;
        seen_num = 0;
    }
};

/**
 * Stores genotypes and provides facilities for generating random ones
 */
class Genotypes {
    public byte geno[];

    // generate random genotypes
    // make sure het-hom ration is 1.5


    /**
     * generate random genotypes, make sure het-hom ration is 1.5
     *
     * @param chr Chromosome name in number format
     * @param num_alt Number of alternate alleles possible
     * @param rand    Random number generator
     */
    Genotypes(ChrString chr, GenderType male, int num_alt, Random rand) {
        this(chr, male, num_alt, rand, 0.6);
    }

    /**
     * generate random genotypes
     *
     * @param chr  Chromosome name in number format
     * @param num_alt  Number of alternate alleles possible
     * @param rand     Random number generator
     * @param prop_het proportion heterozygous, rest are homo
     */
    Genotypes(ChrString chr, GenderType male, int num_alt, Random rand, double prop_het) {
        geno = new byte[2];

        if (chr.isHaploid(male)) {
            geno[0] = 1;
            geno[1] = 1;
        } else {
            if (rand.nextDouble() > prop_het) {
                // hom
                geno[0] = (byte) (rand.nextInt(num_alt) + 1);
                geno[1] = geno[0];
            } else {
                // het
                if (rand.nextInt(2) == 0) { // random phasing
                    geno[0] = (byte) rand.nextInt(num_alt + 1);
                    geno[1] = (byte) rand.nextInt(num_alt + 1);
                    while (geno[1] == geno[0]) {
                        geno[1] = (byte) rand.nextInt(num_alt + 1);
                    }
                } else {
                    geno[1] = (byte) rand.nextInt(num_alt + 1);
                    geno[0] = (byte) rand.nextInt(num_alt + 1);
                    while (geno[1] == geno[0]) {
                        geno[0] = (byte) rand.nextInt(num_alt + 1);
                    }
                }
            }


        }
    }

    /**
     * @param geno0 First allele
     * @param geno1 Second allele
     */
    Genotypes(byte geno0, byte geno1) {
        geno = new byte[2];
        this.geno[0] = geno0;
        this.geno[1] = geno1;
    }

    /**
     * @return false if both genotypes are zero
     */
    public boolean isNonRef() {
        return !(geno[0] == 0 && geno[1] == 0);
    }

    public String toString() {
        return "(" + String.valueOf(geno[0]) + ":" + String.valueOf(geno[1]) + ")";
    }
}


abstract public class randVCFgenerator {

    Random _rand;

    /**
     * This sets a default seed, ideally different between different runs
     */
    randVCFgenerator() {
        _rand = new Random(3333);
    }

    /**
     * @param seed User specified seed for random numbers
     */
    randVCFgenerator(long seed) {
        _rand = new Random(seed);
    }


    /**
     * This is for sampling without replacement. If the genotype is not sampled, it is replaced with zero
     *
     * @param geno       previous genotype
     * @param seen_added params for sampling
     * @param num_sample number that we want to sample
     * @param num_total  total number of such variants
     * @param output_all Don't sample, just output everything
     * @return
     */
    byte sample_genotype(byte geno, Sample_params seen_added, int num_sample,
                         int num_total, boolean output_all) {
        if (!output_all) {
            if (seen_added.added_num < num_sample) {
                double rand_num = _rand.nextDouble();
                if ((num_total - seen_added.seen_num) * rand_num >= (num_sample - seen_added.added_num)) {
                    seen_added.seen_num++;
                    return 0;
                } else {
                    seen_added.seen_num++;
                    seen_added.added_num++;
                    return geno;
                }
            }
        } else {
            return geno;
        }
        return 0;
    }


    /**
     * This is for sampling without replacement. If the genotype is not sampled, it is replaced with zero
     *
     * @param geno       previous genotype
     * @param seen_added params for sampling
     * @param num_sample number that we want to sample
     * @param num_total  total number of such variants
     * @return
     */
    byte sample_genotype(byte geno, Sample_params seen_added, int num_sample,
                         int num_total) {
        return sample_genotype(geno, seen_added, num_sample, num_total, false);
    }

    /**
     * Fills in the insertions sequence with random sequence sampled from known insertions
     *
     * @param var        Variant to be filled in
     * @param insert_seq Sequence on known insertions all concatented together
     * @param geno       allele to be filled in
     */
    void fill_in_seq(Variant var, byte[] insert_seq, int geno) {
        FlexSeq alt = var.getAlt(geno);
        if (alt != null) {
            if (alt.getType() == FlexSeq.Type.INS) {
                // if insertion sequence is not given, we fill it in
                int len = alt.length();
                byte new_seq[] = new byte[len];
                if (len > insert_seq.length) {
                    // need to randomly duplicate insertion sequence
                    int seg_len = (int) Math.ceil(insert_seq.length / 10.0);
                    for (int i = 0; i < len; i += seg_len) {
                        // choose random start loc
                        int rand_start = _rand.nextInt(new_seq.length - seg_len);
                        for (int j = i; (j < len && j < (i + seg_len)); j++) {
                            new_seq[j] = insert_seq[rand_start + j - i];
                        }
                    }
                } else {
                    int rand_start = _rand.nextInt(insert_seq.length - len);
                    for (int j = 0; j < len; j++) {
                        new_seq[j] = insert_seq[rand_start + j];
                    }
                }
                alt = new FlexSeq(new_seq);
            }

            var.setAlt(geno, alt);
        }
    }


    /**
     * Outputs a VCF record of the variant, supports structural variations
     *
     * @param bw    BufferedWriter for writing
     * @param var   The variant to output
     * @param geno0 First allele
     * @param geno1 Second allele
     * @throws IOException
     */
    void output_vcf_record(BufferedWriter bw, Variant var, int geno0, int geno1)
            throws IOException {

        // ignore ACTGN
        String ref = var.getOrig_Ref().toUpperCase();
        String alt = var.alt_string().toString().toUpperCase();

        if (!ref.matches("[ACTGN]*") || (!var.isAltACTGN())) {
            return; // don't output if it is not ACTGN
        }

        // chromosome name
        bw.write(var.getChr().toString());
        bw.write("\t");
        // start position
        bw.write(String.valueOf(var.getPos() - var.getRef_deleted().length()));
        bw.write("\t");
        // variant id
        bw.write(var.getVar_id());
        bw.write("\t");
        // ref allele
        bw.write(ref);
        bw.write("\t");
        // alt alleles
        bw.write(alt);
        bw.write("\t");
        // variant quality
        bw.write(".\t");
        // pass label
        bw.write(var.getFilter());
        bw.write("\t");
        // INFO
        // the SV info is written here
        StringBuilder sbStr = new StringBuilder();
        // assume the SV type is consistent across alleles
        if (var.getType() == Variant.OverallType.Tandem_Duplication) {
            sbStr.append("SVTYPE=DUP;");
            sbStr.append("SVLEN=");
            sbStr.append(var.getLength());
        } else if (var.getType() == Variant.OverallType.Inversion) {
            sbStr.append("SVTYPE=INV;");
            sbStr.append("SVLEN=");
            sbStr.append(var.getLength());
        } else {
            sbStr.append("SVLEN=");
            sbStr.append(var.getLength());
        }
        bw.write(sbStr.toString());
        bw.write("\t");
        // label (GT)
        bw.write("GT:CN");
        bw.write("\t");
        // the genotype
        // for this one we need to work out which one is added
        sbStr = new StringBuilder();
        sbStr.append(geno0);
        sbStr.append("|");
        sbStr.append(geno1);
        sbStr.append(":");
        sbStr.append(var.getCN(geno0));
        sbStr.append("|");
        sbStr.append(var.getCN(geno1));

        bw.write(sbStr.toString());

        bw.newLine();
    }

    /**
     * The genotype is inferred from the variant in this case
     *
     * @param bw  BufferedWriter for writing
     * @param var The variant to output
     * @throws IOException
     */
    void output_vcf_record(BufferedWriter bw, Variant var)
            throws IOException {
        output_vcf_record(bw, var, var.getGeno().geno[0], var.getGeno().geno[1]);
    }


}
