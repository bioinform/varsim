package com.binatechnologies.seqalto.varsim;

import org.apache.log4j.Logger;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Random;
/**
 *
 */

/**
 * @author johnmu
 */


public class RandVCF2VCF extends randVCFgenerator {
    private final static Logger log = Logger.getLogger(RandVCF2VCF.class.getName());

    int num_novel_added;

    RandVCF2VCF() {
        super();
        num_novel_added = 0;
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub
        RandVCF2VCF runner = new RandVCF2VCF();
        runner.run(args);
    }

    // outputs VCF record with random phase
    void rand_output_vcf_record(BufferedWriter bw, Variant var,
                                SimpleReference ref, double ratio_novel, Genotypes geno)
            throws IOException {

        int chr_idx = var.chromosome();

        // determine whether this one is novel
        double rand_num = rand.nextDouble();
        if (rand_num <= ratio_novel) {
            // make the variant novel, simply modify it
            // TODO maybe modifying it is bad

            int chr_len = ref.getRefLen(chr_idx);
            int buffer = Math.max(
                    100000,
                    Math.max(var.max_len(geno.geno[0]),
                            var.max_len(geno.geno[1])));
            int start_val = Math.min(buffer, Math.max(chr_len - buffer, 0));
            int end_val = Math.max(chr_len - buffer, Math.min(buffer, chr_len));

            int time_out = 0;
            while (!var.setNovelPosition(rand.nextInt(end_val - start_val + 1)
                    + start_val + 1, ref)) {
                if (time_out > 100) {
                    log.warn("Error, cannot set novel position: " + (end_val - start_val + 1));
                    log.warn(var.deletion());
                    log.warn(var);
                    break;
                }

                log.info(time_out + " : " + var.deletion());

                time_out++;
            }
            num_novel_added++;
            var.setVarID("Novel_" + num_novel_added);
        }

        output_vcf_record(bw, var, geno.geno[0], geno.geno[1]);

    }

    public void run(String[] args) {
        String usage = "RandVCF2VCF seed num_SNP num_INS num_DEL num_MNP num_COMPLEX ratio_novel"
                + " min_length_lim max_length_lim reference_file file.vcf\n"
                + "     seed            -- Seed for random number generator\n"
                + "     num_SNP         -- Number of SNPs to generate\n"
                + "     num_INS         -- Number of insertions to generate\n"
                + "     num_DEL         -- Number of deletions to generate\n"
                + "     num_MNP         -- Number of MNPs to generate\n"
                + "     num_COMPLEX     -- Number of complex/mixed variants to generate\n"
                + "     ratio_novel     -- Average ratio of SV that are novel [0,1]\n"
                + "     min_length_lim  -- Minimum length variant to generate (inclusive)\n"
                + "     max_length_lim  -- Maximum length variant to generate (inclusive)\n"
                + "     reference_file  -- Reference genome sequence, eg. b37\n"
                + "     file.vcf        -- Input VCF file to sample from\n"
                + "Outputs VCF to stdout. Randomly samples variants from VCF file.";

        if (args.length != 11) {
            System.err.println(usage);
            System.exit(1);
        }

        int seed = Integer.parseInt(args[0]);
        int num_SNP = Integer.parseInt(args[1]);
        int num_INS = Integer.parseInt(args[2]);
        int num_DEL = Integer.parseInt(args[3]);
        int num_MNP = Integer.parseInt(args[4]);
        int num_COMPLEX = Integer.parseInt(args[5]);
        double ratio_novel = Double.parseDouble(args[6]);
        int min_length_lim = Integer.parseInt(args[7]);
        int max_length_lim = Integer.parseInt(args[8]);
        String reference_filename = args[9];
        String vcf_filename = args[10];

        if (ratio_novel > 1 || ratio_novel < 0) {
            System.err.println(usage);
            System.exit(1);
        }

        SimpleReference ref = new SimpleReference(reference_filename);

        rand = new Random(seed);

        // read through VCF file and count the
        log.info("Counting variants and assigning genotypes");
        ArrayList<Genotypes> selected_geno = new ArrayList<Genotypes>();

        int total_num_SNP = 0;
        int total_num_INS = 0;
        int total_num_DEL = 0;
        int total_num_MNP = 0;
        int total_num_COMPLEX = 0;
        int total_num_other = 0;
        int total_num = 0;
        VCFparser parser_one = new VCFparser(vcf_filename, null, false);
        Variant prev_var = new Variant();

        // read though once to count the totals, this is so we don't have
        // to store an array of the variants for sampling without replacement
        while (parser_one.hasMoreInput()) {
            Variant var = parser_one.parseLine();
            if (var == null) {
                continue;
            }

            // select genotypes here
            int chr_idx = var.chromosome();
            int num_alt = var.get_num_alt();

            Genotypes geno = new Genotypes(chr_idx, num_alt, rand);
            selected_geno.add(geno);

            if (prev_var.equals(var)) {
                continue;
            }

            // this is ok because var is not changed
            prev_var = var;

            if (var.max_len() > max_length_lim
                    || var.min_len() < min_length_lim) {
                continue;
            }

            boolean same_geno = false;
            if (geno.geno[0] == geno.geno[1]) {
                same_geno = true;
            }

            for (int i = 0; i < 2; i++) {
                if (i == 1 && same_geno) {
                    break;
                }
                switch (var.getType(geno.geno[i])) {
                    case SNP:
                        total_num_SNP++;
                        break;
                    case Insertion:
                        total_num_INS++;
                        break;
                    case Deletion:
                        total_num_DEL++;
                        break;
                    case MNP:
                        total_num_MNP++;
                        break;
                    case Complex:
                        total_num_COMPLEX++;
                        break;
                    default:
                        log.error(i + ":" + geno.geno[i] + ":" + var.getType(geno.geno[i]) + " : OTHER: " + var);
                        total_num_other++;
                }

            }
            total_num++;
        }

        log.info("total_num_SNP: " + total_num_SNP);
        log.info("total_num_INS: " + total_num_INS);
        log.info("total_num_DEL: " + total_num_DEL);
        log.info("total_num_MNP: " + total_num_MNP);
        log.info("total_num_COMPLEX: " + total_num_COMPLEX);
        log.info("total_num_skipped: " + total_num_other);
        log.info("total_num: " + total_num);

        // read through VCF file again and output the new sampled VCF file
        log.info("Writing sampled variant file");

        // write the header
        System.out.print("##fileformat=VCFv4.0\n");
        System.out
                .print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        System.out
                .print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsv\n");

        Sample_params SNP_params = new Sample_params();
        Sample_params INS_params = new Sample_params();
        Sample_params DEL_params = new Sample_params();
        Sample_params MNP_params = new Sample_params();
        Sample_params COMPLEX_params = new Sample_params();

        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new OutputStreamWriter(System.out));
        } catch (Exception e) {
            e.printStackTrace();
        }

        int geno_idx = 0;

        parser_one = new VCFparser(vcf_filename, null, false);
        prev_var = new Variant();

        // Read through it a second time, this time we do the sampling
        while (parser_one.hasMoreInput()) {
            Variant var = parser_one.parseLine();
            if (var == null) {
                // System.err.println("Bad variant or not a variant line");
                continue;
            }

            Genotypes geno = selected_geno.get(geno_idx);
            geno_idx++;

            if (prev_var.equals(var)) {
                // duplicate
                continue;
            }

            prev_var = var;

            if (var.max_len() > max_length_lim
                    || var.min_len() < min_length_lim) {
                continue;
            }

            // sample genotypes
            boolean same_geno = false;
            if (geno.geno[0] == geno.geno[1]) {
                same_geno = true;
            }

            // this samples genotypes from each allele individually
            for (int i = 0; i < 2; i++) {
                if (i == 1 && same_geno) {
                    geno.geno[1] = geno.geno[0];
                    break;
                }
                switch (var.getType(geno.geno[i])) {
                    case SNP:
                        geno.geno[i] = sample_genotype(geno.geno[i], SNP_params, num_SNP,
                                total_num_SNP);
                        break;
                    case Insertion:
                        geno.geno[i] = sample_genotype(geno.geno[i], INS_params, num_INS,
                                total_num_INS);
                        break;
                    case Deletion:
                        geno.geno[i] = sample_genotype(geno.geno[i], DEL_params, num_DEL,
                                total_num_DEL);
                        break;
                    case MNP:
                        geno.geno[i] = sample_genotype(geno.geno[i], MNP_params, num_MNP,
                                total_num_MNP);
                        break;
                    case Complex:
                        geno.geno[i] = sample_genotype(geno.geno[i], COMPLEX_params, num_COMPLEX,
                                total_num_COMPLEX);
                        break;
                    default:
                        break;
                }

            }

            // write out variant
            try {
                if (geno.isNonRef()) {
                    rand_output_vcf_record(out, var, ref, ratio_novel, geno);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        try {
            out.flush();
            out.close();
        } catch (IOException e) {
            log.error(e);
        }
    }

}
