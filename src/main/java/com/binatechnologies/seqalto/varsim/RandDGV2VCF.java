package com.binatechnologies.seqalto.varsim;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.Random;

/**
 * @author johnmu
 */

public class RandDGV2VCF extends randVCFgenerator {
    private final static Logger log = Logger.getLogger(RandDGV2VCF.class.getName());


    int num_novel_added = 0;

    RandDGV2VCF() {
        super();
        num_novel_added = 0;
    }

    /**
     * @param args command line arguments
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub
        RandDGV2VCF runner = new RandDGV2VCF();
        runner.run(args);
    }

    // outputs VCF record with random phase
    void rand_output_vcf_record(BufferedWriter bw, Variant var,
                                SimpleReference ref, byte[] insert_seq, double ratio_novel,
                                Genotypes geno) throws IOException {

        int chr_idx = var.chromosome();
        int num_alt = var.get_num_alt();

        // determine whether this one is novel
        double rand_num = rand.nextDouble();
        if (rand_num <= ratio_novel) {
            // make the variant novel, simply modify it
            // TODO maybe modifying it is bad

            int chr_len = ref.getRefLen(chr_idx);
            int buffer = Math.max(
                    100000,
                    Math.max(var.max_len(geno.geno[0]),
                            var.max_len(geno.geno[1]))
            );
            int start_val = Math.min(buffer, Math.max(chr_len - buffer, 0));
            int end_val = Math.max(chr_len - buffer, Math.min(buffer, chr_len));

            int time_out = 0;
            int new_pos = rand.nextInt(end_val - start_val + 1) + start_val + 1;
            while (!var.setNovelPosition(new_pos, ref)) {
                if (time_out > 100) {
                    log.warn("Error, cannot set novel position: " + (end_val - start_val + 1));
                    log.warn(var.deletion());
                    log.warn(var);
                    //System.exit(1);
                }

                log.info(time_out + " : " + new_pos + " : " + var.deletion());

                new_pos = rand.nextInt(end_val - start_val + 1) + start_val + 1;
                time_out++;
            }

            num_novel_added++;
            var.setVarID("Novel_" + num_novel_added);
        }

        // this is ok if both the same genotype
        // the second call will return
        fill_in_seq(var, insert_seq, geno.geno[0]);
        fill_in_seq(var, insert_seq, geno.geno[1]);

        output_vcf_record(bw, var, geno.geno[0], geno.geno[1]);

    }


    public void run(String[] args) {
        String usage = "RandDGV2VCF seed num_INS num_DEL num_DUP num_INV ratio_novel min_length_lim max_length_lim "
                + "reference_file insert_seq.txt dgv_file.txt\n"
                + "     seed            -- Seed for random number generator\n"
                + "     num_INS         -- Number of insertion SV to generate\n"
                + "     num_DEL         -- Number of deletion SV to generate\n"
                + "     num_DUP         -- Number of tandem duplication SV to generate\n"
                + "     num_INV         -- Number of inversion SV to generate\n"
                + "     ratio_novel     -- Average ratio of SV that are novel [0,1]\n"
                + "     min_length_lim  -- Minimum length variant to generate (inclusive)\n"
                + "     max_length_lim  -- Maximum length variant to generate (inclusive)\n"
                + "     reference_file  -- Reference genome sequence, eg. b37\n"
                + "     insert_seq      -- Single line concatenation of known insertion sequences\n"
                + "     dgv_file        -- flat file from DGV database\n"
                + "Outputs VCF to stdout. Randomly samples variants from DGV flat file.\n";
        if (args.length != 11) {
            System.err.println(usage);
            System.exit(1);
        }

        int seed = Integer.parseInt(args[0]);
        int num_INS = Integer.parseInt(args[1]);
        int num_DEL = Integer.parseInt(args[2]);
        int num_DUP = Integer.parseInt(args[3]);
        int num_INV = Integer.parseInt(args[4]);
        double ratio_novel = Double.parseDouble(args[5]);
        int min_length_lim = Integer.parseInt(args[6]);
        int max_length_lim = Integer.parseInt(args[7]);
        String reference_filename = args[8];
        String insert_filename = args[9];
        String dgv_filename = args[10];

        if (ratio_novel > 1 || ratio_novel < 0) {
            System.err.println(usage);
            System.exit(1);
        }

        rand = new Random(seed);

        log.info("Reading reference");
        SimpleReference ref = new SimpleReference(reference_filename);

        // read in the insert sequences
        log.info("Reading insert sequences");

        byte[] insert_seq = null;

        try {
            FileReader fileReader = new FileReader(insert_filename);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            StringBuilder sb = new StringBuilder();
            String line = null;
            while ((line = bufferedReader.readLine()) != null) {
                line = line.trim();
                sb.append(line);
            }
            bufferedReader.close();
            insert_seq = sb.toString().getBytes("US-ASCII");
        } catch (IOException e) {
            e.printStackTrace();
        }

        // count the number of variants
        log.info("Counting variants and assigning genotypes");
        ArrayList<Genotypes> selected_geno = new ArrayList<Genotypes>();

        int total_num_INS = 0;
        int total_num_DEL = 0;
        int total_num_DUP = 0;
        int total_num_INV = 0;
        int total_num_other = 0;
        int total_num = 0;
        int total_lines = 0;
        int total_duplicate = 0;
        int total_out_of_range = 0;
        DGVparser parser_one = new DGVparser(dgv_filename, ref);
        Variant prev_var = new Variant();

        // Read through a first time to generate the counts for sampling without replacement
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
            total_lines++;

            if (prev_var.position() == var.position()) {
                // duplicate
                total_duplicate++;
                continue;
            }

            prev_var = var;


            if (var.max_len() > max_length_lim
                    || var.min_len() < min_length_lim) {
                total_out_of_range++;
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
                    case Insertion:
                        total_num_INS++;
                        break;
                    case Deletion:
                        total_num_DEL++;
                        break;
                    case Tandem_Duplication:
                        total_num_DUP++;
                        break;
                    case Inversion:
                        total_num_INV++;
                        break;
                    default:
                        total_num_other++;
                }

            }

            total_num++;

        }

        log.info("total_num_INS: " + total_num_INS);
        log.info("total_num_DEL: " + total_num_DEL);
        log.info("total_num_DUP: " + total_num_DUP);
        log.info("total_num_INV: " + total_num_INV);
        log.info("total_num_skipped: " + total_num_other);
        log.info("total_num: " + total_num);
        log.info("total_duplicate: " + total_duplicate);
        log.info("total_out_of_range: " + total_out_of_range);
        log.info("total_lines: " + total_lines);

        // read through DGV file again and output the new sampled VCF file
        log.info("Writing sampled variant file");

        // write the header
        System.out.print("##fileformat=VCFv4.0\n");
        System.out
                .print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        System.out
                .print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsv\n");

        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new OutputStreamWriter(System.out));
        } catch (Exception e) {
            e.printStackTrace();
        }

        Sample_params INS_params = new Sample_params();
        Sample_params DEL_params = new Sample_params();
        Sample_params DUP_params = new Sample_params();
        Sample_params INV_params = new Sample_params();

        int geno_idx = 0;
        parser_one = new DGVparser(dgv_filename, ref);
        prev_var = new Variant();


        // Read through and do the sampling
        while (parser_one.hasMoreInput()) {
            Variant var = parser_one.parseLine();
            if (var == null) {
                continue;
            }

            Genotypes geno = selected_geno.get(geno_idx);
            geno_idx++;

            if (prev_var.position() == var.position()) {
                // duplicate
                continue;
            }

            prev_var = new Variant(var);

            if (var.max_len() > max_length_lim
                    || var.min_len() < min_length_lim) {
                continue;
            }

            // sample genotypes
            boolean same_geno = false;
            if (geno.geno[0] == geno.geno[1]) {
                same_geno = true;
            }

            for (int i = 0; i < 2; i++) {
                if (i == 1 && same_geno) {
                    geno.geno[1] = geno.geno[0];
                    break;
                }
                switch (var.getType(geno.geno[i])) {
                    case Insertion:
                        geno.geno[i] = sample_genotype(geno.geno[i], INS_params,
                                num_INS, total_num_INS);
                        break;
                    case Deletion:
                        geno.geno[i] = sample_genotype(geno.geno[i], DEL_params,
                                num_DEL, total_num_DEL);
                        break;
                    case Tandem_Duplication:
                        geno.geno[i] = sample_genotype(geno.geno[i], DUP_params,
                                num_DUP, total_num_DUP);
                        break;
                    case Inversion:
                        geno.geno[i] = sample_genotype(geno.geno[i], INV_params,
                                num_INV, total_num_INV);
                        break;
                    default:
                        break;
                }
            }

            // write out variant
            try {
                if (geno.isNonRef()) {
                    rand_output_vcf_record(out, var, ref, insert_seq,
                            ratio_novel, geno);
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
