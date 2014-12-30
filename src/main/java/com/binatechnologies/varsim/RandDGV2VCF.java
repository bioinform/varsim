package com.binatechnologies.varsim;

import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import java.io.*;
import java.util.ArrayList;
import java.util.Random;

/**
 * @author johnmu
 */

public class RandDGV2VCF extends randVCFgenerator {
    private final static Logger log = Logger.getLogger(RandDGV2VCF.class.getName());

    @Option(name = "-all", usage = "Output all variants")
    boolean output_all;

    static final int SEED_ARG = 333;
    @Option(name = "-seed", usage = "Seed for random sampling ["+SEED_ARG+"]")
    static int seed = SEED_ARG;

    static final int NUM_INS_ARG = 2000;
    @Option(name = "-num_ins", usage = "Number of insertion SV to sample ["+NUM_INS_ARG+"]")
    int num_INS = NUM_INS_ARG;

    static final int NUM_DEL_ARG = 2000;
    @Option(name = "-num_del", usage = "Number of deletion SV to sample ["+NUM_DEL_ARG+"]")
    int num_DEL = NUM_DEL_ARG;

    static final int NUM_DUP_ARG = 500;
    @Option(name = "-num_dup", usage = "Number of duplications to sample ["+NUM_DUP_ARG+"]")
    int num_DUP = NUM_DUP_ARG;

    static final int NUM_INV_ARG = 500;
    @Option(name = "-num_inv", usage = "Number of inversions to sample ["+NUM_INV_ARG+"]")
    int num_INV = NUM_INV_ARG;

    static final double NOVEL_RATIO_ARG = 0.01;
    @Option(name = "-novel", usage = "Average ratio of novel variants["+NOVEL_RATIO_ARG+"]")
    double ratio_novel = NOVEL_RATIO_ARG;

    static final int MIN_LEN_ARG = 100;
    @Option(name = "-min_len", usage = "Minimum variant length ["+MIN_LEN_ARG+"], inclusive")
    int min_length_lim = MIN_LEN_ARG;

    static final int MAX_LEN_ARG = 1000000;
    @Option(name = "-max_len", usage = "Maximum variant length ["+MAX_LEN_ARG+"], inclusive")
    int max_length_lim = MAX_LEN_ARG;

    @Option(name = "-ref", usage = "Reference Genome [Required]",metaVar = "file",required = true)
    String reference_filename;

    @Option(name = "-ins", usage = "Known Insertion Sequences [Required]",metaVar = "file",required = true)
    String insert_filename;

    @Option(name = "-dgv", usage = "DGV database flat file [Required]",metaVar = "file",required = true)
    String dgv_filename;

    int num_novel_added = 0;

    RandDGV2VCF() {
        super();
        num_novel_added = 0;
    }

    RandDGV2VCF(long seed) {
        super(seed);
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
        double rand_num = _rand.nextDouble();
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
            int new_pos = _rand.nextInt(end_val - start_val + 1) + start_val + 1;
            while (!var.setNovelPosition(new_pos, ref)) {
                if (time_out > 100) {
                    log.warn("Error, cannot set novel position: " + (end_val - start_val + 1));
                    log.warn(var.deletion());
                    log.warn(var);
                    //System.exit(1);
                }

                log.info(time_out + " : " + new_pos + " : " + var.deletion());

                new_pos = _rand.nextInt(end_val - start_val + 1) + start_val + 1;
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
        String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
        String usage = "Outputs VCF to stdout. Randomly samples variants from DGV flat file.\n";

        CmdLineParser parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java -jar randdgv2vcf.jar [options...]");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }

        if (ratio_novel > 1 || ratio_novel < 0) {
            System.err.println("Novel ratio out of range [0,1]");
            System.exit(1);
        }

        _rand = new Random(seed);

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
        DGVparser parser_one = new DGVparser(dgv_filename, ref,_rand);
        Variant prev_var = new Variant(_rand);

        // Read through a first time to generate the counts for sampling without replacement
        while (parser_one.hasMoreInput()) {
            Variant var = parser_one.parseLine();
            if (var == null) {
                continue;
            }

            // select genotypes here
            int chr_idx = var.chromosome();
            int num_alt = var.get_num_alt();

            Genotypes geno = new Genotypes(chr_idx, num_alt, _rand);
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
        parser_one = new DGVparser(dgv_filename, ref,_rand);
        prev_var = new Variant(_rand);


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

            if (var.max_len() > max_length_lim || var.min_len() < min_length_lim) {
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
                                num_INS, total_num_INS,output_all);
                        break;
                    case Deletion:
                        geno.geno[i] = sample_genotype(geno.geno[i], DEL_params,
                                num_DEL, total_num_DEL,output_all);
                        break;
                    case Tandem_Duplication:
                        geno.geno[i] = sample_genotype(geno.geno[i], DUP_params,
                                num_DUP, total_num_DUP,output_all);
                        break;
                    case Inversion:
                        geno.geno[i] = sample_genotype(geno.geno[i], INV_params,
                                num_INV, total_num_INV,output_all);
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
