package com.bina.varsim.tools.simulation;

import com.bina.varsim.constants.Constant;
import com.bina.varsim.types.*;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.util.SimpleReference;
import com.bina.varsim.util.VCFparser;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Random;
/**
 * TODO ignores the input genotypes for now
 */

/**
 * @author johnmu
 */


public class RandVCF2VCF extends randVCFgenerator {
    static final int SEED_ARG = 333;
    static final int NUM_SNP_ARG = 3000000;
    static final int NUM_INS_ARG = 100000;
    static final int NUM_DEL_ARG = 100000;
    static final int NUM_MNP_ARG = 20000;
    static final int NUM_COMPLEX_ARG = 20000;
    static final double NOVEL_RATIO_ARG = 0.01;
    static final int MIN_LEN_ARG = 1;
    static final int MAX_LEN_ARG = Constant.SVLEN - 1;
    static final double PROP_HET_ARG = 0.6;
    private final static Logger log = Logger.getLogger(RandVCF2VCF.class.getName());
    @Option(name = "-seed", usage = "Seed for random sampling [" + SEED_ARG + "]")
    int seed = SEED_ARG;
    @Option(name = "-num_snp", usage = "Number of SNPs to sample [" + NUM_SNP_ARG + "]")
    int num_SNP = NUM_SNP_ARG;
    @Option(name = "-num_ins", usage = "Number of simple insertions to sample [" + NUM_INS_ARG + "]")
    int num_INS = NUM_INS_ARG;
    @Option(name = "-num_del", usage = "Number of simple deletions to sample [" + NUM_DEL_ARG + "]")
    int num_DEL = NUM_DEL_ARG;
    @Option(name = "-num_mnp", usage = "Number of MNPs to sample [" + NUM_MNP_ARG + "]")
    int num_MNP = NUM_MNP_ARG;
    @Option(name = "-num_complex", usage = "Number of complex variants (other ones) to sample [" + NUM_COMPLEX_ARG + "]")
    int num_COMPLEX = NUM_COMPLEX_ARG;
    @Option(name = "-novel", usage = "Average ratio of novel variants[" + NOVEL_RATIO_ARG + "]")
    double ratio_novel = NOVEL_RATIO_ARG;
    @Option(name = "-min_len", usage = "Minimum variant length [" + MIN_LEN_ARG + "], inclusive")
    int min_length_lim = MIN_LEN_ARG;
    @Option(name = "-max_len", usage = "Maximum variant length [" + MAX_LEN_ARG + "], inclusive")
    int max_length_lim = MAX_LEN_ARG;
    @Option(name = "-prop_het", usage = "Average ratio of novel variants[" + PROP_HET_ARG + "]")
    double prop_het = PROP_HET_ARG;

    @Option(name = "-ref", usage = "Reference Genome [Required]", metaVar = "file", required = true)
    String reference_filename;

    @Option(name = "-vcf", usage = "Known VCF file, eg. dbSNP [Required]", metaVar = "file", required = true)
    String vcf_filename;

    @Option(name = "-t", usage = "Gender of individual [MALE]")
    GenderType gender = GenderType.MALE;


    int num_novel_added;

    public RandVCF2VCF() {
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

        ChrString chr = var.getChr();

        // determine whether this one is novel
        double rand_num = _rand.nextDouble();
        if (rand_num <= ratio_novel) {
            // make the variant novel, simply modify it
            // TODO maybe modifying it is bad

            int chr_len = ref.getRefLen(chr);
            int buffer = Math.max(
                    10,
                    Math.max(var.maxLen(geno.geno[0]),
                            var.maxLen(geno.geno[1])));
            int start_val = Math.min(buffer, Math.max(chr_len - buffer, 0));
            int end_val = Math.max(chr_len - buffer, Math.min(buffer, chr_len));

            int time_out = 0;
            int rand_pos = _rand.nextInt(end_val - start_val + 1) + start_val + 1;
            while (!var.setNovelPosition(rand_pos, ref)) {
                if (time_out > 100) {
                    log.warn("Error: cannot set novel position: " + (end_val - start_val + 1));
                    log.warn(var.deletion());
                    log.warn(var);
                    break;
                }

                rand_pos = _rand.nextInt(end_val - start_val + 1) + start_val + 1;
                //log.info(time_out + " : " + var.deletion());

                time_out++;
            }
            num_novel_added++;

            var.setVarID("Novel_" + num_novel_added);
        }

        output_vcf_record(bw, var, geno.geno[0], geno.geno[1]);

    }

    public void run(String[] args) {
        String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
        String usage = "Outputs VCF to stdout. Randomly samples variants from VCF file.";

        CmdLineParser parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java -jar randvcf2vcf.jar [options...]");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }

        if (ratio_novel > 1 || ratio_novel < 0) {
            System.err.println("Novel ratio out of range [0,1]");
            System.exit(1);
        }

        SimpleReference ref = new SimpleReference(reference_filename);

        _rand = new Random(seed);

        // read through VCF file and count the
        log.info("Counting variants and assigning genotypes");
        ArrayList<Genotypes> selected_geno = new ArrayList<>();

        int total_num_SNP = 0;
        int total_num_INS = 0;
        int total_num_DEL = 0;
        int total_num_MNP = 0;
        int total_num_COMPLEX = 0;
        int total_num_other = 0;
        int total_num = 0;
        VCFparser parser_one = new VCFparser(vcf_filename, null, false, _rand);
        Variant prev_var = new Variant(_rand);

        // read though once to count the totals, this is so we don't have
        // to store an array of the variants for sampling without replacement
        while (parser_one.hasMoreInput()) {
            Variant var = parser_one.parseLine();
            if (var == null) {
                continue;
            }

            // select genotypes here
            ChrString chr = var.getChr();
            int num_alt = var.get_num_alt();

            Genotypes geno = new Genotypes(chr, gender, num_alt, _rand, prop_het);
            selected_geno.add(geno);

            if (prev_var.equals(var)) {
                continue;
            }

            // this is ok because var is not changed later
            prev_var = var;

            if (var.maxLen() > max_length_lim
                    || var.minLen() < min_length_lim) {
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
                        //log.error(i + " : " + geno.geno[i] + " : " + var.getType(geno.geno[i]) + " : OTHER : " + var);
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

        parser_one = new VCFparser(vcf_filename, null, false, _rand);
        prev_var = new Variant(_rand);

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

            if (var.maxLen() > max_length_lim
                    || var.minLen() < min_length_lim) {
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
