package com.binatechnologies.varsim;

import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import java.io.*;
import java.util.Random;

/**
 *
 */

/**
 * @author johnmu
 */

public class RandBED2VCF extends randVCFgenerator {
    private final static Logger log = Logger.getLogger(RandBED2VCF.class.getName());

    SimpleReference ref;
    int num_novel_added = 0;
    int var_idx = 0;
    byte[] insert_seq = null;

    // parameters
    static final int MIN_LEN_ARG = 50;
    @Option(name = "-min_len", usage = "Minimum variant length ["+MIN_LEN_ARG+"], inclusive")
    int min_length_lim = MIN_LEN_ARG;

    static final int MAX_LEN_ARG = 1000000;
    @Option(name = "-max_len", usage = "Maximum variant length ["+MAX_LEN_ARG+"], inclusive")
    int max_length_lim = MAX_LEN_ARG;

    static final int SEED_ARG = 333;
    @Option(name = "-seed", usage = "Seed for random sampling ["+SEED_ARG+"]")
    int seed = 333;

    @Option(name = "-ref", usage = "Reference Genome",metaVar = "file",required = true)
    String reference_filename;

    @Option(name = "-ins", usage = "Known Insertion Sequences",metaVar = "file",required = true)
    String insert_filename;

    @Option(name = "-ins_bed", usage = "Known Insertion BED file",metaVar = "BED_file",required = true)
    String ins_bed_filename;

    @Option(name = "-del_bed", usage = "Known Deletion BED file",metaVar = "BED_file",required = true)
    String del_bed_filename;

    RandBED2VCF() {
        super();
        num_novel_added = 0;
        var_idx = 0;
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        RandBED2VCF runner = new RandBED2VCF();
        runner.run(args);
    }


    // outputs VCF record with random phase
    void rand_output_vcf_record(BufferedWriter bw, Variant var) throws IOException {

        // this is ok if both the same genotype
        // the second call will return
        fill_in_seq(var, insert_seq, var.getGeno().geno[0]);
        fill_in_seq(var, insert_seq, var.getGeno().geno[1]);

        output_vcf_record(bw, var);

    }

    // remember BED is 0-based
    Variant parse_bed_line(String line, Variant.Type type) {
        String[] ll = line.split("\t");
        if (ll.length < 4) {
            return new Variant();
        }

//        String chr_name, int chr, int pos, int del, byte[] ref,
//        FlexSeq[] alts, byte[] phase, String var_id, String filter,
//                String ref_deleted

        int chr_idx = variantFileParser.getChromIndex(ll[0]);
        int pos = Integer.parseInt(ll[1]);
        int end = Integer.parseInt(ll[2]);
        int len = Integer.parseInt(ll[3]);

        if (len > max_length_lim || len < min_length_lim) {
            return null;
        }

        FlexSeq[] alts = new FlexSeq[1];
        String var_idx_str = "";
        byte[] ref_seq;
        if (type == Variant.Type.Deletion) {
            alts[0] = new FlexSeq();
            var_idx_str = "del_";
            ref_seq = ref.byteRange(chr_idx, pos, end);
        } else if (type == Variant.Type.Insertion) {
            alts[0] = new FlexSeq(FlexSeq.Type.INS, len);
            var_idx_str = "ins_";
            ref_seq = new byte[0];
        } else {
            log.error("Bad type!");
            return null;
        }
        var_idx_str += var_idx;
        var_idx++;

        Genotypes geno = new Genotypes(chr_idx, 1, rand);

        return new Variant(ll[0], chr_idx, pos, ref_seq.length, ref_seq, alts,
                geno.geno, false, var_idx_str, "PASS", "");

    }

    private void process_bed(BufferedWriter out, BufferedReader bed_reader, Variant.Type type)
            throws IOException {

        String line;
        while ((line = bed_reader.readLine()) != null) {
            Variant var = parse_bed_line(line, type);
            if (var == null) {
                // System.err.println("Bad variant or not a variant line");
                continue;
            }

            // write out variant
            try {
                if (!var.isRef()) {
                    rand_output_vcf_record(out, var);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }


        }

    }

    public void run(String[] args) {
        String usage = "Generates a VCF file (to stdout) from an insertion and a deletion BED file. Insertions sequences\n"
                + "are randomly sampled from the insert_seq file. This is designed for the Venter SV BED files. \n";

        CmdLineParser parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("java -jar randbed2vcf [options...] arguments...");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println(usage);
            System.err.println();
            return;
        }

        if (max_length_lim < min_length_lim) {
            log.error("Bad lengths");
        }

        rand = new Random(seed);

        log.info("Reading reference");
        ref = new SimpleReference(reference_filename);

        // read in the insert sequences
        log.info("Reading insert sequences");

        insert_seq = null;

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

        BufferedReader ins_bed_reader = null;
        BufferedReader del_bed_reader = null;

        try {
            ins_bed_reader = new BufferedReader(new FileReader(ins_bed_filename));
            del_bed_reader = new BufferedReader(new FileReader(del_bed_filename));
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        // read through DGV file again and output the new sampled VCF file
        log.info("Writing VCF file");

        // write the header
        System.out.print("##fileformat=VCFv4.0\n");
        System.out
                .print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        System.out
                .print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVenter\n");

        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new OutputStreamWriter(System.out));
        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            process_bed(out, del_bed_reader, Variant.Type.Deletion);
            process_bed(out, ins_bed_reader, Variant.Type.Insertion);

            out.flush();
            out.close();
        } catch (IOException e) {
            log.error(e);
        }

    }

}
