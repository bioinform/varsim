package com.binatechnologies.varsim;

import org.apache.log4j.Logger;

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
    int min_length_lim = 50;
    int max_length_lim = 1000000;

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
        String usage = "RandBED2VCF seed min_length_lim max_length_lim reference_file insert_seq.txt " +
                "ins_file.bed del_file.bed\n"
                + "     seed            -- Seed for random number generator (for random genotypes)\n"
                + "     min_length_lim  -- Minimum length variant to generate (inclusive)\n"
                + "     max_length_lim  -- Maximum length variant to generate (inclusive)\n"
                + "     reference_file  -- Reference genome sequence, eg. b37\n"
                + "     insert_seq      -- Single line concatenation of known insertion sequences\n"
                + "     ins_file        -- BED file of insertions\n"
                + "     del_file        -- BED file of deletions\n"
                + "Generates a VCF file (to stdout) from an insertion and a deletion BED file. Insertions sequences\n"
                + "are randomly sampled from the insert_seq file. This is designed for the Venter SV BED files. \n";
        if (args.length != 7) {
            System.err.println(usage);
            System.exit(1);
        }

        int seed = Integer.parseInt(args[0]);
        min_length_lim = Integer.parseInt(args[1]);
        max_length_lim = Integer.parseInt(args[2]);
        String reference_filename = args[3];
        String insert_filename = args[4];
        String ins_bed_filename = args[5];
        String del_bed_filename = args[6];

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
