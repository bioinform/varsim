package com.bina.varsim.tools.simulation;

import com.bina.varsim.VarSimToolNamespace;
import com.bina.varsim.constants.Constant;
import com.bina.varsim.types.*;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.types.variant.VariantType;
import com.bina.varsim.types.variant.alt.Alt;
import com.bina.varsim.util.SimpleReference;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.Random;

/**
 *
 */

/**
 * @author johnmu
 */

public class RandBED2VCF extends RandVCFgenerator {
    // parameters
    static final int MIN_LEN_ARG = Constant.SVLEN;
    static final int MAX_LEN_ARG = 1000000;
    private final static Logger log = Logger.getLogger(RandBED2VCF.class.getName());

    SimpleReference ref;
    int num_novel_added = 0;
    int var_idx = 0;
    byte[] insert_seq = null;
    @Option(name = "-min_len", usage = "Minimum variant length [" + MIN_LEN_ARG + "], inclusive")
    int min_length_lim = MIN_LEN_ARG;
    @Option(name = "-max_len", usage = "Maximum variant length [" + MAX_LEN_ARG + "], inclusive")
    int max_length_lim = MAX_LEN_ARG;
    @Option(name = "-ref", usage = "Reference Genome [Required]", metaVar = "file", required = true)
    String reference_filename;

    @Option(name = "-t", usage = "Gender of individual [MALE]")
    GenderType gender = GenderType.MALE;

    @Option(name = "-ins", usage = "Known Insertion Sequences [Required]", metaVar = "file", required = true)
    String insert_filename;

    @Option(name = "-ins_bed", usage = "Known Insertion BED file [Required]", metaVar = "BED_file", required = true)
    String ins_bed_filename;

    @Option(name = "-del_bed", usage = "Known Deletion BED file [Required]", metaVar = "BED_file", required = true)
    String del_bed_filename;

    @Option(name = "-dup_bed", usage = "Known Duplication BED file [Required]", metaVar = "BED_file", required = true)
    String dup_bed_filename;

    @Option(name = "-inv_bed", usage = "Known Inversions BED file [Required]", metaVar = "BED_file", required = true)
    String inv_bed_filename;

    public RandBED2VCF(final String command, final String description) {
        super(command, description);
        num_novel_added = 0;
        var_idx = 0;
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        new RandBED2VCF("", VarSimToolNamespace.RandBED2VCF.description).run(args);
    }


    // outputs VCF record with random phase
    void randOutputVcfRecord(BufferedWriter bw, Variant var) throws IOException {

        // this is ok if both the same genotype
        // the second call will return
        fillInSeq(var, insert_seq, var.getGenotypes().geno[0]);
        fillInSeq(var, insert_seq, var.getGenotypes().geno[1]);

        outputVcfRecord(bw, var);

    }

    // Parses a BED line into a Variant object, based on the input line and
    // VariantType
    private String[] splitLine(String line) {
        line = line.trim();
        return line.split("\t");
    }

    // Checks whether the first element of an array of strings is a comment (starts
    // with '#')
    private boolean isComment(String[] ll) {
        return ll[0].charAt(0) == '#';
    }

    // Constructs a new ChrString object with the given input string
    private ChrString getChrString(String chrString) {
        return new ChrString(chrString);
    }

    // Parses an integer from the input string and adds 1 to it
    private int getPos(String posString) {
        return Integer.parseInt(posString) + 1;
    }

    // Gets the byte array of the insertion sequence from the given array of strings
    private byte[] getInsertionSequence(String[] meta) {
        try {
            return meta[1].getBytes("US-ASCII");
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
            return null;
        }
    }

    // Creates an array of Alt objects based on the given VariantType, length, and
    // insertion sequence
    private Alt[] createAlts(VariantType type, int len, byte[] ins_seq) {
        Alt[] alts = new Alt[1];
        if (type == VariantType.Deletion) {
            alts[0] = new Alt(new FlexSeq());
        } else if (type == VariantType.Insertion) {
            if (ins_seq != null) {
                alts[0] = new Alt(new FlexSeq(ins_seq));
            } else {
                alts[0] = new Alt(new FlexSeq(FlexSeq.Type.INS, len));
            }
        } else if (type == VariantType.Tandem_Duplication) {
            alts[0] = new Alt(new FlexSeq(FlexSeq.Type.TANDEM_DUP, len, 2));
        } else if (type == VariantType.Inversion) {
            alts[0] = new Alt(new FlexSeq(FlexSeq.Type.INV, len));
        } else {
            log.error("Bad type!");
            return null;
        }
        return alts;
    }

    // Gets the byte array of the reference sequence based on the given ChrString,
    // position, length, and VariantType
    private byte[] getReferenceSequence(ChrString chr, int pos, int len, VariantType type) {
        if (type == VariantType.Deletion) {
            return ref.byteRange(chr, pos, pos + len);
        } else {
            return new byte[0];
        }
    }

    // Generates a unique identifier for the Variant based on the given VariantType
    private String getVarId(VariantType type) {
        String var_idx_str;
        if (type == VariantType.Deletion) {
            var_idx_str = "del_";
        } else if (type == VariantType.Insertion) {
            var_idx_str = "ins_";
        } else if (type == VariantType.Tandem_Duplication) {
            var_idx_str = "dup_";
        } else if (type == VariantType.Inversion) {
            var_idx_str = "inv_";
        } else {
            log.error("Bad type!");
            return null;
        }
        var_idx_str += var_idx;
        var_idx++;
        return var_idx_str;
    }

    // Parses a BED line into a Variant object, based on the input line and
    // VariantType
    Variant parseBedLine(String line, VariantType type) {
        String[] ll = splitLine(line);
        if (ll.length < 4)
            return new Variant(rand);

        if (isComment(ll)) {
            return null;
        }

        ChrString chr = getChrString(ll[0]);
        int pos = getPos(ll[1]);
        String[] meta = ll[3].split(",");
        byte[] ins_seq = null;
        int len = 1;
        if (meta[0].equals("seq")) {
            ins_seq = getInsertionSequence(meta);
        } else {
            len = Integer.parseInt(meta[0]);
        }

        if (pos == 0 || len > max_length_lim || len < min_length_lim)
            return null;

        Alt[] alts = createAlts(type, len, ins_seq);

        if (alts == null) {
            return null;
        }

        byte[] ref_seq = getReferenceSequence(chr, pos, len, type);

        if (ref_seq == null && type == VariantType.Deletion) {
            log.error("Range error: " + line);
            return null;
        }

        String var_idx_str = getVarId(type);

        if (var_idx_str == null) {
            return null;
        }

        Genotypes geno = new Genotypes(chr, gender, 1, rand);

        /*return new Variant(chr, pos, ref_seq.length, ref_seq, alts,
                geno.geno, false, var_idx_str, "PASS", String.valueOf(ref
                .charAt(chr, pos - 1)), _rand);*/
        return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(ref_seq.length)
                .ref(ref_seq).alts(alts).phase(geno.geno).isPhased(false).varId(var_idx_str)
                .filter("PASS").refDeleted(String.valueOf(ref.charAt(chr, pos - 1))).randomNumberGenerator(rand)
                .build();
    }

    private void process_bed(BufferedWriter out, BufferedReader bed_reader, VariantType type)
            throws IOException {

        String line;
        while ((line = bed_reader.readLine()) != null) {
            Variant var = parseBedLine(line, type);
            if (var == null) {
                log.error("Bad variant or not a variant line: " + line);
                continue;
            }

            // write out variant
            try {
                if (!var.isRef()) {
                    randOutputVcfRecord(out, var);
                } else {
                    //log.error("Reference variant: " + line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }


        }

    }

    public void run(String[] args) {
        if (!parseArguments(args)) {
            return;
        }

        if (max_length_lim < min_length_lim) {
            log.error("Bad lengths, max < min");
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
            String line;
            while ((line = bufferedReader.readLine()) != null) {
                line = line.trim();
                sb.append(line);
            }
            bufferedReader.close();
            insert_seq = sb.toString().getBytes("US-ASCII");
        } catch (IOException e) {
            e.printStackTrace();
        }

        BufferedReader ins_bed_reader;
        BufferedReader del_bed_reader;
        BufferedReader dup_bed_reader;
        BufferedReader inv_bed_reader;

        // try to open each of them to make sure the file can open
        try {
            ins_bed_reader = new BufferedReader(new FileReader(ins_bed_filename));
            del_bed_reader = new BufferedReader(new FileReader(del_bed_filename));
            dup_bed_reader = new BufferedReader(new FileReader(dup_bed_filename));
            inv_bed_reader = new BufferedReader(new FileReader(inv_bed_filename));
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
            return;
        }

        try {
            process_bed(out, del_bed_reader, VariantType.Deletion);
            process_bed(out, ins_bed_reader, VariantType.Insertion);
            process_bed(out, dup_bed_reader, VariantType.Tandem_Duplication);
            process_bed(out, inv_bed_reader, VariantType.Inversion);

            out.flush();
            out.close();
        } catch (IOException e) {
            log.error(e);
        }

    }

}
