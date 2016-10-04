package com.bina.varsim.util;

//--- Java imports ---

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.variant.Variant;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.IllegalFormatException;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VCFparser extends GzFileParser<Variant> {
    private final static Logger log = Logger.getLogger(VCFparser.class.getName());

    private Random _rand = null;

    private int _id_ind = -1;
    private String _id = null;
    private boolean _pass = false;
    private boolean chrom_exists = false;

    private VCFparser() {
        _id_ind = 10; // the first sample
    }

    /**
     * Reads a VCF file line by line
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param id       ID of individual, essentially selects a column, null to use first ID column
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, String id, boolean pass, Random rand) {
        _rand = rand;
        try {
            _br = new BufferedReader(new InputStreamReader(decompressStream(fileName)));
            readLine();
        } catch (Exception ex) {
            log.error("Can't open file " + fileName);
            log.error(ex.toString());
        }
        _id = id;
        _pass = pass;

        if (_id == null) {
            _id_ind = 10; // the first sample
        }
    }

    /**
     * Reads a VCF file line by line
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param id       ID of individual, essentially selects a column, null to use first ID column
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, String id, boolean pass) {
        this(fileName, id, pass, null);
    }

    /**
     * Reads a VCF file line by line, if there are multiple individuals, takes the first one
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, boolean pass, Random rand) {
        this(fileName, null, pass, rand);
    }

    /**
     * Reads a VCF file line by line, if there are multiple individuals, takes the first one
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, boolean pass) {
        this(fileName, null, pass, null);
    }

    public static void main(String args[]) {
        VCFparser runner = new VCFparser();
        Variant v = runner.process_line("12\t29557989\t.\tACAAAAGAAATGATCATGTTTGTAGGT\tAAAAAGAAATGATCATGTTTGTAGGT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        System.err.println(v);
        v = runner.process_line("12\t29557989\t.\tACTTT\tACGTTTT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        System.err.println(v);
        v = runner.process_line("12\t29557989\t.\tACT\tAAAACT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        System.err.println(v);
        v = runner.process_line("15\t85825565\tnssv534459\tT\t<DUP:TANDEM>\t.\tPASS\tSVTYPE=DUP;SVLEN=284016\tGT:CN\t0|1:2|2");
        System.err.println(v); // TODO this one fails for now
    }

    /**
     * gets the SVLEN value
     *
     * @param info string from the INFO field of a VCF file
     * @return array of SVLENs that are read, returns empty array if there are none
     */
    int[] getSVLen(String info) {
        Pattern svLen_pattern = Pattern.compile("(.*)SVLEN=([-]?[0-9,]+)(.*)");
        Matcher matcher = svLen_pattern.matcher(info);
        int sv_len[] = {};
        if (matcher.find()) {
            String[] sv_lenString = matcher.group(2).split(",");
            sv_len = new int[sv_lenString.length];
            for (int i = 0; i < sv_lenString.length; i++) {
                sv_len[i] = Integer.parseInt(sv_lenString[i]);
            }
        }
        return sv_len;
    }

    /**
     * gets the END value
     *
     * @param info string from the INFO field of a VCF file
     * @return Returns the END value from INFO field, returns 0 if there is none
     */
    int getEndInfo(String info) {
        Pattern endLoc_pattern = Pattern.compile("(.*)END=([0-9]+)(.*)");
        Matcher matcher = endLoc_pattern.matcher(info);
        if (matcher.find())
            return Integer.parseInt(matcher.group(2));
        else
            return 0; // No end information found
    }

    /**
     * finds where "GT" or similar is in the VCF string so that the genotype can be read
     *
     * @param record The column in VCF that contains GT, CN, etc..
     * @param key    The key to be found "GT" or "CN" or other
     * @return the index of the key
     */
    private int getFormatKeyIndex(String record, String key) {
        StringTokenizer words = new StringTokenizer(record, ":");
        int ret = 0;
        while (words.hasMoreTokens()) {
            if (words.nextToken().equals(key))
                return ret;
            else
                ret++;
        }
        return -1;
    }

    /**
     * Takes genotype string and splits it into alleles, supports a maximum of two
     *
     * @param geno genotype string corresponding to the GT tag (field index 9)
     * @param vals Return the genotypes read here [paternal,maternal]
     * @param chr  chromosome we are dealing with, some are haploid, need to assign parent
     * @return true if the variant is phased
     */
    boolean splitGeno(String geno, byte[] vals, ChrString chr) {
        boolean is_phased = false;

        geno = geno.trim();
        boolean strangePhase = false;
        if (geno.matches("^[0-9]+$")) {
            // phase is only a single number, for haploid chromosomes
            byte val = (byte) Integer.parseInt(geno);
            if (chr.isX()) {
                vals[1] = val; // maternal
                is_phased = true;
            } else if (chr.isY()) {
                vals[0] = val; // paternal
                is_phased = true;
            } else if (chr.isMT()) {
                vals[1] = val;
                is_phased = true;
            } else {
                vals[0] = vals[1] = val;
            }
        } else if (geno.length() >= 3) {
                    // this is the case where phase looks like "1|0" or "10|4"
                    String[] ll = geno.split("[\\|/]");
                    int c1 = -1;
                    int c2 = -1;
                    char phasing = '/';
                    if (ll.length == 2) {
                        try {
                            c1 = Integer.parseInt(ll[0]);
                            c2 = Integer.parseInt(ll[1]);
                            phasing = geno.charAt(ll[0].length());
                        } catch (NumberFormatException e) {
                            strangePhase = true;
                }
            } else {
                strangePhase = true;
            }

            if (c1 >= 0 && c2 >= 0) {
                vals[0] = (byte) c1;
                vals[1] = (byte) c2;
                if (phasing == '|') {
                    is_phased = true;
                }
            } else {
                strangePhase = true;
            }
        } else {
            strangePhase = true;
        }

        if (strangePhase) {
            // System.err.println("Unrecognized phasing '" + phase + "'.");
            vals[0] = -1;
            vals[1] = -1;
            is_phased = false;
        }

        return is_phased;
    }

    /**
     * takes a line from a VCF file, parse it,
     * return a Variant object
     * @param line
     * @return
     */
    public Variant process_line(String line) {

        // try to determine the column we should read for the genotype
        StringTokenizer toks = new StringTokenizer(line);
        if (line.startsWith("#")) {
            if (_id != null && line.startsWith("#CHROM")) {
                chrom_exists = true;
                int index = 0;
                while (toks.hasMoreTokens()) {
                    index++;
                    String tok = toks.nextToken();
                    if (tok.equals(_id))
                        _id_ind = index;
                }
            } else if (_id == null) {
                _id_ind = 10; // the first sample
            }
            return null;
        }

        // If we cannot determine, then use the first one
        if (_id_ind < 0 && !chrom_exists) {
            _id_ind = 10;
        } else if (_id_ind < 0) {
            _id_ind = 10;
            log.warn("Warning!!! ID (" + _id + ") does not exist... ");
        }


        int index = 0, genotype_ind = -1, copynum_ind = -1;
        int pos = -1;
        ChrString chr = null;
        String REF = "", FILTER = "", ALT = "", var_id = "";
        String phase = ".", copy_num = "0/0", INFO = "", FORMAT;
        String[] sampleInfo;
        while (toks.hasMoreTokens()) {
            index++;
            if (index == 1) { // Parsing chromosome
                chr = new ChrString(toks.nextToken());
            } else if (index == 2) // Parsing position
                pos = Integer.parseInt(toks.nextToken());
            else if (index == 3) // Parsing position
                var_id = toks.nextToken();
            else if (index == 4) // Parsing reference allele
                REF = toks.nextToken();
            else if (index == 5) // Parsing alternative allele
                ALT = toks.nextToken();
            else if (index == 7) // FILTER field
                FILTER = toks.nextToken();
            else if (index == 8) // INFO field
                INFO = toks.nextToken();
            else if (index == 9) { // Output format
                FORMAT = toks.nextToken();
                genotype_ind = getFormatKeyIndex(FORMAT, "GT");
                copynum_ind = getFormatKeyIndex(FORMAT, "CN");
            } else if (index == _id_ind) { // Phasing
                sampleInfo = (toks.nextToken()).split(":");
                if (genotype_ind >= 0) {
                    phase = sampleInfo[genotype_ind];
                }
                if (copynum_ind >= 0) {
                    copy_num = sampleInfo[copynum_ind];
                }

                break;
            } else {
                toks.nextToken();
            }
        }

        // unknown chromosome
        if (chr == null) {
            log.warn("Bad chromosome name: " + line);
            return null;
        }

        if (_pass && !(FILTER.contains("PASS") || FILTER.equals("."))) {
            //log.warn("not pass line" + line);
            return null; // Filtered out
        }
        // parse the phase
        byte[] phase_val = new byte[2]; // paternal-maternal
        boolean is_phased = splitGeno(phase, phase_val, chr);


        if (genotype_ind >= 0 && phase_val[0] == 0 && phase_val[1] == 0) {
            //return null; // reference alleles... ignore them for now....
        }

        // determine copy-number
        // TODO need to be able to deal with unphased copy-numbers?
        byte[] copy_num_val = new byte[2]; // paternal-maternal
        boolean is_cn_phased;

        if (copynum_ind >= 0) {
            is_cn_phased = splitGeno(copy_num, copy_num_val, chr);
            if (is_cn_phased != is_phased) {
                // TODO maybe don't throw error, this is not standard format
                // anyways
                log.error("Inconsistent copy number:");
                log.error(line);
                return null;
            }
        }

        // Upper casing
        REF = REF.toUpperCase();
        ALT = ALT.toUpperCase();

        String ref_deleted = "";
        FlexSeq alts[];

        /*if symbolic alleles are present,
        make sure # of alleles equal # of
        SV lengths
         */
        /*!!!!!!!!!!
        CAUTION: we assume symbolic alleles are not mixed
        with non-symbolic alleles.
         */
        if (ALT.indexOf('<') != -1) {
            String[] alts_str = ALT.split(",");
            int[] inv_lens = getSVLen(INFO);
            if (alts_str.length != inv_lens.length) {
                throw new IllegalArgumentException("ERROR: number of symbolic alleles is unequal to number of SV lengths.");
            }
            for (int i = 0; i < alts_str.length; i++) {
                if (!alts_str[i].startsWith("<")) {
                    throw new IllegalArgumentException("ERROR: symbolic alleles are mixed with non-symbolic alleles");
                }
            }
        }
        if (ALT.startsWith("<INV>")) {
            // inversion SV

            int inv_lens[] = getSVLen(INFO);
            int end_loc = getEndInfo(INFO); // TODO may need to check

            ref_deleted = REF;
            byte[] refs = new byte[0];
            pos++;

            if (inv_lens.length > 0) {
                alts = new FlexSeq[inv_lens.length];
                for (int i = 0; i < inv_lens.length; i++) {
                    int len_val = Math.max(Math.abs(inv_lens[i]), 1);
                    alts[i] = new FlexSeq(FlexSeq.Type.INV, len_val);
                }
                // TODO this assumes only one alt
                return new Variant(chr, pos, Math.abs(inv_lens[0]), refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else if (end_loc > 0) {
                int inv_len = Math.max(Math.abs(end_loc - pos + 1), 1);
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.INV, inv_len);
                return new Variant(chr, pos, inv_len, refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else {
                log.error("No length information for INV:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.startsWith("<DUP>") || ALT.startsWith("<DUP:TANDEM>")) {
            // duplication or tandem duplication SV
            int dup_lens[] = getSVLen(INFO);
            int end_loc = getEndInfo(INFO); // may need to check inconsistency
            // btw length and end location

            ref_deleted = REF;
            byte[] refs = new byte[0];
            pos++;

            if (dup_lens.length > 0) {
                alts = new FlexSeq[dup_lens.length];
                for (int i = 0; i < dup_lens.length; i++) {
                    // TODO this is temporary, how to encode copy number?
                    int copy_val = 1;
                    for (int j = 0; j < 2; j++) {
                        if ((i + 1) == phase_val[j]) {
                            /*
                            if i = 0, phase_val[0] = 1, phase_val[1] = 1
                            copy_num_val[0] = 3, copy_num_val[1] = 2
                            then copy_val = 2.
                            what does copy_val mean in real world?
                             */
                            if (copy_num_val[j] > 0) {
                                copy_val = copy_num_val[j];
                            }
                        }
                    }

                    int len_val = Math.max(Math.abs(dup_lens[i]), 1);

                    alts[i] = new FlexSeq(FlexSeq.Type.DUP, len_val, copy_val);
            }

                return new Variant(chr, pos, Math.abs(dup_lens[0]), refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else if (end_loc > 0) {
                int dup_len = Math.max(Math.abs(end_loc - pos + 1), 1);
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.DUP, dup_len, Math.max(
                        copy_num_val[0], copy_num_val[1]));

                return new Variant(chr, pos, dup_len, refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else {
                log.error("No length information for DUP:");
                log.error(line);
                log.error("skipping...");
                return null;
            }

        } else if (ALT.startsWith("<INS>")) {
            // insertion SV

            int ins_lens[] = getSVLen(INFO);
            int end_loc = getEndInfo(INFO); // TODO may need to check

            ref_deleted = REF;
            byte[] refs = new byte[0];
            pos++;

            if (ins_lens.length > 0) {
                alts = new FlexSeq[ins_lens.length];
                for (int i = 0; i < ins_lens.length; i++) {
                    int len_val;
                    if (ins_lens[i] == 0) {
                        len_val = Integer.MAX_VALUE;
                    } else {
                        len_val = Math.max(Math.abs(ins_lens[i]), 1);
                    }
                    alts[i] = new FlexSeq(FlexSeq.Type.INS, len_val);
                }
                return new Variant(chr, pos, 0, refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else if (end_loc > 0) {
                int ins_len = Math.max(Math.abs(end_loc - pos), 1);
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.INS, ins_len);
                return new Variant(chr, pos, 0, refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else {
                log.error("No length information for INS:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.startsWith("<DEL>")) {
            // deletion SV
            // but... we don't have the reference... so we add some random sequence?

            int del_lens[] = getSVLen(INFO);
            int end_loc = getEndInfo(INFO); // TODO may need to check

            ref_deleted = REF;
            byte[] refs = new byte[0];
            pos++;

            if (del_lens.length > 0) {
                alts = new FlexSeq[del_lens.length];
                for (int i = 0; i < del_lens.length; i++) {
                    // deletion has no alt
                    alts[i] = new FlexSeq(FlexSeq.Type.DEL, 0);
                }

                return new Variant(chr, pos, Math.abs(del_lens[0]), refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else if (end_loc > 0) {
                int del_len = end_loc - pos + 1;
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.DEL, 0);
                return new Variant(chr, pos, del_len, refs, alts,
                        phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
            } else {
                log.error("No length information for DEL:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.indexOf('<') >= 0) {
            // inprecise variant
            log.warn("Imprecise line: " + line);
            return null;
        } else {

            // Splitting
            String[] alts_str = ALT.split(",");
            int n = alts_str.length; // number of alts

            alts = new FlexSeq[n];
            for (int i = 0; i < n; i++) {
                byte[] temp = new byte[alts_str[i].length()];
                for (int j = 0; j < alts_str[i].length(); j++) {
                    temp[j] = (byte) alts_str[i].charAt(j);
                }
                alts[i] = new FlexSeq(temp);
            }

            // Check
            for (int i = 0; i < n; i++) {
                if (REF.length() == 1 && alts[i].length() == 1) {
                    // SNP
                } else if (REF.length() == 0 || alts[i].length() == 0) {
                    log.warn("Skipping invalid record:");
                    log.warn(line);
                    return null;
                }
            }

            /* Adjustment of first base
             basically if first base of ref and alt match, first base of
             ref and alt will both be removed, pos will increment by 1 to
             account for the removal.
            */
             // TODO: This needs to be updated to account for multiple matching reference bases

            if (REF.length() > 0) {
                boolean same = true;
                for (int i = 0; i < n; i++) {
                    if (alts[i].length() == 0
                            || REF.charAt(0) != alts[i].charAt(0)) {
                        same = false;
                        break;
                    }
                }
                if (same) {
                    pos++;
                    ref_deleted = String.valueOf(REF.charAt(0));
                    REF = REF.substring(1);

                    //System.err.println(var_id + " before :" + ref_deleted);

                    for (int i = 0; i < n; i++) {
                        alts[i] = new FlexSeq(alts[i].substring(1));
                    }
                }
            }

            // TODO this needs to be done
            // but if we want to preserve the original VCF record, then this
            // needs modification
            if (REF.length() > 0) {
                int ref_len = REF.length();

                int min_clip_len = Integer.MAX_VALUE;
                for (int i = 0; i < n; i++) {
                    int alt_len = alts[i].length();

                    //what does clip_len represent?
                    int clip_len = 0;
                    for (int j = 0; j < alt_len; j++) {

                        // make sure there is at least something in alt
                        if (ref_len - j <= 0 || alt_len - j <= 0) {
                            clip_len = j;
                            break;
                        }
                        if (REF.charAt(ref_len - j - 1) != alts[i].charAt(alt_len - j - 1)) {
                            clip_len = j;
                            break;
                        }
                        clip_len = j + 1;
                    }

                    if (min_clip_len > clip_len) {
                        min_clip_len = clip_len;
                    }
                }

                /*
                apparently this code block is part of normalization.
                is this working properly, though? e.g. it converts
                CGTG,CG => GT,""
                 */
                if (min_clip_len > 0) {
                    REF = REF.substring(0, ref_len - min_clip_len);
                    for (int i = 0; i < n; i++) {
                        alts[i] = new FlexSeq(alts[i].substring(0,
                                alts[i].length() - min_clip_len));
                    }
                }
            }

            byte[] refs = new byte[REF.length()];

            for (int i = 0; i < REF.length(); i++) {
                refs[i] = (byte) REF.charAt(i);
            }


            return new Variant(chr, pos, refs.length, refs, alts,
                    phase_val, is_phased, var_id, FILTER, ref_deleted, _rand);
        }
    }

    public Variant parseLine() {
        String line = _line;
        readLine();

        if (line == null || line.length() == 0) {
            log.info("blank line");
            return null;
        }

        return process_line(line);
    }

}
