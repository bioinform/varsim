package com.bina.varsim.util;

//--- Java imports ---

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.VCFInfo;
import com.bina.varsim.types.variant.Variant;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.rmi.UnexpectedException;
import java.util.Random;
import java.util.StringTokenizer;

public class VCFparser extends GzFileParser<Variant> {
    private final static Logger log = Logger.getLogger(VCFparser.class.getName());

    private Random random = null;

    private int sampleIndex = -1;
    private String sampleId = null;
    private boolean isPassFilterRequired = false;
    private boolean chromLineSeen = false;

    public VCFparser() {
        sampleIndex = 10; // the first sample
    }

    /**
     * Reads a VCF file line by line
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param id       ID of individual, essentially selects a column, null to use first ID column
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, String id, boolean pass, Random rand) {
        random = rand;
        try {
            bufferedReader = new BufferedReader(new InputStreamReader(decompressStream(fileName)));
            readLine();
        } catch (Exception ex) {
            log.error("Can't open file " + fileName);
            log.error(ex.toString());
        }
        sampleId = id;
        isPassFilterRequired = pass;

        if (sampleId == null) {
            sampleIndex = 10; // the first sample
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

    //TODO: remove unused constructor
    /**
     * Reads a VCF file line by line, if there are multiple individuals, takes the first one
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, boolean pass, Random rand) {
        this(fileName, null, pass, rand);
    }

    //TODO: remove unused constructor
    /**
     * Reads a VCF file line by line, if there are multiple individuals, takes the first one
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, boolean pass) {
        this(fileName, null, pass, null);
    }

    /**
     * finds where "GT" or similar is in the VCF string so that the genotype can be read
     *
     * @param record The column in VCF that contains GT, CN, etc..
     * @param key    The key to be found "GT" or "CN" or other
     * @return the index of the key
     */
    private int getFormatKeyIndex(final String record, final String key) {
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
    boolean isPhased(String geno, byte[] vals, ChrString chr) {
        boolean isPhased = false;

        geno = geno.trim();
        boolean strangePhase = false;
        if (geno.matches("^[0-9]+$")) {
            // phase is only a single number, for haploid chromosomes
            byte val = (byte) Integer.parseInt(geno);
            if (chr.isX()) {
                vals[1] = val; // maternal
                isPhased = true;
            } else if (chr.isY()) {
                vals[0] = val; // paternal
                isPhased = true;
            } else if (chr.isMT()) {
                vals[1] = val;
                isPhased = true;
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
                    isPhased = true;
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
            isPhased = false;
        }

        return isPhased;
    }

    /**
     * takes a line from a VCF file, parse it,
     * return a Variant object
     *
     * right now the meta-info lines (beginning with ##) are
     * not tied with data line parsing. This will be corrected
     * in the future (perhaps with help of HTSJDK).
     * @param line
     * @return
     */
    public Variant processLine(String line) throws UnexpectedException {

        // try to determine the column we should read for the genotype
        StringTokenizer toks = new StringTokenizer(line);
        if (line.startsWith("#")) {
            if (sampleId != null && line.startsWith("#CHROM")) {
                chromLineSeen = true;
                int index = 0;
                while (toks.hasMoreTokens()) {
                    index++;
                    String tok = toks.nextToken();
                    if (tok.equals(sampleId))
                        sampleIndex = index;
                }
            } else if (sampleId == null) {
                sampleIndex = 10; // the first sample
            }
            return null;
        }

        // If we cannot determine, then use the first one
        if (sampleIndex < 0 && !chromLineSeen) {
            sampleIndex = 10;
        } else if (sampleIndex < 0) {
            sampleIndex = 10;
            log.warn("Warning!!! ID (" + sampleId + ") does not exist... ");
        }


        int index = 0, genotypeIndex = -1, copyNumberIndex = -1;
        int pos = -1;
        ChrString chr = null;
        String REF = "", FILTER = "", ALT = "", variantId = "";
        String phase = ".", copyNumber = "0/0", infoString = "", FORMAT;
        String[] sampleInfo;
        while (toks.hasMoreTokens()) {
            index++;
            if (index == 1) { // Parsing chromosome
                chr = new ChrString(toks.nextToken());
            } else if (index == 2) // Parsing position
                pos = Integer.parseInt(toks.nextToken());
            else if (index == 3) // Parsing position
                variantId = toks.nextToken();
            else if (index == 4) // Parsing reference allele
                REF = toks.nextToken();
            else if (index == 5) // Parsing alternative allele
                ALT = toks.nextToken();
            else if (index == 7) // FILTER field
                FILTER = toks.nextToken();
            else if (index == 8) // INFO field
                infoString = toks.nextToken();
            else if (index == 9) { // Output format
                FORMAT = toks.nextToken();
                genotypeIndex = getFormatKeyIndex(FORMAT, "GT");
                copyNumberIndex = getFormatKeyIndex(FORMAT, "CN");
            } else if (index == sampleIndex) { // phased or unphased genotype
                sampleInfo = (toks.nextToken()).split(":");
                if (genotypeIndex >= 0) {
                    phase = sampleInfo[genotypeIndex];
                }
                if (copyNumberIndex >= 0) {
                    copyNumber = sampleInfo[copyNumberIndex];
                }

                break;
            } else {
                toks.nextToken();
            }
        }

        // unknown chromosome
        // TODO: throw an exception for unknown chromosome name
        if (chr == null) {
            log.warn("Bad chromosome name: " + line);
            return null;
        }

        if (isPassFilterRequired && !(FILTER.contains("PASS") || FILTER.equals("."))) {
            //log.warn("not pass line" + line);
            return null; // Filtered out
        }
        // parse the phased or unphased genotype
        byte[] genotypeArray = new byte[2]; // paternal-maternal
        boolean isGenotypePhased = isPhased(phase, genotypeArray, chr);


        if (genotypeIndex >= 0 && genotypeArray[0] == 0 && genotypeArray[1] == 0) {
            //return null; // reference alleles... ignore them for now....
        }

        // determine copy-number
        // TODO need to be able to deal with unphased copy-numbers?
        byte[] copyNumberArray = new byte[2]; // paternal-maternal
        boolean isCopyNumberPhased;

        if (copyNumberIndex >= 0) {
            isCopyNumberPhased = isPhased(copyNumber, copyNumberArray, chr);
            if (isCopyNumberPhased != isGenotypePhased) {
                // TODO maybe don't throw error, this is not standard format
                // anyways
                log.error("Inconsistent copy number:");
                log.error(line);
                return null;
            }
            if ((genotypeArray[0] == genotypeArray[1]) != (copyNumberArray[0] == copyNumberArray[1])) {
                throw new IllegalArgumentException("ERROR: genotype does not agree with copy number.\n" + line);
            }
        }

        // Upper casing
        REF = REF.toUpperCase();
        ALT = ALT.toUpperCase();

        String deletedReference = "";
        FlexSeq alts[];
        VCFInfo info = new VCFInfo(infoString);

        /*if symbolic alleles are present,
        make sure # of alleles equal # of
        SV lengths. For non-symbolic alleles, SV lengths
        are not really used or checked.
         */
        /*!!!!!!!!!!
        CAUTION: we assume symbolic alleles are not mixed
        with non-symbolic alleles.
         */
        if (ALT.indexOf('<') != -1) {
            String[] alternativeAlleles = ALT.split(",");
            int[] svlen = (int[]) info.getValue("SVLEN");
            if (alternativeAlleles.length != svlen.length) {
                throw new IllegalArgumentException("ERROR: number of symbolic alleles is unequal to number of SV lengths.\n" + line);
            }
            for (int i = 0; i < alternativeAlleles.length; i++) {
                if (!alternativeAlleles[i].startsWith("<")) {
                    throw new IllegalArgumentException("ERROR: symbolic alleles are mixed with non-symbolic alleles.\n" + line);
                }
            }
        }
        if (ALT.startsWith("<INV>")) {
            // inversion SV

            int[] end = (int[]) info.getValue("END");
            int[] svlen = (int[]) info.getValue("SVLEN");

            deletedReference = REF;
            byte[] refs = new byte[0];
            pos++;

            if (svlen.length > 0) {
                alts = new FlexSeq[svlen.length];
                for (int i = 0; i < svlen.length; i++) {
                    int alternativeAlleleLength = Math.max(Math.abs(svlen[i]), 1);
                    alts[i] = new FlexSeq(FlexSeq.Type.INV, alternativeAlleleLength);
                }
                // TODO this assumes only one alt
                /*return new Variant(chr, pos, Math.abs(svlen[0]), refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(Math.abs(svlen[0])).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
                //TODO: this assumes only one alt, might not be true
            } else if (end[0] > 0) {
                int alternativeAlleleLength = Math.max(Math.abs(end[0] - pos + 1), 1);
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.INV, alternativeAlleleLength);
                /*return new Variant(chr, pos, alternativeAlleleLength, refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(alternativeAlleleLength).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
            } else {
                log.error("No length information for INV:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.startsWith("<DUP>") || ALT.startsWith("<DUP:TANDEM>")) {
            // duplication or tandem duplication SV
            int[] svlens = (int[]) info.getValue("SVLEN");
            // may need to check inconsistency
            int[] ends = (int[]) info.getValue("END");
            // btw length and end location

            deletedReference = REF;
            byte[] refs = new byte[0];
            pos++;

            if (svlens.length > 0) {
                alts = new FlexSeq[svlens.length];
                for (int i = 0; i < svlens.length; i++) {
                    // TODO this is temporary, how to encode copy number?
                    int currentCopyNumber = 1;
                    for (int j = 0; j < 2; j++) {
                        if ((i + 1) == genotypeArray[j]) {
                            /*
                            if i = 0, genotype[0] = 1, genotype[1] = 1
                            copyNumberArray[0] = 3, copyNumberArray[1] = 2
                            then currentCopyNumber = 2.
                            what does currentCopyNumber mean in real world?
                             */
                            if (copyNumberArray[j] > 0) {
                                currentCopyNumber = copyNumberArray[j];
                            }
                        }
                    }

                    int alternativeAlleleLength = Math.max(Math.abs(svlens[i]), 1);

                    alts[i] = new FlexSeq(FlexSeq.Type.DUP, alternativeAlleleLength, currentCopyNumber);
            }

                /*return new Variant(chr, pos, Math.abs(svlens[0]), refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(Math.abs(svlens[0])).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
                //TODO: this assumes only one alt, which might not be true
            } else if (ends[0] > 0) {
                int alternativeAlleleLength = Math.max(Math.abs(ends[0] - pos + 1), 1);
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.DUP, alternativeAlleleLength, Math.max(
                        copyNumberArray[0], copyNumberArray[1]));

                /*return new Variant(chr, pos, alternativeAlleleLength, refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(alternativeAlleleLength).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
            } else {
                log.error("No length information for DUP:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.startsWith("<INS>")) {
            // insertion SV

            int[] ends = (int[]) info.getValue("END");
            // TODO may need to check
            int[] svlens = (int[]) info.getValue("SVLEN");

            deletedReference = REF;
            byte[] refs = new byte[0];
            pos++;

            if (svlens.length > 0) {
                alts = new FlexSeq[svlens.length];
                for (int i = 0; i < svlens.length; i++) {
                    int alternativeAlleleLength;
                    if (svlens[i] == 0) {
                        alternativeAlleleLength = Integer.MAX_VALUE;
                    } else {
                        alternativeAlleleLength = Math.max(Math.abs(svlens[i]), 1);
                    }
                    alts[i] = new FlexSeq(FlexSeq.Type.INS, alternativeAlleleLength);
                }
                /*return new Variant(chr, pos, 0, refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(0).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
            } else if (ends[0] > 0) {
                int alternativeAlleleLength = Math.max(Math.abs(ends[0] - pos), 1);
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.INS, alternativeAlleleLength);
                /*return new Variant(chr, pos, 0, refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(0).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
            } else {
                log.error("No length information for INS:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.startsWith("<DEL>")) {
            // deletion SV
            // but... we don't have the reference... so we add some random sequence?

            int[] ends = (int[]) info.getValue("END");
            // TODO may need to check
            int[] svlens = (int[]) info.getValue("SVLEN");

            deletedReference = REF;
            byte[] refs = new byte[0];
            pos++;

            if (svlens.length > 0) {
                alts = new FlexSeq[svlens.length];
                for (int i = 0; i < svlens.length; i++) {
                    // deletion has no alt
                    alts[i] = new FlexSeq(FlexSeq.Type.DEL, 0);
                }

                /*return new Variant(chr, pos, Math.abs(svlens[0]), refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(Math.abs(svlens[0])).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
            } else if (ends[0] > 0) {
                int alternativeAlleleLength = ends[0] - pos + 1;
                alts = new FlexSeq[1];
                alts[0] = new FlexSeq(FlexSeq.Type.DEL, 0);
                /*return new Variant(chr, pos, alternativeAlleleLength, refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(alternativeAlleleLength).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).build();
            } else {
                log.error("No length information for DEL:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.startsWith("<TRA>")) {
            //translocation SV

            //1-based end for reference allele
            int[] end = (int[]) info.getValue("END");
            int[] end2 = (int[]) info.getValue("END2");
            int[] pos2 = (int[]) info.getValue("POS2");
            //SVLEN for alternative allele length
            int[] svlen = (int[]) info.getValue("SVLEN");
            String[] chr2 = (String[]) info.getValue("CHR2");
            String[] subtype = (String[]) info.getValue("TRASUBTYPE");

            deletedReference = REF;
            byte[] refs = new byte[0];
            pos++; //1-based start

            if (svlen.length > 0) {
                //0 is for reference allele
                alts = new FlexSeq[svlen.length];
                //alternative allele is numbered 1,2,... per VCFSpec
                for (int altAlleleIndex = 1; altAlleleIndex <= svlen.length; altAlleleIndex++) {
                    int currentCopyNumber = 1;
                    /*
                    implicit assumption here: genotype[0] == genotype[1] => copyNumberArray[0] == copyNumberArray[1]
                     */
                    //check paternal
                    if (altAlleleIndex == genotypeArray[0]) {
                        currentCopyNumber = copyNumberArray[0];
                    }
                    //check maternal
                    if (altAlleleIndex == genotypeArray[1]) {
                        currentCopyNumber = copyNumberArray[1];
                    }
                    currentCopyNumber = Math.max(1, currentCopyNumber);
                    int altAllelelength = Math.max(Math.abs(svlen[altAlleleIndex - 1]), 1);
                    /*Not sure if INS is the most appropriate FlexSeq.Type
                    use it unless it breaks things.
                     */
                    if (subtype[altAlleleIndex - 1].equals("ACCEPT")) {
                        alts[altAlleleIndex - 1] = new FlexSeq(FlexSeq.Type.TRA, altAllelelength, currentCopyNumber);
                    } else if (subtype[altAlleleIndex - 1].equals("REJECT")){
                        alts[altAlleleIndex - 1] = new FlexSeq(FlexSeq.Type.DEL, 0);
                    } else {
                        throw new IllegalArgumentException("ERROR: only ACCEPT and REJECT allowed.\n" + line);
                    }
                }

                //TODO: there could be multiple SVLEN values, but for now we only use one
                //make sure all SVLEN values are equal if we only use the first one
                //per VCFv4.1 spec, SVLEN should be length of alternative allele rather than
                //reference allele
                for (int i = 1; i < svlen.length; i++) {
                    if (svlen[i] != svlen[0]) {
                        throw new IllegalArgumentException("ERROR: SVLEN values not equal.\n" + line);
                    }
                }
                //pos is incremented by 1, so it becomes 1-based start
                /*return new Variant(chr, pos, Math.abs(end[0] - pos + 1), refs, alts,
                        genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random, ChrString.string2ChrString(chr2), pos2, end2, end[0], subtype);
                        */
                return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(Math.abs(end[0] - pos + 1)).
                        ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                        varId(variantId).filter(FILTER).refDeleted(deletedReference).
                        randomNumberGenerator(random).chr2(ChrString.string2ChrString(chr2)).
                        pos2(pos2).end2(end2).end(end[0]).translocationSubtype(subtype).build();
                //TODO: this assumes only one alt, which might not be true
            } else {
                log.error("No length information for TRA:");
                log.error(line);
                log.error("skipping...");
                return null;
            }
        } else if (ALT.indexOf('<') >= 0) {
            // imprecise variant
            log.warn("Imprecise line: " + line);
            return null;
        } else {

            // Splitting
            String[] alternativeAlleles = ALT.split(",");
            int n = alternativeAlleles.length; // number of alts

            alts = new FlexSeq[n];
            for (int i = 0; i < n; i++) {
                byte[] temp = new byte[alternativeAlleles[i].length()];
                for (int j = 0; j < alternativeAlleles[i].length(); j++) {
                    temp[j] = (byte) alternativeAlleles[i].charAt(j);
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
            //e.g. ref=AT alt=ATC (VCFv4.1 spec requires only 1-bp before event, but it might
            //not be the case all the time.

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
                    deletedReference = String.valueOf(REF.charAt(0));
                    REF = REF.substring(1);

                    //System.err.println(varId + " before :" + deletedReference);

                    for (int i = 0; i < n; i++) {
                        alts[i] = new FlexSeq(alts[i].substring(1));
                    }
                }
            }

            // TODO this needs to be done
            // but if we want to preserve the original VCF record, then this
            // needs modification
            if (REF.length() > 0) {
                int referenceAlleleLength = REF.length();

                int minClipLength = Integer.MAX_VALUE;
                for (int i = 0; i < n; i++) {
                    int alternativeAlleleLength = alts[i].length();

                    //what does clipLength represent?
                    int clipLength = 0;
                    for (int j = 0; j < alternativeAlleleLength; j++) {

                        // make sure there is at least something in alt
                        if (referenceAlleleLength - j <= 0 || alternativeAlleleLength - j <= 0) {
                            clipLength = j;
                            break;
                        }
                        if (REF.charAt(referenceAlleleLength - j - 1) != alts[i].charAt(alternativeAlleleLength - j - 1)) {
                            clipLength = j;
                            break;
                        }
                        clipLength = j + 1;
                    }

                    if (minClipLength > clipLength) {
                        minClipLength = clipLength;
                    }
                }

                /*
                apparently this code block is part of normalization.
                is this working properly, though? e.g. it converts
                CGTG,CG => GT,""
                 */
                if (minClipLength > 0) {
                    REF = REF.substring(0, referenceAlleleLength - minClipLength);
                    for (int i = 0; i < n; i++) {
                        alts[i] = new FlexSeq(alts[i].substring(0,
                                alts[i].length() - minClipLength));
                    }
                }
            }

            byte[] refs = new byte[REF.length()];

            for (int i = 0; i < REF.length(); i++) {
                refs[i] = (byte) REF.charAt(i);
            }


            /*return new Variant(chr, pos, refs.length, refs, alts,
                    genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                    */
            return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(refs.length).
                    ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                    varId(variantId).filter(FILTER).refDeleted(deletedReference).
                    randomNumberGenerator(random).build();
        }
    }

    public Variant parseLine() {
        /*
        TODO: line reading should be handled in a loop calling
         */
        String line = this.line;
        readLine();

        if (line == null || line.length() == 0) {
            log.info("blank line");
            return null;
        }

        try {
            return processLine(line);
        } catch (Exception e) {
            //TODO: right now just be lazy, die on any error
            log.error(e.getMessage());
            System.exit(255);
        }
        return null;
    }

}
