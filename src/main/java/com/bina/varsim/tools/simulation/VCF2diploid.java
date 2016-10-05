package com.bina.varsim.tools.simulation;

//--- Java imports ---

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.GenderType;
import com.bina.varsim.types.Sequence;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.types.variant.VariantOverallType;
import com.bina.varsim.types.variant.VariantType;
import com.bina.varsim.util.SimpleReference;
import com.bina.varsim.util.VCFparser;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.*;

/**
 * Class to construct diploid genome from genome reference and genome variants
 * calls in VCF format.
 *
 * @author Alexej Abyzov, John C. Mu
 */

public class VCF2diploid {

    static final long SEED_ARG = 3333;
    private final static Logger log = Logger.getLogger(VCF2diploid.class.getName());
    private final static char DELETED_BASE = '~';
    private final static String[] DIPLOID_CHRS = {"maternal", "paternal"};

    private File outputMap;

    // arguments
    @Option(name = "-t", usage = "Gender of individual [MALE]")
    GenderType gender = GenderType.MALE;
    @Option(name = "-chr", usage = "Comma separated list of reference genome files [Required]", metaVar = "FASTA_file", required = true)
    List<String> chrfiles = null;
    @Option(name = "-vcf", usage = "Comma separated list of VCF files [Required]", metaVar = "VCF_file", required = true)
    List<String> vcfFiles = null;
    @Option(name = "-seed", usage = "Seed for random sampling [" + SEED_ARG + "]")
    long seed = 3333;
    private Random rand = null;
    @Option(name = "-id", usage = "ID of individual in VCF file [Optional]")
    private String id = "varsim";
    @Option(name = "-pass", usage = "Only accept the PASS variants")
    private boolean pass = false;
    @Option(name = "-outdir", usage = "Directory to output results in [current directory]")
    File outDir = new File("").getAbsoluteFile();
    private Map<ChrString, List<Variant>> variants = new HashMap<>();

    /**
     * set seed for random number generation
     * seed comes from the command line
     */
    public VCF2diploid() {
        rand = new Random(seed);
    }

    /**
     * main method
     *
     * @param args
     */
    public static void main(String[] args) {
        VCF2diploid runner = new VCF2diploid();
        runner.run(args);
    }

    public File getOutputMap() {
        return outputMap;
    }

    /**
     * parse arguments, aggregate all variants from all VCFs
     * then create perturbed diploid genomes (paternal + maternal)
     *
     * @param args
     */
    public void run(String[] args) {
        String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();

        String usage = "Create a diploid genome as associated files from a reference genome\n"
                + "and some VCF files. \n";

        CmdLineParser cmd_parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        cmd_parser.setUsageWidth(80);

        try {
            cmd_parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java -jar vcf2diploid.jar [options...]");
            // print the list of available options
            cmd_parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }

        if (chrfiles.size() == 0) {
            log.error("No chromosome file(s) is given!\n" + usage);
            return;
        }

        if (vcfFiles.size() == 0) {
            log.error("No VCF file(s) is given!");
        }

        if (!outDir.isDirectory()) {
            log.info("Creating output directory " + outDir);
            outDir.mkdirs();
        }

        parseVCFs(vcfFiles, variants, id, pass);
        makeDiploid();
    }

    /**
     * parse all VCFs for a specific sample with desired filter
     * store them in variants by chromosome
     * @param vcfs
     * @param id sample ID
     * @param pass true only retain PASS variants, false otherwise
     */
    public void parseVCFs(List<String> vcfs, Map<ChrString, List<Variant>> variants, String id, boolean pass) {
        for (String _vcfFile : vcfs) {
            VCFparser parser = new VCFparser(_vcfFile, id, pass, rand);
            /*
            apparently, nVariant is for # of variants
            nVariantBase is for # of nucleotides in variants
            nVariantBase is an approximate indicator for length of variants
            because it involves both reference allele and alternative
            alleles.
             */
            int nVariant = 0, nVariantBase = 0;
            while (parser.hasMoreInput()) {
                Variant var = parser.parseLine();
                if (var == null) {
                    // System.err.println("Bad variant or not a variant line");
                    continue;
                }
                if (var.maternal() == 0 && var.paternal() == 0) {
                    log.warn("Not maternal nor paternal");
                    continue;
                }

                if (var.maternal() < 0 || var.paternal() < 0) {
                    // randomize the genotype
                    var.randomizeGenotype(gender);
                }

                ChrString chr = var.getChr();
                if (!variants.containsKey(chr)) {
                    variants.put(chr, new ArrayList<Variant>());
                }
                variants.get(chr).add(var);
                nVariant++;
                nVariantBase += var.variantBases();
            }
            System.out.println(_vcfFile + ": " + nVariant + " variants, "
                    + nVariantBase + " variant bases");
        }
    }

    /**
     * This is the main function that makes the diploid genome
     * for each chromosome/sequence, all variants are traversed.
     * Each variant is decomposed into a combination of deletion
     * and insertion. Original sequence is flagged at variant loci.
     * Perturbed sequences will be generated based on these flags.
     * Map file format (MFF) records will generated by iterating
     * over flag array of each sequence. VCFs (one for each sequence)
     * containing all sampled variants will be generated as well
     * with optional phasing information (if both maternal and
     * paternal genomes are chosen for output).
     *
     */
    public void makeDiploid() {
        StringBuilder map_string = new StringBuilder();

        // This is the loop if chromosomes exist in separate files
        SimpleReference allSequences = new SimpleReference(chrfiles);

        // This is the loop through each chromosome
        for (ChrString chr : allSequences.keySet()) {
            Sequence referenceSequence = allSequences.getSequence(chr);

            System.out.println("Working on " + referenceSequence.getName() + "...");

            boolean output_paternal = true;
            boolean output_maternal = true;

            if (gender == GenderType.FEMALE) {
                if (chr.isY()) {
                    output_paternal = false;
                    output_maternal = false;
                }
            } else if (gender == GenderType.MALE) {
                // only male and female
                if (chr.isX()) {
                    output_paternal = false;
                }
                if (chr.isY()) {
                    output_maternal = false;
                }
            }
            if (chr.isMT()) {
                output_paternal = false;
            }

            if (!output_paternal && !output_maternal) {
                // skip chromosome
                continue;
            }

            // this is the list of variants for the chromosome of question
            final List<Variant> varList = variants.containsKey(chr) ? variants.get(chr) : Collections.EMPTY_LIST;

            final List<Boolean> maternalIsVariantAdded = new ArrayList<>();
            final List<Boolean> paternalIsVariantAdded = new ArrayList<>();

            int len = referenceSequence.length();
            byte[] maternalMaskedSequence = new byte[len];
            byte[] paternalMaskedSequence = new byte[len];
            // byte[] ins_flag = new byte[len];
            // Flag specification:
            // b -- insertion in both haplotypes
            // p -- insertion in paternal haplotype
            // m -- insertion in maternal haplotype
            // atcgATCG -- original nucleotides from reference sequence

            // fill both maternal and paternal with the original reference
            // sequence
            for (int c = 1; c <= len; c++) {
                maternalMaskedSequence[c - 1] = paternalMaskedSequence[c - 1] = referenceSequence.byteAt(c);
            }

            Hashtable<Integer, FlexSeq> paternalInsertionSeq = new Hashtable<>(150);
            Hashtable<Integer, FlexSeq> maternalInsertionSeq = new Hashtable<>(150);

            int nPaternalVariant = 0, nMaternalVariant = 0;
            int nPaternalVariantBase = 0, nMaternalVariantBase = 0;
            for (Variant var : varList) {
                // iterate over the variants in the chromosome
                if (!var.isPhased()) {
                    var.randomizeHaplotype();
                }

                if (var.paternal() > 0) {
                    if (addVariant(paternalMaskedSequence, referenceSequence, var, var.paternal(),
                            paternalInsertionSeq)) {
                        nPaternalVariant++;
                        nPaternalVariantBase += var.variantBases();
                        paternalIsVariantAdded.add(true);
                    } else {
                        paternalIsVariantAdded.add(false);
                    }
                } else {
                    paternalIsVariantAdded.add(false);
                }

                if (var.maternal() > 0) {
                    if (addVariant(maternalMaskedSequence, referenceSequence, var, var.maternal(),
                            maternalInsertionSeq)) {
                        nMaternalVariant++;
                        nMaternalVariantBase += var.variantBases();
                        maternalIsVariantAdded.add(true);
                    } else {
                        maternalIsVariantAdded.add(false);
                    }
                } else {
                    maternalIsVariantAdded.add(false);
                }

            }

            log.info("number of variants: " + varList.size());

            //VCF is written on a per-chromosome basis
            writeVCF(referenceSequence, varList, paternalIsVariantAdded,
                    maternalIsVariantAdded, output_paternal, output_maternal);

            if (output_paternal) {
                String paternalSequenceName = referenceSequence.getName() + "_" + DIPLOID_CHRS[1];
                String paternalSequenceFileName = referenceSequence.getName() + "_" + id + "_" + DIPLOID_CHRS[1] + ".fa";
                makePosMap(map_string, paternalSequenceName, referenceSequence, paternalMaskedSequence, paternalInsertionSeq);
                writeHaploid(paternalMaskedSequence, paternalInsertionSeq, paternalSequenceName, paternalSequenceFileName);
                System.out.println("Applied " + nPaternalVariant + " variants "
                        + nPaternalVariantBase + " bases to " + "paternal genome.");
            }

            if (output_maternal) {
                String maternalSequenceName = referenceSequence.getName() + "_" + DIPLOID_CHRS[0];
                String maternalSequenceFileName = referenceSequence.getName() + "_" + id + "_" + DIPLOID_CHRS[0] + ".fa";
                makePosMap(map_string, maternalSequenceName, referenceSequence, maternalMaskedSequence, maternalInsertionSeq);
                writeHaploid(maternalMaskedSequence, maternalInsertionSeq, maternalSequenceName, maternalSequenceFileName);
                System.out.println("Applied " + nMaternalVariant + " variants "
                        + nMaternalVariantBase + " bases to " + "maternal genome.");
            }
        }

            try {
                FileWriter fw = new FileWriter(new File(outDir, id + ".map"));
                BufferedWriter bw = new BufferedWriter(fw);
                bw.write(map_string.toString());
                bw.newLine();
                bw.close();
                fw.close();
                outputMap = new File(outDir, id + ".map");
            } catch (IOException ex) {
                log.error(ex.toString());
            }
    }

    /**
     * add variant onto the perturbed genome
     * maskedSequence consists of base pairs (if no deletion occurs)
     * and flags for deletions. for each insertion,
     * position and inserted sequence will be put into InsertPosition2Sequence.
     * Note here, a SNP is modeled as 1-bp del+1-bp insertion,
     * MVN is modeled as multiple 1-bp del + multiple 1-bp
     * insertion. Inversion modeled as deletion + insertion of
     * reverse complemented seq. Tandem duplication as deletion
     * and insertion of multiple copies of sequences.
     *
     * variants out of bound or overlap with previous variants
     * will be discarded.
     *
     * @param maskedSequence -- sequence to be modified, it's a masked version of original sequence in array data structure
     * @param referenceSequence -- reference sequence
     * @param InsertPosition2Sequence --  records insertion locations and sequences
     * @return true if the variant is incorporated, false otherwise
     */
    private boolean addVariant(byte[] maskedSequence, Sequence referenceSequence,
                               Variant variant, int allele, Hashtable<Integer, FlexSeq> InsertPosition2Sequence) {

        //position     -- position of the variant
        //start is 1-based
        int position = variant.getPos();
        //referenceAlleleLength     -- length of reference sequence of the getReferenceAlleleLength
        //for a SNP, referenceAlleleLength == 1, so instead of substitution, a SNP
        //is modeled as one-bp referenceAlleleLength + one-bp insertion
        int referenceAlleleLength = variant.getReferenceAlleleLength();
        //insertions -- the insertion string
        byte[] insertions = variant.insertion(allele);
        VariantType variantType = variant.getType(allele);

        boolean overlap = false;

        /*
        assuming 1-based index, (position + referenceAlleleLength - 1 ) - position + 1 = referenceAlleleLength
        so position + referenceAlleleLength - 1 (instead of position + referenceAlleleLength) is the correct
        index for end, and it should not exceed length of original
        sequence.
         */
        if (position > maskedSequence.length || position + referenceAlleleLength - 1> maskedSequence.length) {
            log.warn("Variant out of chromosome bounds at "
                    + referenceSequence.getName() + ":" + position + ", (referenceAlleleLength,insertions) of (" + referenceAlleleLength
                    + "," + Arrays.toString(insertions) + "). Skipping.");
            return false;
        }

        /*
        variants in the same haplotype should never overlap.
        Otherwise, they should be merged at first place.
         */
        for (int p = position; p < position + referenceAlleleLength; p++) {
            // if any location of this variant overlap a deleted base or a SNP, we skip it
            // DELETED_BASE is used as a flag to mark that this position has been processed/modified
            if (maskedSequence[p - 1] == DELETED_BASE
                    // insertions may not be next to SNPs or deletions, this is to avoid problem of DUP/INV
                    //why p - 1 for maskedSequence? is it 1-based or 0-based?
                    || maskedSequence[p - 1] != referenceSequence.byteAt(p)) {
                overlap = true;
            }
        }

        if (overlap) {
            try {
                if (insertions != null) {
                    log.warn("Variant overlap at " + referenceSequence.getName() + ":"
                            + position + ", (referenceAlleleLength,insertions) of (" + referenceAlleleLength + "," + new String(insertions, "US-ASCII") + "). Skipping.");
                } else {
                    log.warn("Variant overlap at " + referenceSequence.getName() + ":"
                            + position + ", (referenceAlleleLength,insertions) of (" + referenceAlleleLength + ",<imprecise>). Skipping.");
                }
            } catch (UnsupportedEncodingException e) {
                e.printStackTrace();
            }
            return false;
        }

        if (referenceAlleleLength == 1 && (insertions != null && insertions.length == 1)) { // SNP
            // this is to maintain the reference repeat annotations
            if (Character.isLowerCase((char) referenceSequence.byteAt(position))) {
                maskedSequence[position - 1] = (byte) Character.toLowerCase((char) insertions[0]);
            } else {
                maskedSequence[position - 1] = (byte) Character.toUpperCase((char) insertions[0]);
            }
        } else if (variantType == VariantType.MNP) {
            // add this as a bunch of SNPs
            // TODO this may result in other variants getting added in between
            for (int i = position; i < position + referenceAlleleLength; i++) {
                // add each SNP
                assert insertions != null;
                if (Character.isLowerCase((char) referenceSequence.byteAt(i))) {
                    maskedSequence[i - 1] = (byte) Character.toLowerCase((char) insertions[i - position]);
                } else {
                    maskedSequence[i - 1] = (byte) Character.toUpperCase((char) insertions[i - position]);
                }
            }

        } else { // Indel, SV
            // check whether the insertion location has been inserted before
            if (InsertPosition2Sequence.get(position) != null) {
                log.warn("Multiple insertions at "
                        + referenceSequence.getName() + ":" + position);
                try {
                    if (insertions != null) {
                        log.warn("Skipping variant with (referenceAlleleLength,insertions) of (" + referenceAlleleLength
                                + "," + new String(insertions, "US-ASCII") + ").");
                    } else {
                        log.warn("Skipping variant with (referenceAlleleLength,insertions) of (" + referenceAlleleLength
                                + ", <imprecise> ).");
                    }
                } catch (UnsupportedEncodingException e) {
                    e.printStackTrace();
                }
                return false;
            }

            if (variantType == VariantType.Inversion) {
                // Treat this as getReferenceAlleleLength then insertion of reverse complement
                //System.err.println("Insert INV");
                referenceAlleleLength = variant.insertion_len(allele);
                insertions = referenceSequence.revComp(position, position + referenceAlleleLength);
            }

            if (variantType == VariantType.Tandem_Duplication) {
                // treat this as getReferenceAlleleLength then insertion of several copies
                // this prevents the original sequence to be altered and
                // become less like a duplication
                // TODO make length correct
                // System.err.println("Insert DUP");

                referenceAlleleLength = variant.insertion_len(allele);
                int single_ins_len = variant.insertion_len(allele);
                byte[] orig_seq = referenceSequence.subSeq(position, position + single_ins_len);
                insertions = new byte[single_ins_len * variant.getCN(allele)];
                System.arraycopy(orig_seq, 0, insertions, 0, orig_seq.length);
                for (int i = 0; i < insertions.length; i++) {
                    insertions[i] = orig_seq[i % (orig_seq.length)];
                }
            }

            // set to deleted base, so that we don't change those in future
            for (int p = position; p < position + referenceAlleleLength; p++) {
                maskedSequence[p - 1] = DELETED_BASE;
            }

            // matter
            // TODO if insertions is null we still need to add??
            if (insertions != null && insertions.length > 0) {
                // convert to flexseq
                FlexSeq s = new FlexSeq(insertions);
                s.setType(variant.getAlt(allele).getType());
                s.setCopy_num(variant.getAlt(allele).getCopy_num());
                s.setLength(variant.getAlt(allele).length());
                s.setVar_id(variant.getVar_id());

                /*
                so although we model SNP as deletion+insertions, we do not record
                their insert locus. We only record insert loci of indel
                and SVs.
                 */
                InsertPosition2Sequence.put(new Integer(position), s);
            }
        }

        return true;
    }

    /**
     * adjust indeces(positions) of host genome and reference genome
     * based on variant type
     * @param curr_rec current record in map file
     * @param hf_idx host to reference genome index
     */
    private void adjust_idx(MapRecord curr_rec, HostRefIdx hf_idx) {
        switch (curr_rec.feature) {
            case "SEQ":
                hf_idx.host_idx += curr_rec.len;
                hf_idx.ref_idx += curr_rec.len;
                break;
            case "DEL":
                hf_idx.ref_idx += curr_rec.len;
                break;
            case "INS":
                hf_idx.host_idx += curr_rec.len;
                break;
            case "DUP_TANDEM":
                hf_idx.host_idx += curr_rec.len;
                break;
            case "INV":
                hf_idx.host_idx += curr_rec.len;
                break;
        }
    }

    /**
     * create a new map file record based on the variant or sequence object stored at ins_seq[idx]
     * from the perspective of good programming practice, it's better to iterate over all variants
     * in ins_seq rather than using idx.
     * after the record is created, advance host genome and reference indices so we know where next
     * event occurs. For some variants, e.g. TANDEM_DUP, several records will be created on the fly
     * and got appended to output string sb.
     * TODO: separate creation of new map file records with appending to output string (for simplification)
     *
     * @param sb output string
     * @param idx position of examination
     * @param chr_name
     * @param ref_chr_name
     * @param hf_idx
     * @param genome
     * @param ins_seq
     * @return
     */
    private MapRecord new_curr_rec(StringBuilder sb, int idx, String chr_name, String ref_chr_name, HostRefIdx hf_idx,
                                   byte[] genome, Hashtable<Integer, FlexSeq> ins_seq) {
        MapRecord curr_rec = new MapRecord();
        curr_rec.host_chr = chr_name;
        curr_rec.ref_chr = ref_chr_name;

        // compute what the first line is and store as object

        // if it is inserted, we copy the var_id to the getReferenceAlleleLength
        boolean inserted = false;
        String var_id = ".";
        if (ins_seq.containsKey(idx + 1)) {
            inserted = true;
            // insertion at this location
            FlexSeq ins = ins_seq.get(idx + 1);
            var_id = ins.getVar_id();

            //System.err.println("Check type: " + ins.getType());

            // need to check what kind of insertion it is
            switch (ins.getType()) {
                case SEQ:
                case INS:
                    curr_rec.host_pos = hf_idx.host_idx;
                    curr_rec.ref_pos = hf_idx.ref_idx - 1;
                    curr_rec.feature = "INS";
                    curr_rec.dir = true;
                    curr_rec.len = ins.var_length();
                    curr_rec.var_id = var_id;

                    break;
                case INV:
                    // TODO: treat inversion like MNP
                    curr_rec.host_pos = hf_idx.host_idx;
                    /*position of inversion on reference genome is questionable
                    inversion does not change length of reference, so is it
                    still necessary to set its position on reference as 1bp
                    before the event?
                    */
                    curr_rec.ref_pos = hf_idx.ref_idx - 1;
                    curr_rec.feature = "INV";
                    //why direction is false (negative strand)?
                    curr_rec.dir = false;
                    curr_rec.len = ins.var_length();
                    curr_rec.var_id = var_id;

                    break;
                case DUP:
                    // need to replicate several blocks
                    int cn = ins.getCopy_num();

                    // first build one
                    curr_rec.host_pos = hf_idx.host_idx;
                    curr_rec.ref_pos = hf_idx.ref_idx - 1;
                    curr_rec.feature = "DUP_TANDEM";
                    curr_rec.dir = true;
                    curr_rec.len = ins.length();
                    curr_rec.var_id = var_id;

                    // iterate
                    for (int i = 1; i < cn; i++) {
                        adjust_idx(curr_rec, hf_idx);
                        sb.append(curr_rec);
                        sb.append('\n');
                        curr_rec = new MapRecord();
                        curr_rec.host_chr = chr_name;
                        curr_rec.ref_chr = ref_chr_name;
                        curr_rec.host_pos = hf_idx.host_idx;
                        curr_rec.ref_pos = hf_idx.ref_idx - 1;
                        curr_rec.feature = "DUP_TANDEM";
                        curr_rec.dir = true;
                        curr_rec.len = ins.length();
                        curr_rec.var_id = var_id;
                    }

                    break;
            }

            // output it
            adjust_idx(curr_rec, hf_idx);
            sb.append(curr_rec);
            sb.append('\n');

            curr_rec = new MapRecord();
            curr_rec.host_chr = chr_name;
            curr_rec.ref_chr = ref_chr_name;

        }

        if (genome[idx] == DELETED_BASE) {
            // deleted base
            curr_rec.host_pos = hf_idx.host_idx - 1;
            curr_rec.ref_pos = hf_idx.ref_idx;
            curr_rec.feature = "DEL";
            curr_rec.dir = true;
            /*
            although length is 1 here, this is just the beginning
            of a block, length may increase later (after current
            method returns).
             */
            curr_rec.len = 1;
            if (inserted) {
                curr_rec.var_id = var_id;
            }
        } else {
            // regular sequence
            curr_rec.host_pos = hf_idx.host_idx;
            curr_rec.ref_pos = hf_idx.ref_idx;
            curr_rec.feature = "SEQ";
            curr_rec.dir = true;
            /*
            although length is 1 here, this is just the beginning
            of a block, length may increase later (after current
            method returns).
             */
            curr_rec.len = 1;
        }

        return curr_rec;
    }

    /**
     * iterate over the original sequence, create map file records
     * for each block of sequence (each block consists of identical
     * events, e.g. insertion, or no-change), append records to
     * stringBuilder. For details about map file format, look into
     * the comments below.
     *
     * @param sb output string
     * @param chr_name name of haploid perturbed sequence
     * @param ref_seq name of original sequence
     * @param genome flags for each position of the original sequence
     * @param ins_seq recording positions of all insertions (here insertions are broader than typical definition)
     */
    private void makePosMap(StringBuilder sb, String chr_name, Sequence ref_seq, byte[] genome,
                            Hashtable<Integer, FlexSeq> ins_seq) {

        // host is the perturbed genome
        // ref is b37 or hg19
        // Len is length of the block
        // the positions are the start locations of the block (inclusive)
        // for Insertion and DUP, etc the start location is the location before the event.
        // Feature is annotated in the reads, SEQ, DEL, Insertion, DUP, INV
        // INV is whether the block is inverted
        //bw.write("#Len\tHOST_chr\tHOST_pos\tREF_chr\tREF_pos\tDIRECTION\tFEATURE\tVAR_ID");
        //bw.newLine();
        /* in principle, if a record corresponds to a sequence not existent on host, then
        we use the coordinate of host before that event; similarly, if a record corresponds
        to a sequence not existent on reference (e.g. insertion, duplication), then we
        use the coordinate of reference before that event. In other words, non-existent sequence
        should always point to a locus upstream.
         */

        // TODO deal with var_id

        // iterate through both genomes

        // these are 1-indexed
        HostRefIdx hf_idx = new HostRefIdx();
        hf_idx.host_idx = 1;
        hf_idx.ref_idx = 1;

        MapRecord curr_rec = new_curr_rec(sb, 0, chr_name, ref_seq.getName(), hf_idx, genome, ins_seq);

        for (int idx = 1; idx < genome.length; idx++) {
            // if still in the same block increment the length
            boolean same_block = true;

            switch (curr_rec.feature) {
                case "DEL":
                    if (genome[idx] != DELETED_BASE || ins_seq.containsKey(idx + 1)) {
                        same_block = false;
                    }
                    break;
                case "SEQ":
                    if (genome[idx] == DELETED_BASE || ins_seq.containsKey(idx + 1)) {
                        same_block = false;
                    }
                    break;
                default:
                    same_block = false;
                    break;
            }

            if (same_block) {
                curr_rec.len++;
            } else {
                // otherwise output the block
                adjust_idx(curr_rec, hf_idx);
                sb.append(curr_rec);
                sb.append('\n');

                curr_rec = new_curr_rec(sb, idx, chr_name, ref_seq.getName(), hf_idx, genome, ins_seq);
            }

        }

        // make sure the block is outputted?
        sb.append(curr_rec);
        sb.append('\n');

    }

    /**
     * write perturbed maternal or paternal genomes into files
     * flag array and sequences to be inserted for each haploid
     * genome are required.
     *
     * here, insertion has broader meaning as all variants are
     * considered combinations of insertion + getReferenceAlleleLength
     *
     * instead of calling writeDiploid and then writeMultiploid
     * since we really only have two haplotypes, we can call
     * writeHaploid twice. Got more haplotypes? Use a loop.
     * This avoids creation of anonymous lists.
     * @param sequence
     * @param ins_seq
     * @param sequenceName
     * @param sequenceFileName
     */
    private void writeHaploid(final byte[] sequence, Hashtable<Integer, FlexSeq> ins_seq,
                              final String sequenceName, final String sequenceFileName) {
        try {
            FileWriter fw = new FileWriter(new File(outDir, sequenceFileName));
            BufferedWriter bw = new BufferedWriter(fw);
            writeGenome(bw, sequenceName, sequence, ins_seq);
            bw.close();
            fw.close();
        } catch (IOException ex) {
            log.error(ex.toString());
        }
    }

    /**
     * write specified perturbed haploid genome into output stream
     * with proper line wrapping
     * inserted sequences will be inserted on the fly
     * @param bw
     * @param name
     * @param genome
     * @param ins_seq
     * @throws IOException
     */
    private void writeGenome(BufferedWriter bw, String name, byte[] genome,
                             Hashtable<Integer, FlexSeq> ins_seq) throws IOException {
        //TODO: escalate this to instance variable
        final int line_width = 50; // this is default for FASTA files

        // look-up table for where on the chromosome there are insertions
        boolean[] ins_flag = new boolean[genome.length];
        //TODO: refactor the iteration to enhanced for loop using keySet
        Enumeration<Integer> enm = ins_seq.keys();
        while (enm.hasMoreElements()) {
            Integer key = enm.nextElement();
            ins_flag[key - 1] = true;
        }

        // write header
        bw.write(">" + name);
        bw.newLine();

        StringBuilder line = new StringBuilder();
        for (int p = 0; p < genome.length; p++) {
            if (ins_flag[p]) {
                //so inserted sequences are indexed by 1-based coordinates
                byte[] new_seq = ins_seq.get(p + 1).getSeq();
                if (new_seq != null && new_seq.length > 0) {
                    line.append(new String(new_seq));
                }
            }
            if (genome[p] != DELETED_BASE) {
                line.append((char) genome[p]);
            }
            //TODO: refactor output writing with line wrapping into a method
            while (line.length() >= line_width) {
                bw.write(line.toString(), 0, line_width);
                bw.newLine();
                line.delete(0, line_width);
            }
        }
        while (line.length() > 0) {
            int n = line.length();
            if (line_width < n) {
                n = line_width;
            }
            bw.write(line.toString(), 0, n);
            bw.newLine();
            line.delete(0, n);
        }

    }

    /*
     * Writes out the vcf record for all variants that have added_variants =
     * true
     */
    private void writeVCF(final Sequence ref_seq, final List<Variant> varList,
                          final List<Boolean> paternal_added_variants,
                          final List<Boolean> maternal_added_variants,
                          final boolean output_paternal, final boolean output_maternal) {
        String file_name = ref_seq.getName() + "_" + id + ".vcf";
        List<String> idList = new ArrayList<>();
        idList.add(id);
        log.info("Writing out the true variants for " + ref_seq.getName());
        try {
            FileWriter fw = new FileWriter(new File(outDir, file_name));
            BufferedWriter bw = new BufferedWriter(fw);

            // write header
            bw.write(generateVCFHeader(chrfiles.get(0), idList));

            int num_vars = varList.size();
            //it seems indexMap is not necessary here
            //TODO: remove indexMap
            final Map<Variant, Integer> indexMap = new HashMap<>();
            for (int i = 0; i < num_vars; i++) {
                indexMap.put(varList.get(i), i);
            }
            Collections.sort(varList);

            for (int i_ = 0; i_ < num_vars; i_++) {
                final int i = indexMap.get(varList.get(i_));

                if (!output_maternal && output_paternal && !paternal_added_variants.get(i)) {
                    /*
                    output_maternal==false, output_paternal==true, parternal_added_variants[i] == false
                    paternal was supposed to be output, however, variant was discarded
                     */
                    continue;
                } else if (output_maternal && !output_paternal && !maternal_added_variants.get(i)) {
                    //similar to above, but maternal and paternal are switched
                    continue;
                }

                if (paternal_added_variants.get(i)
                        || maternal_added_variants.get(i)) {

                    // System.err.println("write var: " + i);

                    Variant curr_var = varList.get(i_);
                    curr_var.calculateExtraBase(ref_seq);

                    // chromosome name
                    bw.write(curr_var.getChr().toString());
                    bw.write("\t");
                    // start position
                    bw.write(String.valueOf(curr_var.getPos()
                            - curr_var.getRef_deleted().length()));
                    bw.write("\t");
                    // variant id
                    bw.write(curr_var.getVar_id());
                    bw.write("\t");
                    // ref allele
                    bw.write(curr_var.getOrig_Ref() + curr_var.getExtraBase());
                    bw.write("\t");
                    // alt alleles
                    bw.write(curr_var.alt_string().toString());
                    bw.write("\t");
                    // variant quality
                    bw.write(".\t");
                    // pass label
                    bw.write(curr_var.getFilter());
                    bw.write("\t");
                    // INFO
                    // TODO write len here
                    StringBuilder sbStr = new StringBuilder();
                    if (curr_var.getType() == VariantOverallType.Tandem_Duplication) {
                        sbStr.append("SVTYPE=DUP;");
                        sbStr.append("SVLEN=");
                        sbStr.append(curr_var.getLength());
                    } else if (curr_var.getType() == VariantOverallType.Inversion) {
                        sbStr.append("SVTYPE=INV;");
                        sbStr.append("SVLEN=");
                        sbStr.append(curr_var.getLength());
                    } else {
                        sbStr.append("SVLEN=");
                        sbStr.append(curr_var.getLength());
                    }
                    bw.write(sbStr.toString());
                    bw.write("\t");
                    // label (GT)
                    boolean hasCN = curr_var.hasCN();
                    if (hasCN) {
                        bw.write("GT:CN\t");
                    } else {
                        bw.write("GT\t");
                    }
                    // the genotype
                    // for this one we need to work out which one is added
                    boolean output_both = output_paternal && output_maternal;

                    sbStr = new StringBuilder();
                    if (output_paternal) {
                        if (paternal_added_variants.get(i)) {
                            sbStr.append(String.valueOf(curr_var.paternal()));
                        } else {
                            sbStr.append("0");
                        }
                    }
                    if (output_both) {
                        sbStr.append("|");
                    }
                    if (output_maternal) {
                        if (maternal_added_variants.get(i)) {
                            sbStr.append(String.valueOf(curr_var.maternal()));
                        } else {
                            sbStr.append("0");
                        }
                    }

                    if (hasCN) {
                        sbStr.append(":");
                        if (output_paternal) {
                            sbStr.append(String.valueOf(curr_var.getCN(curr_var
                                    .paternal())));
                        }
                        if (output_both) {
                            sbStr.append("|");
                        }
                        if (output_maternal) {
                            sbStr.append(String.valueOf(curr_var.getCN(curr_var
                                    .maternal())));
                        }
                    }

                    bw.write(sbStr.toString());
                    bw.newLine();
                }
            }

            bw.close();
            fw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    /**
     * generate VCF file header
     * @param referenceFileName reference file name
     * @param sampleNames list of sample names
     * @return
     */
    private String generateVCFHeader(String referenceFileName, List<String> sampleNames) {
        String VCFHeader = "##fileformat=VCFv4.1\n" +
                "##reference=" + referenceFileName + "\n" +
                "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of variant\">\n" +
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##ALT=<ID=DEL,Description=\"Deletion\">\n" +
                "##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">\n" +
                "##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">\n" +
                "##ALT=<ID=DUP,Description=\"Duplication\">\n" +
                "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n" +
                "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n" +
                "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n" +
                "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n" +
                "##ALT=<ID=INV,Description=\"Inversion\">\n" +
                "##ALT=<ID=CNV,Description=\"Copy number variable region\">\n" +
                "##ALT=<ID=ITX,Description=\"Intra-chromosomal translocation\">\n" +
                "##ALT=<ID=CTX,Description=\"Inter-chromosomal translocation\">\n";
        StringJoiner joiner = new StringJoiner("\t");
        for (String id : sampleNames) {
            joiner.add(id);
        }
        return VCFHeader + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" +
                (sampleNames.isEmpty() ? "" : "\t") + joiner.toString() + "\n";
    }

    //"#Len\tHOST_chr\tHOST_pos\tREF_chr\tREF_pos\tDIRECTION\tFEATURE\tVAR_ID"
    // this is more like a struct :)
    // I think MapRecord stands for records in MFF file (the map file)
    class MapRecord {
        public int len = 0;
        public String host_chr = "";
        public int host_pos = 0;
        public String ref_chr = "";
        public int ref_pos = 0;
        public boolean dir = true; // true is forward
        public String feature = ".";
        public String var_id = ".";

        public String toString() {
            StringJoiner joiner = new StringJoiner("\t");
            joiner.add(Integer.toString(len));
            joiner.add(host_chr);
            joiner.add(Integer.toString(host_pos));
            joiner.add(ref_chr);
            joiner.add(Integer.toString(ref_pos));
            if (dir) {
                joiner.add("+");
            } else {
                joiner.add("-");
            }
            joiner.add(feature);
            joiner.add(var_id);

            return joiner.toString();
        }
    }

    /**
     * class for storing host genome and
     * reference genome index
     */
    private class HostRefIdx {
        public int host_idx = 0;
        public int ref_idx = 0;
    }

}
