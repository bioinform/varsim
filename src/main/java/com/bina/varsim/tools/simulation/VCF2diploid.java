package com.bina.varsim.tools.simulation;

//--- Java imports ---

import com.bina.varsim.types.*;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.types.variant.VariantOverallType;
import com.bina.varsim.types.variant.VariantType;
import com.bina.varsim.util.SimpleReference;
import com.bina.varsim.util.StringUtilities;
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
    final int LineWidth = 50; // this is default for FASTA files
    static final long SEED_ARG = 3333;
    private final static Logger log = Logger.getLogger(VCF2diploid.class.getName());
    public final static char DELETED_BASE = '~';
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
    public static void main(final String[] args) {
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
    public void run(final String[] args) {
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
        //sort variants of each chromosome by coordinates
        for (ChrString chr : variants.keySet()) {
            Collections.sort(variants.get(chr));
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
                    if (addVariant(paternalMaskedSequence, referenceSequence, allSequences, var, var.paternal(),
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
                    if (addVariant(maternalMaskedSequence, referenceSequence, allSequences, var, var.maternal(),
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
    private boolean addVariant(final byte[] maskedSequence, final Sequence referenceSequence, final SimpleReference allSequences,
                               final Variant variant, final int allele, final Hashtable<Integer, FlexSeq> InsertPosition2Sequence) {

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
            //TODO: wrap insertion generation in Variant class
            //NOTE: only translocation with ACCEPT subtype is marked as Translocation
            //translocation with REJECT subject is marked as Deletion
            if (variantType == VariantType.Translocation) {
              if (Variant.TranslocationSubtype.ACCEPT.equals(variant.getTranslocationSubtype(allele))) {
                  if (variant.getPos2(allele) <= variant.getEnd2(allele)) {
                      insertions = allSequences.getSequence(variant.getChr2(allele)).subSeq(variant.getPos2(allele), variant.getEnd2(allele) + 1);
                  } else {
                      insertions = allSequences.getSequence(variant.getChr2(allele)).revComp(variant.getEnd2(allele), variant.getPos2(allele) + 1);
                  }
              } else {
                  //REJECT, essentially deletion
                  insertions = null;
              }
            }

            if (variantType == VariantType.Inversion) {
                // Treat this as getReferenceAlleleLength then insertion of reverse complement
                //System.err.println("Insert INV");
                referenceAlleleLength = variant.insertionLength(allele);
                insertions = referenceSequence.revComp(position, position + referenceAlleleLength);
            }

            if (variantType == VariantType.Tandem_Duplication) {
                // treat this as getReferenceAlleleLength then insertion of several copies
                // this prevents the original sequence to be altered and
                // become less like a duplication
                // TODO make length correct
                // System.err.println("Insert DUP");

                referenceAlleleLength = variant.insertionLength(allele);
                int single_ins_len = variant.insertionLength(allele);
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
                FlexSeq s = null;

                if (FlexSeq.Type.TRA.equals(variant.getAlt(allele).getSeqType())) {
                    s = new FlexSeq.Builder().sequence(insertions).
                            type(variant.getAlt(allele).getSeqType()).
                            copyNumber(variant.getAlt(allele).getCopyNumber()).
                            length(variant.getAlt(allele).length()).
                            variantId(variant.getVariantId()).
                            chr2(variant.getChr2(allele)).
                            pos2(variant.getPos2(allele)).
                            end2(variant.getEnd2(allele)).
                            referenceAlleleLength(variant.getReferenceAlleleLength()).build();
                } else {
                    s = new FlexSeq.Builder().sequence(insertions).
                            type(variant.getAlt(allele).getSeqType()).
                            copyNumber(variant.getAlt(allele).getCopyNumber()).
                            length(variant.getAlt(allele).length()).
                            variantId(variant.getVariantId()).build();
                }
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
    private void makePosMap(final StringBuilder sb, final String chr_name, final Sequence ref_seq, final byte[] genome,
                            final Hashtable<Integer, FlexSeq> ins_seq) {

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

        // TODO deal with varId

        // iterate through both genomes

        // these are 1-indexed
        HostRefIdx hostRefIdx = new HostRefIdx();
        hostRefIdx.hostIdx = 1;
        hostRefIdx.refIdx = 1;

        MapRecord currentMapRecord = MapRecord.generateNewMapRecord(sb, 0, chr_name, ref_seq.getName(), hostRefIdx, genome, ins_seq);

        for (int idx = 1; idx < genome.length; idx++) {
            // if still in the same block increment the length
            boolean same_block = true;

            switch (currentMapRecord.feature) {
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
                currentMapRecord.len++;
            } else {
                // otherwise output the block
                hostRefIdx.adjust_idx(currentMapRecord);
                sb.append(currentMapRecord);
                sb.append('\n');

                currentMapRecord = MapRecord.generateNewMapRecord(sb, idx, chr_name, ref_seq.getName(), hostRefIdx, genome, ins_seq);
            }
        }

        // make sure the block is outputted?
        sb.append(currentMapRecord);
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
     * @param insertPosition2Sequence
     * @throws IOException
     */
    private void writeGenome(final BufferedWriter bw, final String name, final byte[] genome,
                             final Hashtable<Integer, FlexSeq> insertPosition2Sequence) throws IOException {

        // write header
        bw.write(">" + name);
        bw.newLine();

        StringBuilder line = new StringBuilder();
        for (int p = 0; p < genome.length; p++) {
            if (insertPosition2Sequence.containsKey(p + 1)) {
                //so inserted sequences are indexed by 1-based coordinates
                byte[] new_seq = insertPosition2Sequence.get(p + 1).getSequence();
                if (new_seq != null && new_seq.length > 0) {
                    line.append(new String(new_seq));
                }
            }
            if (genome[p] != DELETED_BASE) {
                line.append((char) genome[p]);
            }
            while (line.length() >= LineWidth) {
                bw.write(line.toString(), 0, LineWidth);
                bw.newLine();
                line.delete(0, LineWidth);
            }
        }
        while (line.length() > 0) {
            int n = line.length();
            if (LineWidth < n) {
                n = LineWidth;
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

            for (int i = 0; i < num_vars; i++) {
                if (!output_maternal && output_paternal && !paternal_added_variants.get(i)) {
                    /*
                    output_maternal==false, output_paternal==true, parternal_added_variants[i] == false
                    paternal was supposed to be output, however, variant was discarded due to ,e.g.,
                    overlapping with other variants.
                     */
                    continue;
                } else if (output_maternal && !output_paternal && !maternal_added_variants.get(i)) {
                    //similar to above, but maternal and paternal are switched
                    continue;
                }

                if (paternal_added_variants.get(i)
                        || maternal_added_variants.get(i)) {

                    Variant curr_var = varList.get(i);
                    curr_var.calculateExtraBase(ref_seq);

                    // chromosome name
                    bw.write(curr_var.getChr().toString());
                    bw.write("\t");
                    // start position
                    bw.write(String.valueOf(curr_var.getPos()
                            - curr_var.getRef_deleted().length));
                    bw.write("\t");
                    // variant id
                    bw.write(curr_var.getVariantId());
                    bw.write("\t");
                    // ref allele
                    bw.write(curr_var.getReferenceString() + curr_var.getExtraBase());
                    bw.write("\t");
                    // alt alleles
                    bw.write(curr_var.alternativeAlleleString());
                    bw.write("\t");
                    // variant quality
                    bw.write(".\t");
                    // pass label
                    bw.write(curr_var.getFilter());
                    bw.write("\t");
                    // INFO
                    // TODO write len here
                    //TODO: use StringJoiner to replace StringBuilder
                    StringBuilder sbStr = new StringBuilder();
                    if (curr_var.getType() == VariantOverallType.Tandem_Duplication) {
                        sbStr.append("SVTYPE=DUP;");
                        sbStr.append("SVLEN=");
                        sbStr.append(curr_var.getLength());
                    } else if (curr_var.getType() == VariantOverallType.Inversion) {
                        sbStr.append("SVTYPE=INV;");
                        sbStr.append("SVLEN=");
                        sbStr.append(curr_var.getLength());
                    } else if (curr_var.getType() == VariantOverallType.Translocation) {
                        sbStr.append("SVTYPE=TRA;");
                        sbStr.append("SVLEN=");
                        sbStr.append(curr_var.getLength());
                        sbStr.append(";");
                        sbStr.append("END=");
                        sbStr.append(Integer.toString(curr_var.getEnd()));
                        sbStr.append(";");
                        sbStr.append("TRASUBTYPE=");
                        sbStr.append(StringUtilities.concatenateArray(curr_var.getAllTranslocationSubtypeString(), ","));
                        sbStr.append(";");
                        //chr2,pos2,end2
                        sbStr.append("CHR2=");
                        sbStr.append(StringUtilities.concatenateArray(curr_var.getAllChr2(), ","));
                                sbStr.append(";");
                        sbStr.append("POS2=");
                        sbStr.append(StringUtilities.concatenateArray(curr_var.getAllPos2(), ","));
                        sbStr.append(";");
                        sbStr.append("END2=");
                        sbStr.append(StringUtilities.concatenateArray(curr_var.getAllEnd2(), ","));
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
     *
     * a little explanation about translocation-related meta-info lines:
     * each locus involved in translocations is modeled independently as
     * a cut-paste event with slight variations for different types of
     * translocations. Basically, a region at locus A (the sink) will be
     * cut (deleted), and a region at locus B (the source) will be placed
     * at the sink. The placement may be: complete transfer; complete transfer
     * with inversion; no transfer (one-way or unbalanced translocation).
     *
     * @param referenceFileName reference file name
     * @param sampleNames list of sample names
     * @return
     */
    private String generateVCFHeader(final String referenceFileName, final List<String> sampleNames) {
        String VCFHeader = "##fileformat=VCFv4.1\n" +
                "##reference=" + referenceFileName + "\n" +
                /*
                SVLEN is for alternative allele in truth VCF
                 */
                "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of variant\">\n" +
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n" +
                /*if POS2<=END2, then another sequence is inserted at positive strand
                if POS2>=END2, then reversed sequence is inserted at negative strand (insert with inversion)
                 */
                "##INFO=<ID=POS2,Number=A,Type=Integer,Description=\"1-based Start position of source sequence\">\n" +
                "##INFO=<ID=END2,Number=A,Type=Integer,Description=\"1-based End position of source sequence\">\n" +
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"1-based End position of sink sequence\">\n" +
                "##INFO=<ID=CHR2,Number=A,Type=String,Description=\"Chromosome of source sequence\">\n" +
                "##INFO=<ID=TRASUBTYPE,Number=A,Type=String,Description=\"Subtype of translocation event:" +
                " source sequence deleted (REJECT); source sequence accepted (ACCEPT).\">\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                /*CN is defined as Integer in VCF4.1,4.3, making it impossible to specify multiple CN values
                here we changed it to String to allow such behavior.
                 */
                "##FORMAT=<ID=CN,Number=1,Type=String,Description=\"Copy number genotype. Different copy numbers do not imply same genotypes.\">\n" +
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
                "##ALT=<ID=TRA,Description=\"Translocation\">\n" +
                "##ALT=<ID=CTX,Description=\"Inter-chromosomal translocation\">\n";
        StringUtilities.StringJoiner joiner = new StringUtilities.StringJoiner("\t");
        for (String id : sampleNames) {
            joiner.add(id);
        }
        return VCFHeader + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" +
                (sampleNames.isEmpty() ? "" : "\t") + joiner.toString() + "\n";
    }
}
