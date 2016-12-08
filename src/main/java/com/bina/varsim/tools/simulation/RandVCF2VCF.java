package com.bina.varsim.tools.simulation;

import com.bina.varsim.constants.Constant;
import com.bina.varsim.types.*;
import com.bina.varsim.types.variant.Variant;
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
 * TODO ignores the input genotypes for now
 */

/**
 * @author johnmu
 */


public class RandVCF2VCF extends RandVCFgenerator {
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
    int numSNP = NUM_SNP_ARG;
    @Option(name = "-num_ins", usage = "Number of simple insertions to sample [" + NUM_INS_ARG + "]")
    int numIns = NUM_INS_ARG;
    @Option(name = "-num_del", usage = "Number of simple deletions to sample [" + NUM_DEL_ARG + "]")
    int numDel = NUM_DEL_ARG;
    @Option(name = "-num_mnp", usage = "Number of MNPs to sample [" + NUM_MNP_ARG + "]")
    int numMNP = NUM_MNP_ARG;
    @Option(name = "-num_complex", usage = "Number of complex variants (other ones) to sample [" + NUM_COMPLEX_ARG + "]")
    int numComplex = NUM_COMPLEX_ARG;
    @Option(name = "-novel", usage = "Average ratio of novel variants[" + NOVEL_RATIO_ARG + "]")
    double ratioNovel = NOVEL_RATIO_ARG;
    @Option(name = "-min_len", usage = "Minimum variant length [" + MIN_LEN_ARG + "], inclusive")
    int minLengthLim = MIN_LEN_ARG;
    @Option(name = "-max_len", usage = "Maximum variant length [" + MAX_LEN_ARG + "], inclusive")
    int maxLengthLim = MAX_LEN_ARG;
    @Option(name = "-prop_het", usage = "Average ratio of novel variants[" + PROP_HET_ARG + "]")
    double propHet = PROP_HET_ARG;

    @Option(name = "-ref", usage = "Reference Genome [Required]", metaVar = "file", required = true)
    String referenceFilename;

    @Option(name = "-vcf", usage = "Known VCF file, eg. dbSNP [Required]", metaVar = "file", required = true)
    String vcfFilename;

    @Option(name = "-t", usage = "Gender of individual [MALE]")
    GenderType gender = GenderType.MALE;

    @Option(name = "-out_vcf", usage = "Output VCF to generate [stdout]")
    String outFilename = null;

    private final Set<VariantType> variantTypes = EnumSet.of(VariantType.SNP, VariantType.Insertion,
            VariantType.Deletion, VariantType.MNP, VariantType.Complex);

    int numNovelAdded;

    public RandVCF2VCF() {
        super();
        numNovelAdded = 0;
    }

    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        // TODO Auto-generated method stub
        RandVCF2VCF runner = new RandVCF2VCF();
        runner.run(args);
    }

    public Map<VariantType, Integer> countVariants(final String vcfFilename, final List<Genotypes> selectedGeno) {
        int totalNumOther = 0;
        int totalNum = 0;

        final Map<VariantType, Integer> variantTypeCounts = new EnumMap(VariantType.class);
        for (final VariantType variantType : variantTypes) {
            variantTypeCounts.put(variantType, 0);
        }

        VCFparser vcfParser = new VCFparser(vcfFilename, null, false, rand);
        Variant prevVar = new Variant(rand);

        // read though once to count the totals, this is so we don't have
        // to store an array of the variants for sampling without replacement
        while (vcfParser.hasMoreInput()) {
            Variant var = vcfParser.parseLine();
            if (var == null) {
                continue;
            }

            // select genotypes here
            ChrString chr = var.getChr();
            int numAlt = var.getNumberOfAlternativeAlleles();

            Genotypes geno = new Genotypes(chr, gender, numAlt, rand, propHet);
            selectedGeno.add(geno);

            if (prevVar.equals(var)) {
                continue;
            }

            // this is ok because var is not changed later
            prevVar = var;

            if (var.maxLen() > maxLengthLim
                    || var.minLen() < minLengthLim) {
                continue;
            }

            int numIters = (geno.geno[0] == geno.geno[1]) ? 1 : 2;
            for (int i = 0; i < numIters; i++) {
                final VariantType variantType = var.getType(geno.geno[i]);
                if (variantTypeCounts.containsKey(variantType)) {
                    variantTypeCounts.put(variantType, variantTypeCounts.get(variantType) + 1);
                } else {
                    totalNumOther++;
                }
            }
            totalNum++;
        }

        log.info("total_num_SNP: " + variantTypeCounts.get(VariantType.SNP));
        log.info("total_num_INS: " + variantTypeCounts.get(VariantType.Insertion));
        log.info("total_num_DEL: " + variantTypeCounts.get(VariantType.Deletion));
        log.info("total_num_MNP: " + variantTypeCounts.get(VariantType.MNP));
        log.info("total_num_COMPLEX: " + variantTypeCounts.get(VariantType.Complex));
        log.info("total_num_skipped: " + totalNumOther);
        log.info("total_num: " + totalNum);

        return variantTypeCounts;
    }

    void sampleFromVCF(final String samplingVCF, final SimpleReference reference, final List<Genotypes> selectedGeno,
                       final Map<VariantType, Integer> variantCounts, final OutputStream outputStream)
            throws IOException {
        // read through VCF file again and output the new sampled VCF file
        log.info("Writing sampled variant file");

        final String VCF_HEADER = "##fileformat=VCFv4.0\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsv\n";

        BufferedWriter out = new BufferedWriter(new OutputStreamWriter(outputStream));
        out.write(VCF_HEADER);

        final Map<VariantType, SampleParams> variantTypeParams = new EnumMap(VariantType.class);
        for (final VariantType variantType : variantTypes) {
            variantTypeParams.put(variantType, new SampleParams());
        }

        final Map<VariantType, Integer> maxCounts = new EnumMap(VariantType.class);
        maxCounts.put(VariantType.SNP, numSNP);
        maxCounts.put(VariantType.Insertion, numIns);
        maxCounts.put(VariantType.Deletion, numDel);
        maxCounts.put(VariantType.MNP, numMNP);
        maxCounts.put(VariantType.Complex, numComplex);

        int genoIdx = 0;

        final VCFparser vcfParser = new VCFparser(samplingVCF, null, false, rand);
        Variant prevVar = new Variant(rand);

        // Read through it a second time, this time we do the sampling
        while (vcfParser.hasMoreInput()) {
            Variant var = vcfParser.parseLine();
            if (var == null) {
                continue;
            }

            Genotypes geno = selectedGeno.get(genoIdx);
            genoIdx++;

            if (prevVar.equals(var)) {
                // duplicate
                continue;
            }

            prevVar = var;

            if (var.maxLen() > maxLengthLim
                    || var.minLen() < minLengthLim) {
                continue;
            }

            // this samples genotypes from each allele individually
            int numIter = (geno.geno[0] == geno.geno[1]) ? 1 : 2;
            for (int i = 0; i < numIter; i++) {
                final VariantType variantType = var.getType(geno.geno[i]);
                if (variantTypeParams.containsKey(variantType)) {
                    geno.geno[i] = sampleGenotype(geno.geno[i], variantTypeParams.get(variantType),
                            maxCounts.get(variantType), variantCounts.get(variantType));
                }
            }
            if (numIter == 1) {
                geno.geno[1] = geno.geno[0];
            }

            // write out variant
            if (geno.isNonRef()) {
                randOutputVcfRecord(out, var, reference, ratioNovel, geno);
            }
        }
        out.close();
    }

    // outputs VCF record with random phase
    void randOutputVcfRecord(BufferedWriter bw, Variant var,
                             SimpleReference ref, double ratioNovel, Genotypes geno)
            throws IOException {

        ChrString chr = var.getChr();

        // determine whether this one is novel
        double randNum = rand.nextDouble();
        if (randNum <= ratioNovel) {
            // make the variant novel, simply modify it
            // TODO maybe modifying it is bad

            int chrLen = ref.getRefLen(chr);
            int buffer = Math.max(
                    10,
                    Math.max(var.maxLen(geno.geno[0]),
                            var.maxLen(geno.geno[1])));
            int startVal = Math.min(buffer, Math.max(chrLen - buffer, 0));
            int endVal = Math.max(chrLen - buffer, Math.min(buffer, chrLen));

            int timeOut = 0;
            int randPos = rand.nextInt(endVal - startVal + 1) + startVal + 1;
            while (!var.setNovelPosition(randPos, ref)) {
                if (timeOut > 100) {
                    log.warn("Error: cannot set novel position: " + (endVal - startVal + 1));
                    log.warn(var.getReferenceAlleleLength());
                    log.warn(var);
                    break;
                }

                randPos = rand.nextInt(endVal - startVal + 1) + startVal + 1;
                //log.info(time_out + " : " + var.getReferenceAlleleLength());

                timeOut++;
            }
            numNovelAdded++;

            var.setVarID("Novel_" + numNovelAdded);
        }

        outputVcfRecord(bw, var, geno.geno[0], geno.geno[1]);
    }

    public void run(String[] args) throws IOException {
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

        if (ratioNovel > 1 || ratioNovel < 0) {
            System.err.println("Novel ratio out of range [0,1]");
            System.exit(1);
        }

        final SimpleReference reference = new SimpleReference(referenceFilename);

        rand = new Random(seed);

        final List<Genotypes> selectedGeno = new ArrayList<>();
        final Map<VariantType, Integer> variantTypeCounts = countVariants(vcfFilename, selectedGeno);

        final OutputStream outputStream = (outFilename != null) ? new FileOutputStream(outFilename) : System.out;
        sampleFromVCF(vcfFilename, reference, selectedGeno, variantTypeCounts, outputStream);
    }

}
