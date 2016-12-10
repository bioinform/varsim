package com.bina.varsim.tools.simulation;

import com.bina.varsim.types.Genotypes;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.util.VCFparser;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.Random;

public class RandSequenceVCF extends RandVCFgenerator {
    private final static Logger log = Logger.getLogger(RandSequenceVCF.class.getName());

    static final int SEED_DEFAULT = 3333;
    @Option(name = "-seed", usage = "Seed for random sampling")
    static int seed = SEED_DEFAULT;

    @Option(name = "-in_vcf", usage = "Input VCF to fill in", required = true)
    File inFile = null;

    @Option(name = "-seq", usage = "Sequence file", required=true)
    File sequenceFile;

    @Option(name = "-out_vcf", usage = "Output VCF to generate")
    File outFile = null;

    @Option(name = "-h", usage = "Print help message", help=true, aliases = {"-help"})
    boolean help = false;

    public RandSequenceVCF() {
        super();
    }

    public RandSequenceVCF(long seed) {
        super(seed);
    }

    public static void main(String[] args) throws IOException {
        new RandSequenceVCF().run(args);
    }

    void printUsage(final CmdLineParser parser) {
        final String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
        final String usage = "Outputs VCF to stdout. Randomly fills in missing ALT sequences.\n";

        System.err.println(VERSION);
        System.err.println("java -jar VarSim.jar randsequencevcf [options...]");
        System.err.println(usage);
        // print the list of available options
        parser.printUsage(System.err);
    }

    public void run(String[] args) throws IOException {

        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            printUsage(parser);
            return;
        }

        if (help) {
            printUsage(parser);
            return;
        }

        if (printVersion) {
            System.out.println(VERSION);
            return;
        }

        rand = new Random(seed);

        final byte[] samplingSequence = fileToByteArray(sequenceFile);

        final VCFparser vcfParser = new VCFparser(inFile, false);
        final OutputStream outputStream = (outFile != null) ? new FileOutputStream(outFile) : System.out;
        final BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(outputStream));

        while (vcfParser.hasMoreInput()) {
            final Variant var = vcfParser.parseLine();
            if (var == null) {
                continue;
            }
            final Genotypes geno = var.getGenotypes();
            fillInSeq(var, samplingSequence, geno.geno[0]);
            fillInSeq(var, samplingSequence, geno.geno[1]);

            outputVcfRecord(bw, var, geno.geno[0], geno.geno[1]);
        }
        bw.close();
    }
}
