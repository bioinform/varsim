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
    @Option(name = "-seed", usage = "Seed for random sampling [" + SEED_DEFAULT + "]")
    static int seed = SEED_DEFAULT;

    @Option(name = "-in_vcf", usage = "Input VCF to fill in", required = true)
    String inFilename = null;

    @Option(name = "-seq", usage = "Sequence file", required=true)
    String sequenceFilename = null;

    @Option(name = "-out_vcf", usage = "Output VCF to generate [stdout]")
    String outFilename = null;

    public RandSequenceVCF() {
        super();
    }

    public RandSequenceVCF(long seed) {
        super(seed);
    }

    public static void main(String[] args) throws IOException {
        // TODO Auto-generated method stub
        RandSequenceVCF runner = new RandSequenceVCF();
        runner.run(args);
    }

    public void run(String[] args) throws IOException {
        String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
        String usage = "Outputs VCF to stdout. Randomly fills in missing ALT sequences.\n";

        CmdLineParser parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java -jar randdgv2vcf.jar [options...]");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }

        rand = new Random(seed);

        final byte[] samplingSequence = fileToByteArray(new File(sequenceFilename));

        final VCFparser vcfParser = new VCFparser(inFilename, false);
        final OutputStream outputStream = (outFilename != null) ? new FileOutputStream(outFilename) : System.out;
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
