package com.bina.varsim.tools.simulation;

import com.bina.varsim.VarSimToolNamespace;
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

    @Option(name = "-in_vcf", usage = "Input VCF to fill in", required = true)
    File inFile = null;

    @Option(name = "-seq", usage = "Sequence file", required=true)
    File sequenceFile;

    @Option(name = "-out_vcf", usage = "Output VCF to generate")
    File outFile = null;

    public RandSequenceVCF(final String command, final String description) {
        super(command, description);
    }

    public static void main(String[] args) throws IOException {
        new RandSequenceVCF("", VarSimToolNamespace.RandSequenceVCF.description).run(args);
    }

    public void run(String[] args) throws IOException {
        if (!parseArguments(args)) {
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
