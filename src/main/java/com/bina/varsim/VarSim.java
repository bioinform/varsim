package com.bina.varsim;

import com.bina.varsim.fastqLiftover.FastqLiftOver;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.Arrays;

/**
 * VarSim
 * <p/>
 * This is the main class the calls each of the individual tools
 * <p/>
 * Created by johnmu on 11/25/14.
 */

// parameter parsing will be updated to be nicer later
public class VarSim {
    private final static Logger log = Logger.getLogger(VarSim.class.getName());

    String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();

    public void run(String[] args) throws IOException {
        String usage = "java -jar VarSim.jar <tool> <tool_args>... \n"
                + "      --= Simulation =-- \n"
                + "       randvcf2vcf    -- Randomly samples variants from a VCF file\n"
                + "       randdgv2vcf    -- Randomly samples variants from a DGV database file\n"
                + "      --= Validation =-- \n"
                + "       vcfcompare     --  Generate JSON describing vcf accuracy relative to truth \n"
                + "       samcompare     --  Generate JSON describing alignment accuracy relative to truth \n"
                + "       vcfstats       --  Generate stats on size range and variant types in a VCF\n"
                + "      --= Internal =-- \n"
                + "       vcf2diploid    -- Enhanced version of vcf2diploid from alleleseq \n"
                + "       fastq_liftover -- Lifts over simulated FASTQ files to reference coordinates \n"
                + "\n";
        if (args.length == 0) {
            System.err.println(VERSION);
            System.err.println(usage);
            System.exit(1);
        }

        String[] pass_args = Arrays.copyOfRange(args, 1, args.length);

        switch (args[0]) {
            case "vcf2diploid":
                new VCF2diploid().run(pass_args);
                break;
            case "randvcf2vcf":
                new RandVCF2VCF().run(pass_args);
                break;
            case "randdgv2vcf":
                new RandDGV2VCF().run(pass_args);
                break;
            case "vcfstats":
                new VCFstats().run(pass_args);
                break;
            case "vcfcompare":
                new VCFcompare().run(pass_args);
                break;
            case "samcompare":
                new SAMcompare().run(pass_args);
                break;
            case "randbed2vcf":
                new RandBED2VCF().run(pass_args);
                break;
            case "fastq_liftover":
                new FastqLiftOver().run(pass_args);
                break;
            case "json_inserter":
                new JSONInserter().run(pass_args);
                break;
            default:
                log.error("Unknown tool: " + args[0]);
                System.err.println(usage);
                System.exit(1);
                break;
        }

    }

    public static void main(String[] args) throws IOException {
        VarSim runner = new VarSim();
        runner.run(args);
    }

}
