package com.bina.varsim;

import com.bina.varsim.fastqLiftover.FastqLiftOver;
import com.bina.varsim.fastqLiftover.LongISLNDReadMapLiftOver;
import com.bina.varsim.tools.VCFstats;
import com.bina.varsim.tools.evaluation.JSONInserter;
import com.bina.varsim.tools.evaluation.SAMcompare;
import com.bina.varsim.tools.evaluation.VCFcompare;
import com.bina.varsim.tools.evaluation.VCFCompareResultsParser;
import com.bina.varsim.tools.simulation.*;
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

    final String VERSION = getClass().getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException {
        new VarSim().run(args);
    }

    private void printUsage() {
        final String usage = "java -jar VarSim.jar <tool> <tool_args>... \n"
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

        System.err.println(VarSim.class.getSimpleName() + " " + VERSION);
        System.err.println(usage);
    }

    void runVarSimTool(final VarSimToolNamespace tool, final String[] args) throws IOException {
        String[] pass_args = Arrays.copyOfRange(args, 1, args.length);

        final String command = tool.command;
        final String description = tool.description;

        switch (tool) {
            case VCF2Diploid:
                new VCF2diploid(command, description).run(pass_args);
                break;
            case RandVCF2VCF:
                new RandVCF2VCF(command, description).run(pass_args);
                break;
            case RandDGV2VCF:
                new RandDGV2VCF(command, description).run(pass_args);
                break;
            case VCFStats:
                new VCFstats(command, description).run(pass_args);
                break;
            case VCFCompare:
                new VCFcompare(command, description).run(pass_args);
                break;
            case VCFCompareResultsParser:
                new VCFCompareResultsParser(command, description).run(pass_args);
                break;
            case SAMCompare:
                new SAMcompare(command, description).run(pass_args);
                break;
            case RandBED2VCF:
                new RandBED2VCF(command, description).run(pass_args);
                break;
            case FastqLiftover:
                new FastqLiftOver(command, description).run(pass_args);
                break;
            case JSONInserter:
                new JSONInserter(command, description).run(pass_args);
                break;
            case LongISLNDLiftover:
                new LongISLNDReadMapLiftOver(command, description).run(pass_args);
                break;
            case RandSequenceVCF:
                new RandSequenceVCF(command, description).run(pass_args);
                break;
            case Help:
                printUsage();
                break;
            case Version:
                System.out.println(VERSION);
                break;
            default:
                log.error("Unknown tool: " + args[0]);
                System.exit(1);
                break;
        }
    }


    public void run(String[] args) throws IOException {
        if (args.length == 0) {
            printUsage();
            System.exit(1);
        }

        final VarSimToolNamespace varSimToolName = VarSimToolNamespace.fromName(args[0]);
        runVarSimTool(varSimToolName, args);

        return;
    }
}
