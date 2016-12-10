package com.bina.varsim.tools;

/**
 *
 */

import com.bina.varsim.VarSimTool;
import com.bina.varsim.types.BedFile;
import com.bina.varsim.types.Parent_record;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.util.VCFparser;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

public class VCFstats extends VarSimTool {
    private final static Logger log = Logger.getLogger(VCFstats.class.getName());

    final String DESCRIPTION = "Compute the distribution of variants in a VCF file";
    BedFile bed_file;

    @Option(name = "-bed", usage = "BED file to restrict the analysis [Optional]", metaVar = "BED_file")
    String bed_filename = null;

    @Option(name = "-vcf", usage = "VCF file to analyse [Required]", metaVar = "VCF_file", required = true)
    String vcf_filename;

    /**
     * @param args command line parameters
     */
    public static void main(String[] args) {
        VCFstats runner = new VCFstats();
        runner.run(args);
    }

    public void run(String[] args) {
        CmdLineParser cmd_parser = new CmdLineParser(this);

        try {
            cmd_parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java -jar vcfstats.jar [options...]");
            // print the list of available options
            cmd_parser.printUsage(System.err);
            System.err.println(DESCRIPTION);
            return;
        }

        if (printVersion) {
            System.out.println(VERSION);
            return;
        }

        bed_file = null;

        if (bed_filename != null) {
            log.info("Reading BED file...");
            bed_file = new BedFile(bed_filename);
        }

        Parent_record data = new Parent_record();

        VCFparser parser = new VCFparser(vcf_filename, null, false);

        while (parser.hasMoreInput()) {
            Variant var = parser.parseLine();
            if (var == null) {
                continue;
            }
            if (var.getGoodMaternal() == 0 && var.getGoodPaternal() == 0) {
                continue;
            }
            data.add(var, bed_file);
        }

        System.out.println(data);
    }

}
