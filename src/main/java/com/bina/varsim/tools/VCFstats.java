package com.bina.varsim.tools;

/**
 *
 */

import com.bina.varsim.VarSimTool;
import com.bina.varsim.VarSimToolNamespace;
import com.bina.varsim.types.BedFile;
import com.bina.varsim.types.ParentRecord;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.util.VCFparser;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Option;

public class VCFstats extends VarSimTool {
    private final static Logger log = Logger.getLogger(VCFstats.class.getName());

    BedFile bed_file;

    @Option(name = "-bed", usage = "BED file to restrict the analysis [Optional]", metaVar = "BED_file")
    String bed_filename = null;

    @Option(name = "-vcf", usage = "VCF file to analyse [Required]", metaVar = "VCF_file", required = true)
    String vcf_filename;

    public VCFstats(final String command, final String description) {
        super(command, description);
    }

    /**
     * @param args command line parameters
     */
    public static void main(String[] args) {
        new VCFstats("", VarSimToolNamespace.VCFStats.description).run(args);
    }

    public void run(String[] args) {
        if (!parseArguments(args)) {
            return;
        }

        bed_file = null;

        if (bed_filename != null) {
            log.info("Reading BED file...");
            bed_file = new BedFile(bed_filename);
        }

        ParentRecord data = new ParentRecord();

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
