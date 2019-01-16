package com.bina.varsim;

public enum VarSimToolNamespace {
    VCF2Diploid("vcf2diploid", "Create a diploid genome as associated files from a reference genom and some VCF files.", com.bina.varsim.tools.simulation.VCF2diploid.class),
    RandVCF2VCF("randvcf2vcf", "Randomly samples variants from VCF file.", com.bina.varsim.tools.simulation.RandVCF2VCF.class),
    RandDGV2VCF("randdgv2vcf", "Randomly samples variants from DGV file.", com.bina.varsim.tools.simulation.RandDGV2VCF.class),
    RandBED2VCF("randbed2vcf", "Generates a VCF file (to stdout) from an insertion and a deletion BED file. Insertions sequences are randomly sampled from the insert_seq file. This is designed for the Venter SV BED files.", com.bina.varsim.tools.simulation.RandBED2VCF.class),
    RandSequenceVCF("randsequencevcf", "Fill in missing insertion sequences by randomly sampling from a sequence file.", com.bina.varsim.tools.simulation.RandSequenceVCF.class),
    VCFStats("vcfstats", "Get stats on variants from a VCF", com.bina.varsim.tools.VCFstats.class),
    VCFCompare("vcfcompare", "Generates a JSON with accuracy statistics of a VCF file relative to a truth", com.bina.varsim.tools.evaluation.VCFcompare.class),
    VCFCompareResultsParser("vcfcompareresultsparser", "Generates a JSON with accuracy statistics given TP, FN, FP VCFs", com.bina.varsim.tools.evaluation.VCFCompareResultsParser.class),
    SAMCompare("samcompare", "Analyses the accuracy of the alignments in a SAM/BAM file. bed_file restricts the analysis to the bed regions", com.bina.varsim.tools.evaluation.SAMcompare.class),
    FastqLiftover("fastq_liftover", "Lift FASTQs to the right reference", com.bina.varsim.fastqLiftover.FastqLiftOver.class),
    LongISLNDLiftover("longislnd_liftover", "Lift read map files to the right reference", com.bina.varsim.fastqLiftover.LongISLNDReadMapLiftOver.class),
    JSONInserter("json_inserter", "Inserts n JSON files to one HTML to create n HTML files", com.bina.varsim.tools.evaluation.JSONInserter.class),
    LiftOver("liftover", "Lift over a file to desired coordinates", com.bina.varsim.tools.LiftOver.class),
    Help("-help", null, null, new String[]{"-h"}),
    Version("-version", null, null),
    Unknown("", null, null);

    public String command;
    public String description;
    Class toolClass;
    String[] commandAliases;

    VarSimToolNamespace(final String command, final String description, final Class toolClass) {
        this.command = command;
        this.description = description;
        this.toolClass = toolClass;
        this.commandAliases = null;
    }

    VarSimToolNamespace(final String command, final String description, final Class toolClass, final String[] aliases) {
        this.command = command;
        this.description = description;
        this.toolClass = toolClass;
        this.commandAliases = aliases;
    }

    public static VarSimToolNamespace fromName(final String name) {
        if (name != null) {
            for (final VarSimToolNamespace varSimTool : values()) {
                if (varSimTool == Unknown) {
                    continue;
                }
                if (varSimTool.command.equals(name) || varSimTool.name().equalsIgnoreCase(name)) {
                    return varSimTool;
                }
                if (varSimTool.commandAliases != null) {
                    for (final String commandAlias : varSimTool.commandAliases) {
                        if (commandAlias.equalsIgnoreCase(name)) {
                            return varSimTool;
                        }
                    }
                }
            }
        }
        return Unknown;
    }
}
