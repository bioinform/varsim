package com.bina.varsim;

public enum VarSimToolNamespace {
    VCF2Diploid("vcf2diploid", "Create a diploid genome as associated files from a reference genom and some VCF files."),
    RandVCF2VCF("randvcf2vcf", "Randomly samples variants from VCF file."),
    RandDGV2VCF("randdgv2vcf", "Randomly samples variants from DGV file."),
    RandBED2VCF("randbed2vcf", "Generates a VCF file (to stdout) from an insertion and a deletion BED file. Insertions sequences are randomly sampled from the insert_seq file. This is designed for the Venter SV BED files."),
    RandSequenceVCF("randsequencevcf", "Fill in missing insertion sequences by randomly sampling from a sequence file."),
    VCFStats("vcfstats", "Get stats on variants from a VCF"),
    VCFCompare("vcfcompare", "Generates a JSON with accuracy statistics of a VCF file relative to a truth"),
    SAMCompare("samcompare", "Analyses the accuracy of the alignments in a SAM/BAM file. bed_file restricts the analysis to the bed regions"),
    FastqLiftover("fastq_liftover", "Lift FASTQs to the right reference"),
    LongISLNDLiftover("longislnd_liftover", "Lift read map files to the right reference"),
    JSONInserter("json_inserter", "Inserts n JSON files to one HTML to create n HTML files"),
    Unknown("", null);

    protected String command;
    protected String description;

    VarSimToolNamespace(final String command, final String description) {
        this.command = command;
        this.description = description;
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
            }
        }
        return Unknown;
    }
}
