package com.bina.varsim;

public enum VarSimToolNamespace {
    VCF2Diploid("vcf2diploid", com.bina.varsim.tools.simulation.VCF2diploid.class),
    RandVCF2VCF("randvcf2vcf", com.bina.varsim.tools.simulation.RandVCF2VCF.class),
    RandDGV2VCF("randdgv2vcf", com.bina.varsim.tools.simulation.RandDGV2VCF.class),
    RandBED2VCF("randbed2vcf", com.bina.varsim.tools.simulation.RandBED2VCF.class),
    RandSequenceVCF("randsequencevcf", com.bina.varsim.tools.simulation.RandSequenceVCF.class),
    VCFStats("vcfstats", com.bina.varsim.tools.VCFstats.class),
    VCFCompare("vcfcompare", com.bina.varsim.tools.evaluation.VCFcompare.class),
    SAMCompare("samcompare", com.bina.varsim.tools.evaluation.SAMcompare.class),
    FastqLiftover("fastq_liftover", com.bina.varsim.fastqLiftover.FastqLiftOver.class),
    LongISLNDLiftover("longislnd_liftover", com.bina.varsim.fastqLiftover.LongISLNDReadMapLiftOver.class),
    JSONInserter("json_inserter", com.bina.varsim.tools.evaluation.JSONInserter.class),
    Unknown("", null);

    protected String toolName;
    Class toolClass;

    VarSimToolNamespace(final String toolName, Class toolClass) {
        this.toolName = toolName;
        this.toolClass = toolClass;
    }

    public static VarSimToolNamespace fromName(final String name) {
        if (name != null) {
            for (final VarSimToolNamespace varSimTool : values()) {
                if (varSimTool == Unknown) {
                    continue;
                }
                if (varSimTool.toolName.equals(name) || varSimTool.name().equalsIgnoreCase(name)) {
                    return varSimTool;
                }
            }
        }
        return Unknown;
    }
}
