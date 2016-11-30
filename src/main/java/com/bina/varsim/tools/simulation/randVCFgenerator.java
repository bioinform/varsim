package com.bina.varsim.tools.simulation;

import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.variant.VariantOverallType;
import com.bina.varsim.types.SampleParams;
import com.bina.varsim.types.variant.Variant;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Random;


abstract public class RandVCFgenerator {

    protected Random rand;

    /**
     * This sets a default seed, ideally different between different runs
     */
    public RandVCFgenerator() {
        rand = new Random(3333);
    }

    /**
     * @param seed User specified seed for random numbers
     */
    public RandVCFgenerator(long seed) {
        rand = new Random(seed);
    }


    /**
     * This is for sampling without replacement. If the genotype is not sampled, it is replaced with zero
     *
     * @param geno       previous genotype
     * @param seenAdded params for sampling
     * @param numSample number that we want to sample
     * @param numTotal  total number of such variants
     * @param outputAll Don't sample, just output everything
     * @return
     */
    public byte sampleGenotype(byte geno, SampleParams seenAdded, int numSample,
                               int numTotal, boolean outputAll) {
        if (!outputAll) {
            if (seenAdded.addedNum < numSample) {
                double randNum = rand.nextDouble();
                if ((numTotal - seenAdded.seenNum) * randNum >= (numSample - seenAdded.addedNum)) {
                    seenAdded.seenNum++;
                    return 0;
                } else {
                    seenAdded.seenNum++;
                    seenAdded.addedNum++;
                    return geno;
                }
            }
        } else {
            return geno;
        }
        return 0;
    }


    /**
     * This is for sampling without replacement. If the genotype is not sampled, it is replaced with zero
     *
     * @param geno       previous genotype
     * @param seenAdded params for sampling
     * @param numSample number that we want to sample
     * @param numTotal  total number of such variants
     * @return
     */
    public byte sampleGenotype(byte geno, SampleParams seenAdded, int numSample,
                               int numTotal) {
        return sampleGenotype(geno, seenAdded, numSample, numTotal, false);
    }

    /**
     * Fills in the insertions sequence with random sequence sampled from known insertions
     *
     * @param var        Variant to be filled in
     * @param insertSeq Sequence on known insertions all concatented together
     * @param geno       allele to be filled in
     */
    public void fillInSeq(Variant var, byte[] insertSeq, int geno) {
        FlexSeq alt = var.getAlt(geno);
        if (alt != null) {
            if (alt.getType() == FlexSeq.Type.INS) {
                final double NSEGMENTS = 10.0;
                // if insertion sequence is not given, we fill it in
                int len = alt.length();
                byte newSeq[] = new byte[len];
                // randomly duplicate insertion sequence if it is small
                final int segLen = (len > insertSeq.length) ? (int) Math.ceil(insertSeq.length / NSEGMENTS) : len;
                for (int i = 0; i < len; i += segLen) {
                    // choose random start loc
                    int randStart = rand.nextInt(insertSeq.length - segLen);
                    System.arraycopy(insertSeq, randStart, newSeq, i, Math.min(segLen, len - i));
                }
                alt = new FlexSeq(newSeq);
            }

            var.setAlt(geno, alt);
        }
    }


    /**
     * Outputs a VCF record of the variant, supports structural variations
     *
     * @param bw    BufferedWriter for writing
     * @param var   The variant to output
     * @param geno0 First allele
     * @param geno1 Second allele
     * @throws IOException
     */
    public void outputVcfRecord(BufferedWriter bw, Variant var, int geno0, int geno1)
            throws IOException {

        // ignore ACTGN
        String ref = var.getReferenceString().toUpperCase();
        String alt = var.alternativeAlleleString().toUpperCase();

        if (!ref.matches("[ACTGN]*") || (!var.isAltACTGN())) {
            return; // don't output if it is not ACTGN
        }

        // chromosome name
        bw.write(var.getChr().toString());
        bw.write("\t");
        // start position
        bw.write(String.valueOf(var.getPos() - var.getRef_deleted().length()));
        bw.write("\t");
        // variant id
        bw.write(var.getVariantId());
        bw.write("\t");
        // ref allele
        bw.write(ref);
        bw.write("\t");
        // alt alleles
        bw.write(alt);
        bw.write("\t");
        // variant quality
        bw.write(".\t");
        // pass label
        bw.write(var.getFilter());
        bw.write("\t");
        // INFO
        // the SV info is written here
        StringBuilder sbStr = new StringBuilder();
        // assume the SV type is consistent across alleles
        if (var.getType() == VariantOverallType.Tandem_Duplication) {
            sbStr.append("SVTYPE=DUP;");
            sbStr.append("SVLEN=");
            sbStr.append(var.getLength());
        } else if (var.getType() == VariantOverallType.Inversion) {
            sbStr.append("SVTYPE=INV;");
            sbStr.append("SVLEN=");
            sbStr.append(var.getLength());
        } else {
            sbStr.append("SVLEN=");
            sbStr.append(var.getLength());
        }
        bw.write(sbStr.toString());
        bw.write("\t");
        // label (GT)
        bw.write("GT:CN");
        bw.write("\t");
        // the genotype
        // for this one we need to work out which one is added
        sbStr = new StringBuilder();
        sbStr.append(geno0);
        sbStr.append("|");
        sbStr.append(geno1);
        sbStr.append(":");
        sbStr.append(var.getCN(geno0));
        sbStr.append("|");
        sbStr.append(var.getCN(geno1));

        bw.write(sbStr.toString());

        bw.newLine();
    }

    /**
     * The genotype is inferred from the variant in this case
     *
     * @param bw  BufferedWriter for writing
     * @param var The variant to output
     * @throws IOException
     */
    public void outputVcfRecord(BufferedWriter bw, Variant var)
            throws IOException {
        outputVcfRecord(bw, var, var.getGenotypes().geno[0], var.getGenotypes().geno[1]);
    }


}
