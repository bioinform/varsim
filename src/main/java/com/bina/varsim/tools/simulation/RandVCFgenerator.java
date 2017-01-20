package com.bina.varsim.tools.simulation;

import com.bina.varsim.VarSimTool;
import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.variant.VariantOverallType;
import com.bina.varsim.types.SampleParams;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.types.variant.alt.Alt;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.Random;


abstract public class RandVCFgenerator extends VarSimTool {

    protected Random rand;

    static private final long DEFAULT_SEED = 3333;

    @Option(name = "-seed", usage = "Seed for random number generator")
    static protected long seed = DEFAULT_SEED;

    /**
     * This sets a default seed, ideally different between different runs
     */
    public RandVCFgenerator(final String command, final String description) {
        super(command, description);
        rand = new Random(DEFAULT_SEED);
    }

    /**
     *
     * @param  file Sequence file
     * @return Sequence file as a byte array
     */
    byte[] fileToByteArray(final File file) {
        byte[] array = new byte[0];

        try {
            FileReader fileReader = new FileReader(file);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            StringBuilder sb = new StringBuilder();
            String line;
            while ((line = bufferedReader.readLine()) != null) {
                line = line.trim();
                sb.append(line);
            }
            bufferedReader.close();
            array = sb.toString().getBytes("US-ASCII");
        } catch (IOException e) {
            e.printStackTrace();
        }

        return array;
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
        Alt alt = var.getAlt(geno);
        if (alt != null) {
            if (alt.getSeqType() == FlexSeq.Type.INS) {
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
                alt = new Alt(new FlexSeq(newSeq));
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

        //skip checking because it does not apparently harm to output nonATCGN ref and alt alleles
        /*
        // ignore ACTGN
        String ref = var.getReferenceString().toUpperCase();
        String alt = var.alternativeAlleleString().toUpperCase();

        if (!ref.matches("[ACTGN]*") || (!var.isAltACTGN())) {
            return; // don't output if it is not ACTGN
        }
        */
        bw.write(var.toString(geno0, geno1));
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
