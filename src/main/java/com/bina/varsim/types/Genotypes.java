package com.bina.varsim.types;

import java.util.Random;

/**
 * Stores genotypes and provides facilities for generating random ones
 */
public class Genotypes {
    public byte geno[];

    // generate random genotypes
    // make sure het-hom ration is 1.5


    /**
     * generate random genotypes, make sure het-hom ration is 1.5
     *
     * @param chr Chromosome name in number format
     * @param num_alt Number of alternate alleles possible
     * @param rand    Random number generator
     */
    Genotypes(ChrString chr, GenderType male, int num_alt, Random rand) {
        this(chr, male, num_alt, rand, 0.6);
    }

    /**
     * generate random genotypes
     *
     * @param chr  Chromosome name in number format
     * @param num_alt  Number of alternate alleles possible
     * @param rand     Random number generator
     * @param prop_het proportion heterozygous, rest are homo
     */
    Genotypes(ChrString chr, GenderType male, int num_alt, Random rand, double prop_het) {
        geno = new byte[2];

        if (chr.isHaploid(male)) {
            geno[0] = 1;
            geno[1] = 1;
        } else {
            if (rand.nextDouble() > prop_het) {
                // hom
                geno[0] = (byte) (rand.nextInt(num_alt) + 1);
                geno[1] = geno[0];
            } else {
                // het
                if (rand.nextInt(2) == 0) { // random phasing
                    geno[0] = (byte) rand.nextInt(num_alt + 1);
                    geno[1] = (byte) rand.nextInt(num_alt + 1);
                    while (geno[1] == geno[0]) {
                        geno[1] = (byte) rand.nextInt(num_alt + 1);
                    }
                } else {
                    geno[1] = (byte) rand.nextInt(num_alt + 1);
                    geno[0] = (byte) rand.nextInt(num_alt + 1);
                    while (geno[1] == geno[0]) {
                        geno[0] = (byte) rand.nextInt(num_alt + 1);
                    }
                }
            }


        }
    }

    /**
     * @param geno0 First allele
     * @param geno1 Second allele
     */
    Genotypes(byte geno0, byte geno1) {
        geno = new byte[2];
        this.geno[0] = geno0;
        this.geno[1] = geno1;
    }

    /**
     * @return false if both genotypes are zero
     */
    public boolean isNonRef() {
        return !(geno[0] == 0 && geno[1] == 0);
    }

    public String toString() {
        return "(" + String.valueOf(geno[0]) + ":" + String.valueOf(geno[1]) + ")";
    }
}
