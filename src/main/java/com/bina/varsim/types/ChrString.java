package com.bina.varsim.types;

/**
 * Used to represent the chromosome name, able to useful things like determine sex
 * <p/>
 * Created by johnmu on 1/27/15.
 */
public class ChrString {
    final String name;

    public ChrString(String name) {
        this.name = stripChr(name);
    }

    /**
     * Converts hg19 or b37 format chromosomes to b37 format
     *
     * @param chr Chromosome name as a string
     * @return chromosome name as a string in b37 format
     */
    public static String stripChr(String chr) {
        if (chr.length() > 3 && chr.substring(0, 3).equalsIgnoreCase("chr")) {
            return chr.substring(3);
        }
        if (chr.equals("M")) {
            chr = "MT";
        }
        return chr;
    }

    public String getName() {
        return name;
    }

    @Override
    public String toString() {
        return name;
    }

    public boolean isX() {
        return name.equals("X");
    }

    public boolean isY() {
        return name.equals("Y");
    }

    public boolean isMT() {
        return name.equals("MT");
    }

    /**
     * Checks if the chromosome is haploid.
     * TODO ignore the alternate contigs for now
     *
     * @param gender Gender of individual
     * @return True if the chromosome is haploid given the sex
     */
    public boolean isHaploid(GenderType gender) {
        if (isMT()) {
            return true;
        } else if (gender == GenderType.MALE) {
            return (isX() || isY());
        }
        return false;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ChrString chrString = (ChrString) o;

        return name.equals(chrString.name);

    }

    @Override
    public int hashCode() {
        return name.hashCode();
    }
}
