package com.bina.varsim.types;

/**
 * Used to represent the chromosome name, able to useful things like determine sex
 * <p/>
 * Created by johnmu on 1/27/15.
 */
public class ChrString implements Comparable<ChrString>{
    final String name;

    public ChrString(final String name) {
        this.name = stripChr(name);
    }

    /**
     * Converts hg19 or b37 format chromosomes to b37 format
     *
     * @param chr Chromosome name as a string
     * @return chromosome name as a string in b37 format
     */
    public static String stripChr(final String chr) {
        if (chr.length() > 3 && chr.substring(0, 3).equalsIgnoreCase("chr")) {
            return chr.substring(3);
        }
        if (chr.equals("M")) {
            return "MT";
        }
        return chr;
    }

    /**
     * convert string array to ChrString array
     * @param s
     * @return
     */
    public static ChrString[] string2ChrString(final String[] s) {
        if (s == null) {
            return null;
        }
        ChrString[] chrStrings = new ChrString[s.length];
        for (int i = 0; i < s.length; i++) {
            chrStrings[i] = new ChrString(s[i]);
        }
        return chrStrings;
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
    public boolean isHaploid(final GenderType gender) {
        if (isMT()) {
            return true;
        } else if (gender == GenderType.MALE) {
            return (isX() || isY());
        }
        return false;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ChrString chrString = (ChrString) o;

        return name.equals(chrString.name);

    }

    @Override
    public int hashCode() {
        return name.hashCode();
    }

    @Override
    public int compareTo(final ChrString other) {
        return name.compareTo(other.getName());
    }
}
