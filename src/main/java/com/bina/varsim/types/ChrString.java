package com.bina.varsim.types;

import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * Used to represent the chromosome name, able to useful things like determine sex
 * <p/>
 * Created by johnmu on 1/27/15.
 */
public class ChrString implements Comparable<ChrString>{
    private final static Logger log = Logger.getLogger(ChrString.class.getName());
    private final static Pattern restrictedChrX = Pattern.compile("(chr)?X_\\d+_\\d+");
    private final static Pattern restrictedChrY = Pattern.compile("(chr)?Y_\\d+_\\d+");
    private final static Pattern restrictedChrMT = Pattern.compile("(MT|chrM)_\\d+_\\d+");
    final String name;

    public ChrString(final String name) {
        this.name = name;
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
            return new ChrString[0];
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
        if (name.equals("X") || name.equals("chrX")) {
            return true;
        } else if (restrictedChrX.matcher(name).matches()) {
            log.warning("Treating " + name + " as restricted genome.");
            return true;
        }
        return false;
    }

    public boolean isY() {
        if (name.equals("Y") || name.equals("chrY")) {
            return true;
        } else if (restrictedChrY.matcher(name).matches()) {
            log.warning("Treating " + name + " as restricted genome.");
            return true;
        }
        return false;
    }

    public boolean isMT() {
        if (name.equals("MT") || name.equals("chrM")) {
            return true;
        } else if (restrictedChrMT.matcher(name).matches()) {
            log.warning("Treating " + name + " as restricted genome.");
            return true;
        }
        return false;
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
        if (!(o instanceof ChrString)) return false;

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
