package com.bina.varsim.fastqLiftover.types;

import com.bina.varsim.types.ChrString;
import com.google.common.base.Joiner;
import htsjdk.tribble.annotation.Strand;

public class GenomeInterval implements Comparable<GenomeInterval> {
    public static final String SEPARATOR = ",";
    private static final Joiner joiner = Joiner.on(SEPARATOR).skipNulls();

    protected ChrString chromosome;
    protected int start = -1; //0-based start
    protected int end = -1; //1-based end
    protected Strand strand = Strand.NONE;
    protected MapBlock.BlockType feature = MapBlock.BlockType.UNKNOWN;

    public GenomeInterval() {
    }

    public GenomeInterval(final String fields[]) {
        chromosome = new ChrString(fields[0]);
        start = Integer.parseInt(fields[1]);
        end = Integer.parseInt(fields[2]);
        strand = Strand.valueOf(fields[3]);
        feature = MapBlock.BlockType.valueOf(fields[4]);
    }

    public GenomeInterval(final ChrString chromosome, final int start, final int end, final Strand strand, final MapBlock.BlockType feature) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.feature = feature;
    }

    @Override
    public int compareTo(final GenomeInterval rhs) {
        final int chrCmp = this.chromosome.compareTo(rhs.chromosome);
        if (chrCmp != 0) {
            return chrCmp;
        }
        if (start != rhs.start) {
            return Integer.compare(start, rhs.start);
        }
        return Integer.compare(end, rhs.end);
    }

    @Override
    public String toString() {
        return joiner.join(chromosome, start, end, strand.name(), feature.name());
    }

    @Override
    public boolean equals(Object object) {
        if (this == object) return true; //this is purely for performance sake
        if (!(object instanceof GenomeInterval)) return false;
        GenomeInterval rhs = (GenomeInterval) object;
        return chromosome.equals(rhs.chromosome) && start == rhs.start && end == rhs.end && strand == rhs.strand && feature == rhs.feature;
    }

    public ChrString getChromosome() {
        return chromosome;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public Strand getStrand() {
        return strand;
    }

    public MapBlock.BlockType getFeature() {
        return feature;
    }
}
