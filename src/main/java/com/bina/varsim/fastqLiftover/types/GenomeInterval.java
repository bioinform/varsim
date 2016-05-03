package com.bina.varsim.fastqLiftover.types;

import com.google.common.base.Joiner;
import htsjdk.tribble.annotation.Strand;

public class GenomeInterval implements Comparable<GenomeInterval> {
    public static final String SEPARATOR = ",";
    private static final Joiner joiner = Joiner.on(SEPARATOR).skipNulls();

    protected String chromosome;
    protected int start = -1;
    protected int end = -1;
    protected Strand strand = Strand.NONE;
    protected MapBlock.BlockType feature = MapBlock.BlockType.UNKNOWN;

    public GenomeInterval() {
    }

    public GenomeInterval(final String intervalString) {
        final String fields[] = intervalString.split(SEPARATOR, -1);
        chromosome = fields[0];
        start = Integer.parseInt(fields[1]);
        end = Integer.parseInt(fields[2]);
        strand = Strand.valueOf(fields[3]);
        feature = MapBlock.BlockType.valueOf(fields[4]);
    }

    public GenomeInterval(final String fields[]) {
        chromosome = fields[0];
        start = Integer.parseInt(fields[1]);
        end = Integer.parseInt(fields[2]);
        strand = Strand.valueOf(fields[3]);
        feature = MapBlock.BlockType.valueOf(fields[4]);
    }

    public GenomeInterval(final String chromosome, final int start, final int end, final Strand strand, final MapBlock.BlockType feature) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.feature = feature;
    }

    public GenomeInterval(final String chromosome, final int start, final int end, final MapBlock.BlockType feature) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
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

    public String toString() {
        return joiner.join(chromosome, start, end, strand.name(), feature.name());
    }

    @Override
    public boolean equals(Object object) {
        if (!(object instanceof GenomeInterval)) return false;
        GenomeInterval rhs = (GenomeInterval) object;
        return chromosome.equals(rhs.chromosome) && start == rhs.start && end == rhs.end && strand == rhs.strand && feature == rhs.feature;
    }

    public String getChromosome() {
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
