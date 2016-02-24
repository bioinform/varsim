package com.bina.varsim.fastqLiftover.types;

public class GenomeLocation implements Comparable<GenomeLocation> {
    public String chromosome;
    public int location = 0;
    public int direction = 0;
    public MapBlock.BlockType feature;

    public GenomeLocation(final String locationString) {
        final String fields[] = locationString.split("-", -1);
        this.chromosome = fields[0];
        this.location = "".equals(fields[1]) ? 0 : Integer.parseInt(fields[1]);
        this.feature = MapBlock.BlockType.fromName(fields[2]);
        this.direction = (fields.length > 3) ? 1 : 0;
    }

    public GenomeLocation(final String chromosome, final int location) {
        this.chromosome = chromosome;
        this.location = location;
        this.direction = 0;
        this.feature = MapBlock.BlockType.SEQ;
    }

    public GenomeLocation(final String chromosome, final int location, final int direction) {
        this.chromosome = chromosome;
        this.location = location;
        this.direction = direction;
        this.feature = MapBlock.BlockType.SEQ;
    }

    public boolean isClose(final String chromosome, final int location, final int wiggle) {
        return this.chromosome.equals(chromosome) && Math.abs(this.location - location) <= wiggle;
    }

    public boolean isClose(final GenomeLocation other, final int wiggle) {
        return chromosome.equals(other.chromosome) && Math.abs(location - other.location) <= wiggle;
    }

    @Override
    public int compareTo(final GenomeLocation rhs) {
        final int chrCmp = this.chromosome.compareTo(rhs.chromosome);
        if (chrCmp != 0) {
            return chrCmp;
        }
        return Integer.compare(location, rhs.location);
    }

    private String encodeInt(int n) {
        return (n == 0) ? "" : Integer.toString(n);
    }

    public String toString() {
        return chromosome + "-" + encodeInt(location) + "-" + (feature == MapBlock.BlockType.SEQ ? "" : feature.getShortName()) + ((direction == 0) ? "" : "-");
    }

    @Override
    public boolean equals(Object object) {
        if (!(object instanceof GenomeLocation)) return false;
        GenomeLocation rhs = (GenomeLocation) object;
        return chromosome.equals(rhs.chromosome) && location == rhs.location && direction == rhs.direction;
    }
}
