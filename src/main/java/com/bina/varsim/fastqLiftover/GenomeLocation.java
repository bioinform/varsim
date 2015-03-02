package com.bina.varsim.fastqLiftover;

import com.bina.varsim.util.*;

public class GenomeLocation implements Comparable<GenomeLocation> {
    public String chromosome;
    public int location = 0;
    public Position1 read_location;
    public int direction = 0;
    public String feature;

    public GenomeLocation(final String locationString) {
        final String fields[] = locationString.split("-", -1);
        this.chromosome = fields[0];
        String [] location_fields = fields[1].split("/");
        this.location = "".equals(location_fields[0]) ? 0 : Integer.parseInt(location_fields[0]);
        switch (location_fields.length) {
            case 1: this.read_location = null; break;
            case 2: this.read_location = new Position1(Long.parseLong(location_fields[1])); break;
            default: throw new RuntimeException("unexpected location field in " + locationString);
        }
        this.feature = "".equals(fields[2]) ? "S" : fields[2];
        this.direction = (fields.length > 3) ? 1 : 0;
    }

    public GenomeLocation(final String chromosome, final int location) {
        this.chromosome = chromosome;
        this.location = location;
        this.direction = 0;
        this.feature = "";
    }

    public GenomeLocation(final String chromosome, final int location, final int direction) {
        this.chromosome = chromosome;
        this.location = location;
        this.direction = direction;
        this.feature = "";
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
        return chromosome + "-" +
               encodeInt(location) + "/"+ (read_location!=null?Long.toString(read_location.longValue()):"") +"-" +
                ("S".equals(feature) ? "" : feature) + ((direction == 0) ? "" : "-");
    }

    @Override
    public boolean equals(Object object) {
        if (!(object instanceof GenomeLocation)) return false;
        GenomeLocation rhs = (GenomeLocation) object;
        return chromosome.equals(rhs.chromosome) && location == rhs.location && direction == rhs.direction;
    }
}
