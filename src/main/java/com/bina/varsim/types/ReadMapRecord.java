package com.bina.varsim.types;

import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.google.common.base.Joiner;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by mohiyudm on 2/22/16.
 */
public class ReadMapRecord {
    protected String readName;
    private static final Joiner joiner = Joiner.on(",").skipNulls();

    protected Collection<GenomeLocation> alignmentLocations;

    public ReadMapRecord(final String readName, final Collection<GenomeLocation> alignmentLocations) {
        this.readName = readName;
        this.alignmentLocations = alignmentLocations;
    }

    public ReadMapRecord(final String line) {
        final String[] fields = line.split("\t");
        readName = fields[0];
        alignmentLocations = new ArrayList<>();
        final String locationStrings[] = fields[1].split(",", -1);
        for (String locationString : locationStrings) {
            if (!("".equals(locationString))) {
                alignmentLocations.add(new GenomeLocation(locationString));
            }
        }
    }

    public String toString() {
        return readName + "\t" + joiner.join(alignmentLocations);
    }
}
