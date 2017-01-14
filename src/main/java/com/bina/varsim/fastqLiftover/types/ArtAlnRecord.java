package com.bina.varsim.fastqLiftover.types;

import com.bina.varsim.types.ChrString;

public class ArtAlnRecord {
    public ChrString chromosome = new ChrString("");
    public int location = -1;
    public int direction = -1;
    public String name = null;
    public int alignedBases = 0;

    public ArtAlnRecord(final ChrString chromosome, final int location, final int direction, final String name, final int alignedBases) {
        this.chromosome = chromosome;
        this.location = location;
        this.direction = direction;
        this.name = name;
        this.alignedBases = alignedBases;
    }
}
