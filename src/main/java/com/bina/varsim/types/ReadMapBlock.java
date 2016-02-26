package com.bina.varsim.types;

import com.bina.varsim.fastqLiftover.types.GenomeInterval;

/**
 * Created by mohiyudm on 2/25/16.
 */
public class ReadMapBlock {
    public static final String SEPARATOR = GenomeInterval.SEPARATOR;

    protected int readStart;

    protected int readEnd;

    protected GenomeInterval mapInterval;

    public ReadMapBlock() {
    }

    public ReadMapBlock(final int readStart, final int readEnd, final GenomeInterval mapInterval) {
        this.readStart = readStart;
        this.readEnd = readEnd;
        this.mapInterval = mapInterval;
    }

    public ReadMapBlock(final String s) {
        final String[] fields = s.split(SEPARATOR, -1);
        mapInterval = new GenomeInterval(fields);
        readStart = Integer.parseInt(fields[fields.length - 2]);
        readEnd = Integer.parseInt(fields[fields.length - 1]);
    }

    public String toString() {
        return mapInterval + SEPARATOR + readStart + SEPARATOR + readEnd;
    }

    public int getReadStart() {
        return readStart;
    }

    public int getReadEnd() {
        return readEnd;
    }

    public GenomeInterval getMapInterval() {
        return mapInterval;
    }
}
