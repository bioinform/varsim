package com.bina.varsim.types;

import com.bina.varsim.fastqLiftover.types.GenomeInterval;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import org.apache.log4j.Logger;

/**
 * Created by mohiyudm on 2/25/16.
 */
public class ReadMapBlock {
    private final static Logger log = Logger.getLogger(ReadMapBlock.class.getName());

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
        log.trace("Parsing " + s);
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

    public GenomeLocation getUnclippedStart() {
        return new GenomeLocation(mapInterval.getChromosome(), mapInterval.getStart() - readStart, mapInterval.getStrand(), mapInterval.getFeature());
    }
}
