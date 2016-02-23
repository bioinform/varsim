package com.bina.varsim.readers.longislnd;

import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MapBlocks;
import com.bina.varsim.types.ReadMapRecord;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;

import java.util.ArrayList;
import java.util.Collection;

public class LongISLNDReadMapRecord {
    final protected String chromosome;

    final protected int start;

    final protected int end;

    final protected String readName;

    final protected Strand strand;

    LongISLNDReadMapRecord(final BEDFeature bedFeature) {
        chromosome = bedFeature.getContig();
        start = bedFeature.getStart();
        end = bedFeature.getEnd();
        readName = bedFeature.getName();
        strand = bedFeature.getStrand();
    }

    public String getChromosome() {
        return chromosome;
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return end;
    }

    public String getReadName() {
        return readName;
    }

    public Strand getStrand() {
        return strand;
    }

    public ReadMapRecord toReadMapRecord(final MapBlocks mapBlocks) {
        return new ReadMapRecord(readName, mapBlocks.liftOverInterval(chromosome, start, end, 0));
    }
}
