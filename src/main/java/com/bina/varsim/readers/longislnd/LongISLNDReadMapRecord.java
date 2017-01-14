package com.bina.varsim.readers.longislnd;

import com.bina.varsim.fastqLiftover.types.GenomeInterval;
import com.bina.varsim.fastqLiftover.types.MapBlock;
import com.bina.varsim.fastqLiftover.types.MapBlocks;
import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.ReadMapBlock;
import com.bina.varsim.types.ReadMapRecord;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;

import java.util.Collection;
import java.util.Collections;


public class LongISLNDReadMapRecord {
    final protected ChrString chromosome;

    final protected int start;

    final protected int end;

    final protected String readName;

    final protected Strand strand;

    LongISLNDReadMapRecord(final BEDFeature bedFeature) {
        chromosome = new ChrString(bedFeature.getContig());
        start = bedFeature.getStart();
        end = bedFeature.getEnd();
        readName = bedFeature.getName();
        strand = bedFeature.getStrand();
    }

    public ChrString getChromosome() {
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

    public ReadMapRecord toReadMapRecord() {
        return new ReadMapRecord(readName, Collections.singletonList(new ReadMapBlock(0, end - start, new GenomeInterval(chromosome, start, end, strand, MapBlock.BlockType.UNKNOWN))));
    }
}
