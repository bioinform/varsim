package com.bina.varsim.fastqLiftover.types;

import com.bina.varsim.fastqLiftover.readers.MapFileReader;
import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.ReadMapBlock;
import com.bina.varsim.types.ReadMapRecord;
import htsjdk.tribble.annotation.Strand;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class MapBlocks {
    public final static Logger log = Logger.getLogger(MapBlocks.class.getName());
    public final static int MIN_LENGTH_INTERVAL = 10;
    public Map<ChrString, NavigableSet<MapBlock>> chrBlocks;

    public MapBlocks(final File mapFile) throws IOException, IllegalArgumentException {
        MapFileReader mfr = new MapFileReader(mapFile);

        MapBlock mapBlock;
        chrBlocks = new HashMap<>();
        while ((mapBlock = mfr.getNext()) != null) {
            final ChrString chromosome = mapBlock.srcLoc.chromosome;
            if (!chrBlocks.containsKey(chromosome)) {
                chrBlocks.put(chromosome, new TreeSet<MapBlock>());
            }
            chrBlocks.get(chromosome).add(mapBlock);
        }
    }

    // Half-open interval
    public List<GenomeLocation> liftOverInterval(final ChrString chromosome, final int start, final int end, final int direction) {
        final MapBlock keyStart = new MapBlock(new GenomeLocation(chromosome, start));
        final MapBlock keyEnd = new MapBlock(new GenomeLocation(chromosome, end - 1));

        List<GenomeLocation> liftedLocs = new ArrayList<>();
        if (!chrBlocks.containsKey(chromosome)) {
            return liftedLocs;
        }

        final NavigableSet<MapBlock> blocks = chrBlocks.get(chromosome);
        final NavigableSet<MapBlock> subset = blocks.headSet(keyEnd, true).tailSet(blocks.headSet(keyStart, true).last(), true);

        Iterator<MapBlock> it = subset.iterator();
        log.trace("Going to lift over " + chromosome + ":[" + start + "," + (end - 1) + "](" + direction + ")");

        boolean seenIns = false;
        boolean seenDel = false;
        while (it.hasNext()) {
            final MapBlock b = it.next();
            if (b.blockType == MapBlock.BlockType.INS || b.blockType == MapBlock.BlockType.DEL) {
                if ((b.blockType == MapBlock.BlockType.INS && !seenIns) || (b.blockType == MapBlock.BlockType.DEL && !seenDel)) {
                    GenomeLocation liftedLoc = new GenomeLocation(b.dstLoc.chromosome, b.dstLoc.location);
                    liftedLoc.feature = b.blockType;
                    liftedLoc.direction = direction;
                    liftedLocs.add(liftedLoc);
                    if (b.blockType == MapBlock.BlockType.INS) seenIns = true;
                    if (b.blockType == MapBlock.BlockType.DEL) seenDel = true;
                } else {
                    log.trace("Skipping block " + b + " since it's an INS/DEL");
                }
                continue;
            }

            int intervalStart = Math.max(start, b.srcLoc.location);
            int intervalEnd = Math.min(end - 1, b.srcLoc.location + b.size - 1);
            int lengthOfInterval = intervalEnd - intervalStart + 1;

            log.trace("intervalStart = " + intervalStart + " intervalEnd = " + intervalEnd + " lengthOfInterval = " + lengthOfInterval);

            if (lengthOfInterval < MIN_LENGTH_INTERVAL) {
                log.trace("Skipping block " + b + " since the overlap is too small ( < " + MIN_LENGTH_INTERVAL + ")");
                continue;
            }

            GenomeLocation liftedLoc = new GenomeLocation(b.dstLoc.chromosome, b.dstLoc.location);
            liftedLoc.feature = b.blockType;
            if (b.direction == 0) {
                liftedLoc.location = b.dstLoc.location + start - b.srcLoc.location;
                liftedLoc.direction = direction;
            } else {
                liftedLoc.location = b.dstLoc.location + b.size - (end - (b.srcLoc.location + 1));
                liftedLoc.direction = 1 - direction;
            }
            log.trace(chromosome + ":[" + intervalStart + "," + intervalEnd + "](" + direction + ") lifted to " + liftedLoc + " using block " + b);
            liftedLocs.add(liftedLoc);
        }
        return liftedLocs;
    }

    public Collection<ReadMapBlock> liftOverGenomeInterval(final GenomeInterval interval, final int minIntervalLength) {
        final Collection<ReadMapBlock> readMapBlocks = new ArrayList<>();

        final ChrString chromosome = interval.chromosome;
        final int start = interval.start;
        final int end = interval.end;

        final MapBlock keyStart = new MapBlock(new GenomeLocation(chromosome, start));
        final MapBlock keyEnd = new MapBlock(new GenomeLocation(chromosome, end - 1));

        if (!chrBlocks.containsKey(chromosome)) {
            return readMapBlocks;
        }

        final NavigableSet<MapBlock> blocks = chrBlocks.get(chromosome);
        final NavigableSet<MapBlock> subset = blocks.headSet(keyEnd, true).tailSet(blocks.headSet(keyStart, true).last(), true);

        Iterator<MapBlock> it = subset.iterator();
        log.trace("Going to lift over " + interval);

        int intervalOffset = 0;
        while (it.hasNext()) {
            final MapBlock b = it.next();

            int srcStart = Math.max(start, b.srcLoc.location);
            int srcEnd = Math.min(end - 1, b.srcLoc.location + b.size - 1);
            int lengthOfInterval = srcEnd - srcStart + 1;
            final int lengthOfIntervalOnRead = b.blockType != MapBlock.BlockType.DEL ? lengthOfInterval : 0;
            final int lengthOfIntervalOnRef = b.blockType != MapBlock.BlockType.INS ? lengthOfInterval : 0;

            log.trace("intervalStart = " + srcStart + " intervalEnd = " + srcEnd + " lengthOfInterval = " + lengthOfInterval);

            if (lengthOfInterval < minIntervalLength) {
                log.trace("Skipping block " + b + " since the overlap is too small ( < " + MIN_LENGTH_INTERVAL + ")");
            } else {
                final GenomeInterval liftedInterval = new GenomeInterval();
                liftedInterval.chromosome = b.dstLoc.chromosome;
                liftedInterval.feature = b.blockType;
                if (b.direction == 0) {
                    liftedInterval.start = b.dstLoc.location + srcStart - b.srcLoc.location;
                    liftedInterval.end = liftedInterval.start + lengthOfIntervalOnRef;
                    liftedInterval.strand = interval.strand;
                } else {
                    liftedInterval.start = b.dstLoc.location + (b.srcLoc.location + b.size - 1 - srcEnd);
                    liftedInterval.end = liftedInterval.start + lengthOfIntervalOnRef;
                    liftedInterval.strand = interval.strand == Strand.POSITIVE ? Strand.NEGATIVE : Strand.POSITIVE;
                }

                readMapBlocks.add(new ReadMapBlock(intervalOffset, intervalOffset + lengthOfIntervalOnRead, liftedInterval));
            }
            intervalOffset += lengthOfIntervalOnRead;
        }

        return readMapBlocks;
    }

    public Collection<ReadMapBlock> liftOverGenomeInterval(final GenomeInterval interval) {
        return liftOverGenomeInterval(interval, MIN_LENGTH_INTERVAL);
    }

    public ReadMapRecord liftOverReadMapRecord(final ReadMapRecord readMapRecord, final int minIntervalLength) {
        final List<Collection<ReadMapBlock>> liftedReadMaps = new ArrayList<>();
        for (final Collection<ReadMapBlock> readMapBlocks : readMapRecord.getMultiReadMapBlocks()) {
            final Collection<ReadMapBlock> liftedReadMapBlocks = new ArrayList<>();
            for (final ReadMapBlock readMapBlock : readMapBlocks) {
                final int offset = readMapBlock.getReadStart();
                for (final ReadMapBlock liftedReadMapBlock : liftOverGenomeInterval(readMapBlock.getMapInterval(), minIntervalLength)) {
                    liftedReadMapBlocks.add(new ReadMapBlock(offset + liftedReadMapBlock.getReadStart(), offset + liftedReadMapBlock.getReadEnd(), liftedReadMapBlock.getMapInterval()));
                }
            }
            liftedReadMaps.add(liftedReadMapBlocks);
        }
        return new ReadMapRecord(readMapRecord.getReadName(), liftedReadMaps);
    }

    public ReadMapRecord liftOverReadMapRecord(final ReadMapRecord readMapRecord) {
        return liftOverReadMapRecord(readMapRecord, MIN_LENGTH_INTERVAL);
    }
}
