package com.bina.varsim.fastqLiftover.types;

import com.bina.intervalTree.SimpleInterval1D;
import com.bina.intervalTree.ValueInterval1D;
import com.bina.varsim.fastqLiftover.readers.MapFileReader;
import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.ReadMapBlock;
import com.bina.varsim.types.ReadMapRecord;
import com.bina.varsim.util.chrSearchTree;
import htsjdk.tribble.annotation.Strand;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class MapBlocks {
    public final static Logger log = Logger.getLogger(MapBlocks.class.getName());
    public final static int MIN_LENGTH_INTERVAL = 10;
    public Map<ChrString, NavigableMap<MapBlock, List<MapBlock>>> chrBlocks;
    public chrSearchTree<ValueInterval1D<MapBlock>> blockIntervalTree;

    public MapBlocks(final File mapFile) throws IOException, IllegalArgumentException {
        MapFileReader mfr = new MapFileReader(mapFile);

        MapBlock mapBlock;
        chrBlocks = new HashMap<>();
        blockIntervalTree = new chrSearchTree<>(true);
        while ((mapBlock = mfr.getNext()) != null) {
            final ChrString chromosome = mapBlock.srcLoc.chromosome;
            final long start = mapBlock.srcLoc.location;
            final long end = start + mapBlock.size - 1;
            if (!chrBlocks.containsKey(chromosome)) {
                chrBlocks.put(chromosome, new TreeMap<MapBlock, List<MapBlock>>());
            }
            if (!chrBlocks.get(chromosome).containsKey(mapBlock)) {
                chrBlocks.get(chromosome).put(mapBlock, new ArrayList<>());
            }
            chrBlocks.get(chromosome).get(mapBlock).add(mapBlock);
            blockIntervalTree.put(chromosome, new ValueInterval1D<>(start, end, mapBlock));
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

        final NavigableMap<MapBlock, List<MapBlock>> blocks = chrBlocks.get(chromosome);
        final NavigableMap<MapBlock, List<MapBlock>> subset = blocks.headMap(keyEnd, true).tailMap(blocks.headMap(keyStart, true).lastKey(), true);

        Iterator<List<MapBlock>> it = subset.values().iterator();
        log.trace("Going to lift over " + chromosome + ":[" + start + "," + (end - 1) + "](" + direction + ")");

        boolean seenIns = false;
        boolean seenDel = false;
        while (it.hasNext()) {
            for (MapBlock b : it.next()) {
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
        }
        return liftedLocs;
    }

    public Collection<ReadMapBlock> liftOverGenomeInterval(final GenomeInterval interval, final int minIntervalLength) {
        final Collection<ReadMapBlock> readMapBlocks = new ArrayList<>();

        final ChrString chromosome = interval.chromosome;
        final int start = interval.start + 1; //convert to 1-based start
        final int end = interval.end; //1-based end

        if (!blockIntervalTree.containsKey(chromosome)) {
            return readMapBlocks;
        }

        final List<MapBlock> subset = findIntersectingBlocks(blockIntervalTree, chromosome, start, end);

        log.trace("Going to lift over " + interval);

        // since intervalOffset is only used for MapBlock, which is 1-based, we need to set the beginning offset to 1
        int intervalOffset = 1;
        for (MapBlock b : subset) {

                int srcStart = Math.max(start, b.srcLoc.location);
                int srcEnd = Math.min(end, b.srcLoc.location + (
                        b.blockType != MapBlock.BlockType.DEL? b.size: 0) - 1);
                final int lengthOfIntervalOnRead = srcEnd - srcStart + 1;
                int dstStart = b.dstLoc.location + srcStart - b.srcLoc.location; //1-based start on destination (ref)
                int dstEnd = dstStart - 1;
                if (b.blockType == MapBlock.BlockType.DEL) {
                    dstEnd += b.size;
                } else if (b.blockType == MapBlock.BlockType.INS) {
                    dstEnd += 0;
                } else {
                    dstEnd = Math.min(dstStart + lengthOfIntervalOnRead - 1,
                                    b.dstLoc.location + b.size - 1);
                }
                final int lengthOfIntervalOnRef = dstEnd - dstStart + 1;

                log.trace("intervalStart = " + srcStart + " intervalEnd = " + srcEnd + " lengthOfInterval = " + lengthOfIntervalOnRead);

                // lengthOfIntervalOnRef is length of the lifted over interval
                if (lengthOfIntervalOnRef < minIntervalLength) {
                    log.trace("Skipping block " + b + " since the overlap is too small ( < " + minIntervalLength + ")");
                } else {
                    // 0-based start, 1-based end
                    final GenomeInterval liftedInterval = new GenomeInterval();
                    liftedInterval.chromosome = b.dstLoc.chromosome;
                    liftedInterval.feature = b.blockType;
                    if (b.direction == 0) {
                        liftedInterval.start = b.dstLoc.location + srcStart - b.srcLoc.location - 1;
                        liftedInterval.end = liftedInterval.start + lengthOfIntervalOnRef;
                        liftedInterval.strand = interval.strand;
                    } else {
                        liftedInterval.start = b.dstLoc.location + (b.srcLoc.location + b.size - 1 - srcEnd) - 1;
                        liftedInterval.end = liftedInterval.start + lengthOfIntervalOnRef;
                        liftedInterval.strand = interval.strand == Strand.POSITIVE ? Strand.NEGATIVE : Strand.POSITIVE;
                    }

                    readMapBlocks.add(new ReadMapBlock(intervalOffset, intervalOffset + lengthOfIntervalOnRead - 1, liftedInterval));
                }
                intervalOffset += lengthOfIntervalOnRead;
        }
        return readMapBlocks;
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

    private List<MapBlock> findIntersectingBlocks(final chrSearchTree chrSearchTree, final ChrString chr, final int start, final int end) {
        List<MapBlock> overlappingMapBlocks = new ArrayList<>();
        Iterable<ValueInterval1D<MapBlock>> overlaps = chrSearchTree.getOverlaps(chr, new SimpleInterval1D(start, end));
        if (overlaps == null) {
            // nothing found
            return overlappingMapBlocks;
        }
        for (ValueInterval1D<MapBlock> i : overlaps) {
            overlappingMapBlocks.add(i.getContent());
        }
        return overlappingMapBlocks;
    }
}
