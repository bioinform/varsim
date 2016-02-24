package com.bina.varsim.fastqLiftover.types;

import com.bina.varsim.fastqLiftover.readers.MapFileReader;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class MapBlocks {
    public final static Logger log = Logger.getLogger(MapBlocks.class.getName());
    public final static int MIN_LENGTH_INTERVAL = 10;
    public Map<String, NavigableSet<MapBlock>> chrBlocks;

    public MapBlocks(final File mapFile) throws IOException, IllegalArgumentException {
        MapFileReader mfr = new MapFileReader(mapFile);

        MapBlock mapBlock;
        chrBlocks = new HashMap<>();
        while ((mapBlock = mfr.getNext()) != null) {
            final String chromosome = mapBlock.srcLoc.chromosome;
            if (!chrBlocks.containsKey(chromosome)) {
                chrBlocks.put(chromosome, new TreeSet<MapBlock>());
            }
            chrBlocks.get(chromosome).add(mapBlock);
        }
    }

    // Half-open interval
    public List<GenomeLocation> liftOverInterval(final String chromosome, final int start, final int end, final int direction) {
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
}
