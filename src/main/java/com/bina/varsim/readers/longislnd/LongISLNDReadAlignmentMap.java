package com.bina.varsim.readers.longislnd;

import com.bina.varsim.types.ReadMapRecord;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * Created by mohiyudm on 2/22/16.
 */
public class LongISLNDReadAlignmentMap {
    private final static Logger log = Logger.getLogger(LongISLNDReadAlignmentMap.class.getName());

    final protected Map<String, ReadMapRecord> readAlignmentMap = new TreeMap<>();

    public LongISLNDReadAlignmentMap(final Collection<File> readAlignmentMapFiles) throws IOException {
        for (final File readAlignmentMapFile : readAlignmentMapFiles) {
            log.info("Reading in read map from " + readAlignmentMapFile.getName());
            final AbstractFeatureReader<BEDFeature, LineIterator> featureReader = AbstractFeatureReader.getFeatureReader(readAlignmentMapFile.getName(), new BEDCodec(), false);
            try {
                final CloseableTribbleIterator<BEDFeature> featureIterator = featureReader.iterator();
                while (featureIterator.hasNext()) {
                    final BEDFeature feature = featureIterator.next();
                    readAlignmentMap.put(feature.getName(), new LongISLNDReadMapRecord(feature).toReadMapRecord());
                }
            } finally {
                featureReader.close();
            }
        }
    }

    public void writeReadMap(final File outFile) throws IOException {
        final PrintStream ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFile)));
        for (final Map.Entry<String, ReadMapRecord> entry : readAlignmentMap.entrySet()) {
            ps.println(entry.getValue());
        }
        ps.close();
    }

    public Map<String, ReadMapRecord> getReadAlignmentMap() {
        return readAlignmentMap;
    }
}
