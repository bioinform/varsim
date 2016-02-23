package com.bina.varsim.readers.longislnd;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by mohiyudm on 2/22/16.
 */
public class LongISLNDReadAlignmentMap {
    final protected Map<String, LongISLNDReadMapRecord> readAlignmentMap = new TreeMap<>();

    public LongISLNDReadAlignmentMap(final File readAlignmentMapFile) throws IOException {
        final AbstractFeatureReader<BEDFeature, LineIterator> featureReader = AbstractFeatureReader.getFeatureReader(readAlignmentMapFile.getName(), new BEDCodec(), false);
        try {
            final CloseableTribbleIterator<BEDFeature> featureIterator = featureReader.iterator();
            while (featureIterator.hasNext()) {
                final BEDFeature feature = featureIterator.next();
                readAlignmentMap.put(feature.getName(), new LongISLNDReadMapRecord(feature));
            }
        } finally {
            featureReader.close();
        }
    }

    public void writeReadMap(final File outFile) throws IOException {
        final PrintStream ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFile)));
        for (final Map.Entry<String, LongISLNDReadMapRecord> entry : readAlignmentMap.entrySet()) {
            ps.println(entry.getValue());
        }
        ps.close();
    }
}
