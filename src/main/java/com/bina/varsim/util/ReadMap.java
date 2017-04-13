package com.bina.varsim.util;

import com.bina.varsim.types.ReadMapRecord;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by mohiyudm on 2/29/16.
 */
public class ReadMap extends GzFileParser<ReadMapRecord>{
    private final static Logger log = Logger.getLogger(ReadMap.class.getName());

    protected Map<String, ReadMapRecord> readMap = new HashMap<>();

    public ReadMap(final File file) throws IOException {
        log.info("Loading read alignment map from " + file);
        bufferedReader = new BufferedReader(new InputStreamReader(decompressStream(file)));
        while (hasMoreInput()) {
            readLine();
            log.trace("Parsing line " + line);
            if (line == null)
                break;
            if (line.trim().isEmpty())
                continue;
            final ReadMapRecord record = parseLine();
            readMap.put(record.getReadName(), record);
        }
        log.info("Done loading read alignment map");
    }

    public ReadMapRecord parseLine() {
        return hasMoreInput() ? new ReadMapRecord(line.trim()) : null;
    }

    public ReadMapRecord getReadMapRecord(final String readName) {
        return readMap.containsKey(readName) ? readMap.get(readName) : null;
    }
}
