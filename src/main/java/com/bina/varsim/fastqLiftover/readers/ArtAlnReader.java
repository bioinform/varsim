package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.ArtAlnRecord;
import com.bina.varsim.types.ChrString;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import java.util.regex.Pattern;

public class ArtAlnReader {
    private final static Logger log = Logger.getLogger(ArtAlnReader.class.getName());
    private final static Pattern splitterPattern = Pattern.compile("[\\s]+");
    private LineNumberReader br;
    private String currentLine;
    private Map<ChrString, Integer> chromosomeLengths;

    public ArtAlnReader(final LineNumberReader br) throws IOException {
        this.br = br;

        chromosomeLengths = new HashMap<>();
        readHeader();
    }

    private void readHeader() throws IOException {
        while ((currentLine = br.readLine()) != null) {
            if (!currentLine.startsWith("#") && !currentLine.startsWith("@")) {
                break;
            }
            if (currentLine.startsWith("@SQ")) {
                final String[] fields = currentLine.split("[\\s]+");
                chromosomeLengths.put(new ChrString(fields[1]), Integer.valueOf(fields[2]));
            }
        }
    }

    public ArtAlnRecord getNextAln() throws IOException {
        if (currentLine == null) {
            return null;
        }

        final String refAln = br.readLine().trim();
        final String readAln = br.readLine().trim();

        final String[] alnFields = splitterPattern.split(currentLine.trim().substring(1));
        int direction = (alnFields[3].equals("-")) ? 1 : 0;
        final int refLength = refAln.replace("-", "").length();
        ArtAlnRecord record = new ArtAlnRecord(new ChrString(alnFields[0]), Integer.parseInt(alnFields[2]) + 1, direction, alnFields[1], refLength);
        if (!chromosomeLengths.containsKey(record.chromosome)) {
            //minus 2 because there are two more readLine operations
            //this might be fragile
            log.warning("got nonexistent chromosome " + record.chromosome + " at aln file line " + (br.getLineNumber() - 2));
            return null;
        }
        if (record.direction == 1) {
            final int chromosomeLength = chromosomeLengths.get(record.chromosome);
            record.location = chromosomeLength - (record.location - 1) + 1 - refLength;
        }

        currentLine = br.readLine();
        return record;
    }

    public boolean hasNextAln() {
        return currentLine != null;
    }
}
