package com.binatechnologies.varsim.fastq_liftover;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

public class ArtAlnReader {
    private final static Pattern splitterPattern = Pattern.compile("[\\s]+");
    private BufferedReader br;
    private String currentLine;
    private Map<String, Integer> chromosomeLengths;

    public ArtAlnReader(final BufferedReader br) throws IOException {
        this.br = br;

        chromosomeLengths = new HashMap<String, Integer>();
        readHeader();
    }

    private void readHeader() throws IOException {
        while ((currentLine = br.readLine()) != null) {
            if (!currentLine.startsWith("#") && !currentLine.startsWith("@")) {
                break;
            }
            if (currentLine.startsWith("@SQ")) {
                final String[] fields = currentLine.split("[\\s]+");
                chromosomeLengths.put(fields[1], Integer.valueOf(fields[2]));
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
        ArtAlnRecord record = new ArtAlnRecord(alnFields[0], Integer.parseInt(alnFields[2]) + 1, direction, alnFields[1], refLength);
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
