package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.ArtAlnRecord;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.logging.Logger;
import java.util.regex.Pattern;

public class ARTFastqAlnReader {
    private final static Logger log = Logger.getLogger(ARTFastqAlnReader.class.getName());
    private final static Pattern nonACGTPattern = Pattern.compile("[^acgtACGT]");
    private BufferedReader fastqBr;
    private ArtAlnReader artAlnR;
    private boolean forceFiveBaseEncoding;
    private long fqLineCount;

    public ARTFastqAlnReader(final BufferedReader alnR, final BufferedReader fastqBr, final boolean forceFiveBaseEncoding) throws IOException {
        this.fastqBr = fastqBr;
        this.artAlnR = new ArtAlnReader(alnR);
        this.forceFiveBaseEncoding = forceFiveBaseEncoding;
        //assume no one has read from the above BufferedReaders
        //this might not be true...
        fqLineCount = 0;
    }

    public SimulatedRead getNextRead() throws IOException {
        SimulatedRead read;
        final String nameLine = fastqBr.readLine();
        fqLineCount++;
        //TODO: add logging
        if (nameLine == null || nameLine.trim().length() < 1) {
            if (nameLine.trim().length() < 1)
                log.warning("got empty string at FASTQ line " + fqLineCount);
            return null;
        }
        ArtAlnRecord alnRecord;
        if ((alnRecord = artAlnR.getNextAln()) == null) {
            return null;
        }

        final String nameFields[] = nameLine.trim().substring(1).split("[-/]");
        if (nameFields.length < 2) {
            log.warning("expect at least 2 fields at fastq line " + fqLineCount + ": " + nameLine);
            return null;
        }
        read = new SimulatedRead();
        read.fragment = Integer.parseInt(nameFields[nameFields.length - 1]);
        read.setReadId(nameFields[nameFields.length - 2]);

        read.sequence = fastqBr.readLine().trim();
        fqLineCount++;
        if (forceFiveBaseEncoding) {
            read.sequence = nonACGTPattern.matcher(read.sequence).replaceAll("N");
        }

        fastqBr.readLine();
        fqLineCount++;
        read.quality = fastqBr.readLine().trim();
        fqLineCount++;
        if (read.fragment == 1) {
            read.locs1.add(new GenomeLocation(alnRecord.chromosome, alnRecord.location, alnRecord.direction));
            read.origLocs1.add(new GenomeLocation(alnRecord.chromosome, alnRecord.location, alnRecord.direction));
            read.alignedBases1 = alnRecord.alignedBases;
        } else {
            read.locs2.add(new GenomeLocation(alnRecord.chromosome, alnRecord.location, alnRecord.direction));
            read.origLocs2.add(new GenomeLocation(alnRecord.chromosome, alnRecord.location, alnRecord.direction));
            read.alignedBases2 = alnRecord.alignedBases;
        }

        return read;
    }
}
