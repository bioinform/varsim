package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.ArtAlnRecord;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.regex.Pattern;

public class ARTFastqAlnReader {
    private final static Pattern nonACGTPattern = Pattern.compile("[^acgtACGT]");
    private BufferedReader fastqBr;
    private ArtAlnReader artAlnR;
    private boolean forceFiveBaseEncoding;

    public ARTFastqAlnReader(final BufferedReader alnR, final BufferedReader fastqBr, final boolean forceFiveBaseEncoding) throws IOException {
        this.fastqBr = fastqBr;
        this.artAlnR = new ArtAlnReader(alnR);
        this.forceFiveBaseEncoding = forceFiveBaseEncoding;
    }

    public SimulatedRead getNextRead() throws IOException {
        SimulatedRead read;
        final String nameLine = fastqBr.readLine();
        //TODO: add logging
        if (nameLine == null || nameLine.trim().length() < 1) return null;
        ArtAlnRecord alnRecord;
        if ((alnRecord = artAlnR.getNextAln()) == null) return null;

        final String nameFields[] = nameLine.trim().substring(1).split("[-/]");
        //TODO: add logging
        if (nameFields.length < 2) return null;
        read = new SimulatedRead();
        read.fragment = Integer.parseInt(nameFields[nameFields.length - 1]);
        read.setReadId(nameFields[nameFields.length - 2]);

        read.sequence = fastqBr.readLine().trim();
        if (forceFiveBaseEncoding) {
            read.sequence = nonACGTPattern.matcher(read.sequence).replaceAll("N");
        }

        fastqBr.readLine();
        read.quality = fastqBr.readLine().trim();
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
