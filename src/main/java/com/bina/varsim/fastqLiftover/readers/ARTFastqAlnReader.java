package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.ArtAlnRecord;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;

import java.io.IOException;
import java.io.LineNumberReader;
import org.apache.log4j.Logger;
import java.util.regex.Pattern;

public class ARTFastqAlnReader {
    private final static Logger log = Logger.getLogger(ARTFastqAlnReader.class.getName());
    private final static Pattern nonACGTPattern = Pattern.compile("[^acgtACGT]");
    private FastqReader fastqReader;
    private ArtAlnReader artAlnR;
    private boolean forceFiveBaseEncoding;

    public ARTFastqAlnReader(final LineNumberReader alnR, final LineNumberReader fastqBr, final boolean forceFiveBaseEncoding) throws IOException {
        this.fastqReader = new FastqReader(fastqBr);
        this.artAlnR = new ArtAlnReader(alnR);
        this.forceFiveBaseEncoding = forceFiveBaseEncoding;
    }

    public SimulatedRead getNextRead() throws IOException {
        SimulatedRead read;
        String[] fastqEntry = fastqReader.getNextFastqEntry();

        if (fastqEntry == null) {
            return null;
        }

        final String nameLine = fastqEntry[0];
        if (nameLine.trim().length() < 1) {
            log.error("got empty name string at FASTQ line " + fastqReader.getLineNumber());
            return null;
        }
        ArtAlnRecord alnRecord;
        if ((alnRecord = artAlnR.getNextAln()) == null) {
            return null;
        }

        final String nameFields[] = nameLine.trim().substring(1).split("[-/]");
        if (nameFields.length < 2) {
            log.warn("expect at least 2 fields at fastq line " + fastqReader.getLineNumber() + ": " + nameLine);
            return null;
        }

        SimulatedReadFactory factory;
        int fragment = Integer.parseInt(nameFields[nameFields.length - 1]);
        if (fragment == 1) {
            factory = new Fragment1SimulatedReadFactory();
        } else {
            factory = new Fragment2SimulatedReadFactory();
        }

        return factory.createSimulatedRead(alnRecord, fastqEntry, fastqReader.getLineNumber());
    }

}