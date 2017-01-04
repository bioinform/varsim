package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;
import com.google.common.base.Joiner;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

public class DWGSIMFastqReader {
    private final static Pattern nonACGTPattern = Pattern.compile("[^acgtACGT]");
    private final static Joiner joiner = Joiner.on("_").skipNulls();
    private final static Logger log = Logger.getLogger(DWGSIMFastqReader.class.getName());
    private BufferedReader br;
    private boolean forceFiveBaseEncoding;

    public DWGSIMFastqReader(final BufferedReader br, final boolean forceFiveBaseEncoding) {
        this.br = br;
        this.forceFiveBaseEncoding = forceFiveBaseEncoding;
    }

    public SimulatedRead getNextRead() throws IOException {
        SimulatedRead read;
        final String nameLine = br.readLine();
        if (nameLine == null) return null;

        final String nameFields[] = StringUtils.replaceChars(nameLine.trim().substring(1), "/_", "::").split(":");
        read = new SimulatedRead();
        read.fragment = Integer.parseInt(nameFields[nameFields.length - 1]);
        read.setReadId(nameFields[nameFields.length - 2]);
        read.indels2 = Integer.parseInt(nameFields[nameFields.length - 3]);
        read.snps2 = Integer.parseInt(nameFields[nameFields.length - 4]);
        read.seqErrors2 = Integer.parseInt(nameFields[nameFields.length - 5]);
        read.indels1 = Integer.parseInt(nameFields[nameFields.length - 6]);
        read.snps1 = Integer.parseInt(nameFields[nameFields.length - 7]);
        read.seqErrors1 = Integer.parseInt(nameFields[nameFields.length - 8]);
        read.random2 = Integer.parseInt(nameFields[nameFields.length - 9]);
        read.random1 = Integer.parseInt(nameFields[nameFields.length - 10]);

        final int direction2 = Integer.parseInt(nameFields[nameFields.length - 11]);
        final int direction1 = Integer.parseInt(nameFields[nameFields.length - 12]);
        final int start2 = Integer.parseInt(nameFields[nameFields.length - 13]);
        final int start1 = Integer.parseInt(nameFields[nameFields.length - 14]);
        final String chromosome1 = joiner.join(Arrays.copyOfRange(nameFields, 0, nameFields.length - 14));
        final String chromosome2 = chromosome1; // Marghoob was there a reason for this?

        read.locs1.add(new GenomeLocation(chromosome1, start1, direction1));
        read.locs2.add(new GenomeLocation(chromosome2, start2, direction2));
        read.origLocs1.add(new GenomeLocation(chromosome1, start1, direction1));
        read.origLocs2.add(new GenomeLocation(chromosome2, start2, direction2));

        read.sequence = br.readLine().trim();
        if (forceFiveBaseEncoding) {
            read.sequence = nonACGTPattern.matcher(read.sequence).replaceAll("N");
        }

        br.readLine();
        read.quality = br.readLine().trim();
        if (read.fragment == 1) {
            read.alignedBases1 = read.sequence.length();
        } else {
            read.alignedBases2 = read.sequence.length();
        }

        return read;
    }
}
