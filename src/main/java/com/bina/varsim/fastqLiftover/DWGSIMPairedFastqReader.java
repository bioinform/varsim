package com.bina.varsim.fastqLiftover;

import java.io.BufferedReader;
import java.io.IOException;

public class DWGSIMPairedFastqReader implements PairedFastqReader {
    private DWGSIMFastqReader r1;
    private DWGSIMFastqReader r2;

    public DWGSIMPairedFastqReader(final BufferedReader br1, final BufferedReader br2, final boolean forceFiveBaseEncoding) {
        r1 = new DWGSIMFastqReader(br1, forceFiveBaseEncoding);
        r2 = new DWGSIMFastqReader(br2, forceFiveBaseEncoding);
    }

    public SimulatedReadPair getNextReadPair() throws IOException {
        final SimulatedRead read1 = r1.getNextRead();
        final SimulatedRead read2 = r2.getNextRead();

        if (read1 == null || read2 == null) {
            return null;
        }
        return new SimulatedReadPair(read1, read2);
    }
}
