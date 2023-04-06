package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.SimulatedRead;
import com.bina.varsim.fastqLiftover.types.SimulatedReadPair;

import java.io.BufferedReader;
import java.io.IOException;

public class DWGSIMPairedFastqReader implements PairedFastqReader {
    private DWGSIMFastqReader dwgsimFastqReader1;
    private DWGSIMFastqReader dwgsimFastqReader2;

    public DWGSIMPairedFastqReader(final BufferedReader br1, final BufferedReader br2,
            final boolean forceFiveBaseEncoding) {
        dwgsimFastqReader1 = new DWGSIMFastqReader(br1, forceFiveBaseEncoding);
        dwgsimFastqReader2 = new DWGSIMFastqReader(br2, forceFiveBaseEncoding);
    }

    public SimulatedReadPair getNextReadPair() throws IOException {
        final SimulatedRead read1 = dwgsimFastqReader1.getNextRead();
        final SimulatedRead read2 = dwgsimFastqReader2.getNextRead();

        if (read1 == null || read2 == null) {
            return null;
        }
        return new SimulatedReadPair(read1, read2);
    }
}
