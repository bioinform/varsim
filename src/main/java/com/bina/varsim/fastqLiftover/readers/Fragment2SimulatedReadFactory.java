package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.ArtAlnRecord;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;

public class Fragment2SimulatedReadFactory extends SimulatedReadFactory {
    @Override
    public SimulatedRead createSimulatedRead(ArtAlnRecord alnRecord, String[] fastqEntry, int lineNumber) {
        SimulatedRead read = new SimulatedRead();
        read.fragment = 2;
        // rest of the code for handling fragment 2
        return read;
    }
}