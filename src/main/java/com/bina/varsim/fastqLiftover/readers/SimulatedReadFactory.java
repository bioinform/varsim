package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.ArtAlnRecord;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;

public abstract class SimulatedReadFactory {
    public abstract SimulatedRead createSimulatedRead(ArtAlnRecord alnRecord, String[] fastqEntry, int lineNumber);
}
