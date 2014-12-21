package com.binatechnologies.varsim.fastqLiftover;

import java.io.IOException;

public interface PairedFastqReader {
    public SimulatedReadPair getNextReadPair() throws IOException;
}
