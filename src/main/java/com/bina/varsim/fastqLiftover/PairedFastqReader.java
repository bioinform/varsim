package com.bina.varsim.fastqLiftover;

import java.io.IOException;

public interface PairedFastqReader {
    public SimulatedReadPair getNextReadPair() throws IOException;
}
