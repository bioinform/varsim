package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.SimulatedReadPair;

import java.io.IOException;

public interface PairedFastqReader {
    public SimulatedReadPair getNextReadPair() throws IOException;
}
