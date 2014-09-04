package com.binatechnologies.varsim.fastq_liftover;

import java.io.IOException;

public interface PairedFastqReader {
    public SimulatedReadPair getNextReadPair() throws IOException;
}
