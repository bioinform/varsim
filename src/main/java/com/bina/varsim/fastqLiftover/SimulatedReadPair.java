package com.bina.varsim.fastqLiftover;

public class SimulatedReadPair {
    public SimulatedRead read1;
    public SimulatedRead read2;

    public SimulatedReadPair(final SimulatedRead read1, final SimulatedRead read2) {
        this.read1 = read1;
        this.read2 = read2;

        this.read1.locs2 = this.read2.locs2;
        this.read1.origLocs2 = this.read2.origLocs2;
        this.read2.locs1 = this.read1.locs1;
        this.read2.origLocs1 = this.read1.origLocs1;
        this.read1.alignedBases2 = this.read2.alignedBases2;
        this.read2.alignedBases1 = this.read1.alignedBases1;
    }
}
