package com.bina.varsim.fastqLiftover;

import com.bina.varsim.fastqLiftover.types.GenomeInterval;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MapBlock;
import htsjdk.tribble.annotation.Strand;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for GenomeLocation
 */
public class GenomeIntervalTest
        extends TestCase {
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public GenomeIntervalTest(String testName) {
        super(testName);
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite() {
        return new TestSuite(GenomeIntervalTest.class);
    }

    public void testGenomeLocationSimple() {
        final String intervalString = "1,100,200,POSITIVE,INV,";
        final GenomeInterval genomeInterval = new GenomeInterval(intervalString);

        assertEquals("1", genomeInterval.getChromosome());
        assertEquals(100, genomeInterval.getStart());
        assertEquals(200, genomeInterval.getEnd());
        assertEquals(Strand.POSITIVE, genomeInterval.getStrand());
        assertEquals(MapBlock.BlockType.INV, genomeInterval.getFeature());
    }
}
