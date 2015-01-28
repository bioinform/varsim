package com.bina.varsim.fastqLiftover;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for SimulatedRead
 */
public class SimulatedReadTest 
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public SimulatedReadTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( SimulatedReadTest.class );
    }

    public void testSimulatedReadNameParsingEmpty()
    {
        final String readName = ":::::::::::::2/2";
        final SimulatedRead testRead = new SimulatedRead(readName);

        assertEquals(0, testRead.locs1.size());
        assertEquals(0, testRead.locs2.size());
        assertEquals(0, testRead.origLocs1.size());
        assertEquals(0, testRead.origLocs2.size());
        assertEquals(0, testRead.random1);
        assertEquals(0, testRead.random2);
        assertEquals(0, testRead.seqErrors1);
        assertEquals(0, testRead.snps1);
        assertEquals(0, testRead.indels1);
        assertEquals(0, testRead.seqErrors2);
        assertEquals(0, testRead.snps2);
        assertEquals(0, testRead.indels2);
        assertEquals("", testRead.readId);
        assertEquals(2, testRead.laneId);
    }

    public void testSimulatedReadNameParsingNormal()
    {
        final String readName = "1-223167089--,1-223167099-:1-223167425--:1_paternal-223817897-:1_paternal-223818223--:::0:::1:::e9f:2/1";
        final SimulatedRead testRead = new SimulatedRead(readName);

        assertEquals(2, testRead.locs1.size());
        assertEquals(new GenomeLocation("1", 223167089, 1), testRead.locs1.get(0));
        assertEquals(new GenomeLocation("1", 223167099, 0), testRead.locs1.get(1));

        assertEquals(1, testRead.locs2.size());
        assertEquals(new GenomeLocation("1", 223167425, 1), testRead.locs2.get(0));

        assertEquals(1, testRead.origLocs1.size());
        assertEquals(new GenomeLocation("1_paternal", 223817897, 0), testRead.origLocs1.get(0));

        assertEquals(1, testRead.origLocs2.size());
        assertEquals(new GenomeLocation("1_paternal", 223818223, 1), testRead.origLocs2.get(0));

        assertEquals(0, testRead.random1);
        assertEquals(0, testRead.random2);
        assertEquals(0, testRead.seqErrors1);
        assertEquals(0, testRead.snps1);
        assertEquals(0, testRead.indels1);
        assertEquals(1, testRead.seqErrors2);
        assertEquals(0, testRead.snps2);
        assertEquals(0, testRead.indels2);
        assertEquals("e9f", testRead.readId);
        assertEquals(2, testRead.laneId);
    }
}
