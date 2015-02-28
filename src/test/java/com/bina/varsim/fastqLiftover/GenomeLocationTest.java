package com.bina.varsim.fastqLiftover;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for GenomeLocation
 */
public class GenomeLocationTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public GenomeLocationTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( GenomeLocationTest.class );
    }

    public void testGenomeLocationSimple()
    {
        final String locationString = "1-100-INV-";
        final GenomeLocation genomeLocation = new GenomeLocation(locationString);

        assertEquals("1", genomeLocation.chromosome);
        assertEquals(100, genomeLocation.location);
        assertEquals(-1, genomeLocation.read_location1);
        assertEquals("INV", genomeLocation.feature);
        assertEquals(1, genomeLocation.direction);
    }

    public void testGenomeLocationSlash()
    {
        final String locationString = "1-100/-INV-";
        final GenomeLocation genomeLocation = new GenomeLocation(locationString);

        assertEquals("1", genomeLocation.chromosome);
        assertEquals(100, genomeLocation.location);
        assertEquals(-1, genomeLocation.read_location1);
        assertEquals("INV", genomeLocation.feature);
        assertEquals(1, genomeLocation.direction);
    }

    public void testGenomeLocationReadLoc()
    {
        final String locationString = "1-100/7-INV-";
        final GenomeLocation genomeLocation = new GenomeLocation(locationString);

        assertEquals("1", genomeLocation.chromosome);
        assertEquals(100, genomeLocation.location);
        assertEquals(7, genomeLocation.read_location1);
        assertEquals("INV", genomeLocation.feature);
        assertEquals(1, genomeLocation.direction);
    }
}
