package com.binatechnologies.varsim.intervalTree;

import junit.framework.TestCase;

public class SimpleInterval1DTest extends TestCase {
    SimpleInterval1D A;
    SimpleInterval1D B;
    SimpleInterval1D C;

    public void setUp() throws Exception {
        super.setUp();
        A = new SimpleInterval1D(0, 10);
        B = new SimpleInterval1D(11, 20);
        C = new SimpleInterval1D(5, 15);
    }

    public void testLength() throws Exception {
        assertEquals(11, A.length());
        assertEquals(10, B.length());
        assertEquals(11, C.length());
    }

    public void testIntersects() throws Exception {

    }

    public void testIntersects1() throws Exception {

    }

    public void testIntersects2() throws Exception {

    }

    public void testContains() throws Exception {
        assertEquals(true, A.contains(0));
        assertEquals(true, A.contains(10));
        assertEquals(false, A.contains(-1));
        assertEquals(false, A.contains(11));
    }

    public void testCompareTo() throws Exception {
        assertEquals(0,A.compareTo(new SimpleInterval1D(0,10)));
        assertEquals(false,A.compareTo(B) == 0);
    }

    public void testUnion() throws Exception {

    }

    public void testGetCenter() throws Exception {

    }
}