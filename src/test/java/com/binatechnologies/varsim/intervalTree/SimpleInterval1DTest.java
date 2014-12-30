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
        assertEquals(false, A.intersects(B));
        assertEquals(false, B.intersects(A));
        assertEquals(true, A.intersects(C));
        assertEquals(true, C.intersects(A));
    }

    public void testIntersects1() throws Exception {
        assertEquals(false, A.intersects(B,0.5));
        assertEquals(false, B.intersects(A,0.5));
        assertEquals(false, A.intersects(C,0.55));
        assertEquals(false, C.intersects(A,0.55));
    }

    public void testIntersects2() throws Exception {
        assertEquals(true, A.intersects(B,0.5,6));
        assertEquals(true, B.intersects(A,0.5,6));
        assertEquals(true, A.intersects(C,0.55,1));
        assertEquals(true, C.intersects(A,0.55,1));
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
        assertEquals(new SimpleInterval1D(0,20),A.union(B));
    }

    public void testGetCenter() throws Exception {
        assertEquals(5,A.getCenter());
    }
}