package com.binatechnologies.seqalto.varsim;

import org.apache.log4j.Logger;

// Modified from http://algs4.cs.princeton.edu/93intersection/

/**
 * **********************************************************************
 * Compilation: javac Interval1D.java Execution: java Interval1D
 * <p/>
 * Interval ADT for integer coordinates.
 * <p/>
 * ***********************************************************************
 */

// end-points are inclusive
// this class seems to be immutable
public class Interval1D implements Comparable<Interval1D> {
    private final static Logger log = Logger.getLogger(Interval1D.class.getName());
    public final int low; // left endpoint
    public final int high; // right endpoint

    // precondition: left <= right
    public Interval1D(int left, int right) {
        if (left <= right) {
            this.low = left;
            this.high = right;
        } else {
            throw new RuntimeException("Illegal interval: " + left + "-" + right);
        }
    }

    public static void main(String[] args) {
        // TODO Auto-generated method stub

        // Test the interval tree
        Interval1D a = new Interval1D(15, 30);
        Interval1D b = new Interval1D(20, 30);
        Interval1D c = new Interval1D(15, 20);
        Interval1D d = new Interval1D(20, 25);

        System.err.println("a == b " + a.compareTo(b));
        System.err.println("a == c " + a.compareTo(c));
        System.err.println("a == a " + a.compareTo(a));
        System.err.println("a == d " + a.compareTo(d));

        System.err.println("a = " + a);
        System.err.println("b = " + b);
        System.err.println("c = " + c);
        System.err.println("d = " + d);

        System.err.println("b intersects a = " + b.intersects(a));
        System.err.println("a intersects b = " + a.intersects(b));
        System.err.println("a intersects c = " + a.intersects(c));
        System.err.println("c intersects a = " + c.intersects(a));
        System.err.println("a intersects d = " + a.intersects(d));
        System.err.println("b intersects c = " + b.intersects(c));
        System.err.println("b intersects d = " + b.intersects(d));
        System.err.println("c intersects d = " + c.intersects(d));

        double ratio = 0.5;

        System.err.println("ratio b intersects a = " + b.intersects(a, ratio));
        System.err.println("ratio a intersects b = " + a.intersects(b, ratio));
        System.err.println("ratio a intersects c = " + a.intersects(c, ratio));
        System.err.println("ratio c intersects a = " + c.intersects(a, ratio));
        System.err.println("ratio a intersects d = " + a.intersects(d, ratio));
        System.err.println("ratio b intersects c = " + b.intersects(c, ratio));
        System.err.println("ratio b intersects d = " + b.intersects(d, ratio));
        System.err.println("ratio c intersects d = " + c.intersects(d, ratio));

        int wiggle = 5;

        System.err.println("wiggle ratio b intersects a = " + b.intersects(a, ratio, wiggle));
        System.err.println("wiggle ratio a intersects b = " + a.intersects(b, ratio, wiggle));
        System.err.println("wiggle ratio a intersects c = " + a.intersects(c, ratio, wiggle));
        System.err.println("wiggle ratio c intersects a = " + c.intersects(a, ratio, wiggle));
        System.err.println("wiggle ratio a intersects d = " + a.intersects(d, ratio, wiggle));
        System.err.println("wiggle ratio b intersects c = " + b.intersects(c, ratio, wiggle));
        System.err.println("wiggle ratio c intersects b = " + c.intersects(b, ratio, wiggle));
        System.err.println("wiggle ratio b intersects d = " + b.intersects(d, ratio, wiggle));
        System.err.println("wiggle ratio c intersects d = " + c.intersects(d, ratio, wiggle));
    }

    public int length() {
        return this.high - this.low + 1;
    }

    // does this interval intersect that one?
    public boolean intersects(Interval1D that) {
        if (that.high < this.low)
            return false;
        if (this.high < that.low)
            return false;
        return true;
    }

    // does this interval reciprocal intersect that one, with a particular
    // ratio? The ratio is rounded up.
    public boolean intersects(Interval1D that, double ratio) {
        int len_this = (int) Math.ceil(this.length() * ratio);
        int len_that = (int) Math.ceil(that.length() * ratio);
        int overlap = Math.min(this.high, that.high)
                - Math.max(this.low, that.low) + 1;

        if (overlap >= Math.max(len_this, len_that)) {
            return true;
        }
        return false;
    }

    /*
    does this interval reciprocal intersect that one, with a particular
    ratio and wiggle range? The ratio is rounded up.
    That is, for any shift forward or back of the wiggle, is there ratio overlap
    */
    public boolean intersects(Interval1D that, double ratio, int wiggle) {
        int len_this = (int) Math.ceil(this.length() * ratio);
        int len_that = (int) Math.ceil(that.length() * ratio);

        int max_overlap = 0;
        int right_lim = 0;
        int left_lim = 0;

        log.trace("Comparing " + toString() + " with " + that);

        // right limit
        if (high < that.high) {
            right_lim = Math.min(that.high, high + wiggle);
        } else {
            right_lim = Math.max(that.high, high - wiggle);
        }
        left_lim = right_lim - length() + 1;

        //System.err.println("right: " + left_lim + "-" + right_lim);

        int overlap = Math.min(right_lim, that.high)
                - Math.max(left_lim, that.low) + 1;

        //System.err.println("roverlap: " + overlap);

        max_overlap = Math.max(max_overlap, overlap);

        // left limit
        if (low < that.low) {
            left_lim = Math.min(that.low, low + wiggle);
        } else {
            left_lim = Math.max(that.low, low - wiggle);
        }
        right_lim = left_lim + length() - 1;

        //System.err.println("left: " + left_lim + "-" + right_lim);

        overlap = Math.min(right_lim, that.high)
                - Math.max(left_lim, that.low) + 1;

        //System.err.println("loverlap: " + overlap);

        max_overlap = Math.max(max_overlap, overlap);

        //System.err.println("max_overlap: " + max_overlap);
        //System.err.println("len_this: " + len_this);
        //System.err.println("len_that: " + len_that);

        if (max_overlap >= Math.max(len_this, len_that)) {
            return true;
        }
        return false;
    }

    // does this interval a intersect b?
    public boolean contains(int x) {
        return (low <= x) && (x <= high);
    }

    public int compareTo(Interval1D that) {
        if (this.low < that.low)
            return -1;
        else if (this.low > that.low)
            return +1;
        else if (this.high < that.high)
            return -1;
        else if (this.high > that.high)
            return +1;
        else
            return 0;
    }

    // creates union of two intervals
    public Interval1D union(final Interval1D inter) {
        return new Interval1D(Math.min(low, inter.low), Math.max(high, inter.high));
    }

    public String toString() {
        return "[" + low + ", " + high + "]";
    }

}
