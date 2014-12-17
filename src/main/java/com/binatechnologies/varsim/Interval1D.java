package com.binatechnologies.varsim;

import org.apache.log4j.Logger;

// end-points are inclusive, immutable
public class Interval1D implements Comparable<Interval1D> {
    private final static Logger log = Logger.getLogger(Interval1D.class.getName());
    public final long left; // left endpoint, inclusive
    public final long right; // right endpoint, inclusive

    // precondition: left <= right
    public Interval1D(long left, long right) {
        if (left <= right) {
            this.left = left;
            this.right = right;
        } else {
            throw new RuntimeException("Illegal interval: " + left + "-" + right);
        }
    }

    public static void main(String[] args) {
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

    public long length() {
        return this.right - this.left + 1;
    }

    // does this interval intersect that one?
    public boolean intersects(Interval1D that) {
        if (that.right < this.left)
            return false;
        if (this.right < that.left)
            return false;
        return true;
    }

    // does this interval reciprocal intersect that one, with a particular
    // ratio. The ratio is rounded up.
    public boolean intersects(Interval1D that, double ratio) {
        long len_this = (long) Math.ceil(this.length() * ratio);
        long len_that = (long) Math.ceil(that.length() * ratio);
        long overlap = Math.min(this.right, that.right)
                - Math.max(this.left, that.left) + 1;

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
        long len_this = (long) Math.ceil(this.length() * ratio);
        long len_that = (long) Math.ceil(that.length() * ratio);

        long max_overlap = 0;
        long right_lim = 0;
        long left_lim = 0;

        // right limit
        if (right < that.right) {
            right_lim = Math.min(that.right, right + wiggle);
        } else {
            right_lim = Math.max(that.right, right - wiggle);
        }
        left_lim = right_lim - length() + 1;

        long overlap = Math.min(right_lim, that.right)
                - Math.max(left_lim, that.left) + 1;

        max_overlap = Math.max(max_overlap, overlap);

        // left limit
        if (left < that.left) {
            left_lim = Math.min(that.left, left + wiggle);
        } else {
            left_lim = Math.max(that.left, left - wiggle);
        }
        right_lim = left_lim + length() - 1;

        overlap = Math.min(right_lim, that.right)
                - Math.max(left_lim, that.left) + 1;

        max_overlap = Math.max(max_overlap, overlap);

        if (max_overlap >= Math.max(len_this, len_that)) {
            return true;
        }
        return false;
    }

    /**
     * @param point point to test
     * @return True if interval contains the point
     */
    public final boolean contains(final long point) {
        return (left <= point) && (point <= right);
    }

    /**
     * @param that com.binatechnologies.varsim.Interval1D to compare to
     * @return zero if equal
     */
    public int compareTo(final Interval1D that) {
        if (this.left < that.left)
            return -1;
        else if (this.left > that.left)
            return +1;
        else if (this.right < that.right)
            return -1;
        else if (this.right > that.right)
            return +1;
        else
            return 0;
    }

    /**
     * @param inter Another interval
     * @return union of current interval with the provided one
     */
    public Interval1D union(final Interval1D inter) {
        return new Interval1D(Math.min(left, inter.left), Math.max(right, inter.right));
    }

    /**
     *
     * @return The midpoint of the interval floor((right+left)/2)
     */
    public final long getCenter(){
        return ((right+left)/2);
    }

    public String toString() {
        return "[" + left + ", " + right + "]";
    }

}
