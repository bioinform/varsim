package com.binatechnologies.varsim.intervalTree;

//import org.apache.log4j.Logger;

/**
 * Simple interval, End-points inclusive, immutable
 *
 * @author johnmu
 */
public class SimpleInterval1D implements Comparable<Interval1D>,Interval1D {
    //private final static Logger log = Logger.getLogger(SimpleInterval1D.class.getName());
    public final long left; // left endpoint, inclusive
    public final long right; // right endpoint, inclusive

    /**
     * Require left <= right
     *
     * @param left  Left end-point (inclusive)
     * @param right Right end-point (inclusive)
     */
    public SimpleInterval1D(long left, long right) {
        if (left <= right) {
            this.left = left;
            this.right = right;
        } else {
            throw new RuntimeException("Illegal interval: " + left + "-" + right);
        }
    }

    public SimpleInterval1D(Interval1D reg){
        this.left = reg.left;
        this.right = reg.right;
    }

    /**
     * Some simple test cases
     *
     * @param args
     */
    public static void main(String[] args) {
        // Test the interval tree
        SimpleInterval1D a = new SimpleInterval1D(15, 30);
        SimpleInterval1D b = new SimpleInterval1D(20, 30);
        SimpleInterval1D c = new SimpleInterval1D(15, 20);
        SimpleInterval1D d = new SimpleInterval1D(31, 35);

        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("c = " + c);
        System.out.println("d = " + d);

        System.out.println("a == b " + a.compareTo(b));
        System.out.println("a == c " + a.compareTo(c));
        System.out.println("a == a " + a.compareTo(a));
        System.out.println("a == d " + a.compareTo(d));

        System.out.println("b intersects a = " + b.intersects(a));
        System.out.println("a intersects b = " + a.intersects(b));
        System.out.println("a intersects c = " + a.intersects(c));
        System.out.println("c intersects a = " + c.intersects(a));
        System.out.println("a intersects d = " + a.intersects(d));
        System.out.println("b intersects c = " + b.intersects(c));
        System.out.println("b intersects d = " + b.intersects(d));
        System.out.println("c intersects d = " + c.intersects(d));

        double ratio = 0.5;

        System.out.println("ratio b intersects a = " + b.intersects(a, ratio));
        System.out.println("ratio a intersects b = " + a.intersects(b, ratio));
        System.out.println("ratio a intersects c = " + a.intersects(c, ratio));
        System.out.println("ratio c intersects a = " + c.intersects(a, ratio));
        System.out.println("ratio a intersects d = " + a.intersects(d, ratio));
        System.out.println("ratio b intersects c = " + b.intersects(c, ratio));
        System.out.println("ratio b intersects d = " + b.intersects(d, ratio));
        System.out.println("ratio c intersects d = " + c.intersects(d, ratio));

        int wiggle = 5;

        System.out.println("wiggle ratio b intersects a = " + b.intersects(a, ratio, wiggle));
        System.out.println("wiggle ratio a intersects b = " + a.intersects(b, ratio, wiggle));
        System.out.println("wiggle ratio a intersects c = " + a.intersects(c, ratio, wiggle));
        System.out.println("wiggle ratio c intersects a = " + c.intersects(a, ratio, wiggle));
        System.out.println("wiggle ratio a intersects d = " + a.intersects(d, ratio, wiggle));
        System.out.println("wiggle ratio b intersects c = " + b.intersects(c, ratio, wiggle));
        System.out.println("wiggle ratio c intersects b = " + c.intersects(b, ratio, wiggle));
        System.out.println("wiggle ratio b intersects d = " + b.intersects(d, ratio, wiggle));
        System.out.println("wiggle ratio c intersects d = " + c.intersects(d, ratio, wiggle));
    }

    /**
     * @return Length of the interval, end-points are inclusive
     */
    public long length() {
        return this.right - this.left + 1;
    }

    /**
     * Test if an interval intersects this one
     *
     * @param that Interval to test
     * @return True if it intersects
     */
    public boolean intersects(Interval1D that) {
        if (that.right < this.left)
            return false;
        if (this.right < that.left)
            return false;
        return true;
    }

    /**
     * does this interval reciprocal overlap intersect that one, with a particular
     * ratio. The ratio is rounded up. Require at least one base overlap.
     *
     * @param that            Interval to compare to
     * @param reciprocalRatio Ratio between 0 and 1, values outside this range will have undefined behaviour
     * @return True if it intersects, with the required reciprocal overlap
     */
    public boolean intersects(Interval1D that, double reciprocalRatio) {
        if (reciprocalRatio == 0) {
            return intersects(that);
        }
        // Note: The max may be able to be removed in this case, left to be safe
        long thisLen = (long) Math.max(Math.ceil(this.length() * reciprocalRatio), 1);
        long thatLen = (long) Math.max(Math.ceil(that.length() * reciprocalRatio), 1);
        long overlap = Math.min(this.right, that.right)
                - Math.max(this.left, that.left) + 1;
        if (overlap >= Math.max(thisLen, thatLen)) {
            return true;
        }
        return false;
    }

    /**
     * Checks reciprocal overlap with wiggle, each possible shift within the wiggle is checked
     *
     * @param that            Interval to compare with
     * @param reciprocalRatio Ratio between 0 and 1, values outside this range will have undefined behaviour
     * @param wiggle          Shift the interval between this range and for each shift check the overlap
     * @return True if overlaps with given criterion
     */
    public boolean intersects(Interval1D that, double reciprocalRatio, int wiggle) {
        if (wiggle == 0) {
            return intersects(that, reciprocalRatio);
        }
        long len_this = (long) Math.max(Math.ceil(this.length() * reciprocalRatio), 1);
        long len_that = (long) Math.max(Math.ceil(that.length() * reciprocalRatio), 1);

        long maxOverlap = 0;
        long rightLim = 0;
        long leftLim = 0;

        // right limit
        if (right < that.right) {
            rightLim = Math.min(that.right, right + wiggle);
        } else {
            rightLim = Math.max(that.right, right - wiggle);
        }
        leftLim = rightLim - length() + 1;

        long overlap = Math.min(rightLim, that.right)
                - Math.max(leftLim, that.left) + 1;

        maxOverlap = Math.max(maxOverlap, overlap);

        // left limit
        if (left < that.left) {
            leftLim = Math.min(that.left, left + wiggle);
        } else {
            leftLim = Math.max(that.left, left - wiggle);
        }
        rightLim = leftLim + length() - 1;

        overlap = Math.min(rightLim, that.right)
                - Math.max(leftLim, that.left) + 1;

        maxOverlap = Math.max(maxOverlap, overlap);

        if (maxOverlap >= Math.max(len_this, len_that)) {
            return true;
        }
        return false;
    }

    /**
     * @param point point to test
     * @return True if interval contains the point
     */
    public boolean contains(final long point) {
        return (left <= point) && (point <= right);
    }

    /**
     * @param that Interval to compare to
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
     * @param that Another interval
     * @return union of current interval with the provided one
     */
    public SimpleInterval1D union(final Interval1D that) {
        return new SimpleInterval1D(Math.min(this.left, that.left), Math.max(this.right, that.right));
    }

    /**
     * @return The midpoint of the interval floor((right+left)/2)
     */
    public long getCenter() {
        return ((right + left) / 2);
    }

    public String toString() {
        return "[" + left + ", " + right + "]";
    }

}
