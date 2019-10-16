package com.bina.intervalTree;

import org.apache.log4j.Logger;

import java.util.Random;

/**
 * Simple interval, End-points inclusive, immutable
 *
 * @author johnmu
 */
public class SimpleInterval1D implements Comparable<Interval1D>, Interval1D {
    private final static Logger log = Logger.getLogger(SimpleInterval1D.class.getName());
    public final long left; // left endpoint, inclusive
    public final long right; // right endpoint, inclusive

    /**
     * Require left <= right
     *
     * @param left  Left end-point (inclusive)
     * @param right Right end-point (inclusive)
     */
    public SimpleInterval1D(long left, long right) {
        if ((left - 1 <= right) && ((((right / 2) + 1) - ((left / 2) - 1)) < (Long.MAX_VALUE / 2))) {
            this.left = left;
            this.right = right;
        } else {
            throw new RuntimeException("Illegal SimpleInterval1D (negative range or too large range): " + left + "-" + right);
        }
    }

    public SimpleInterval1D(Interval1D reg) {
        if (reg.length() < 0) {
            throw new RuntimeException("Illegal SimpleInterval1D (negative range or too large range): " + reg);
        }
        this.left = reg.getLeft();
        this.right = reg.getRight();
    }

    /**
     * Initialise with a random interval
     *
     * @param r
     */
    public SimpleInterval1D(Random r) {
        long left;
        long right;
        do {
            long a = r.nextLong();
            long b = r.nextLong();
            left = Math.min(a, b);
            right = Math.max(a, b);
        } while (right - left + 1 < 0);
        this.left = left;
        this.right = right;
    }

    public SimpleInterval1D(Random r, int lower, int upper) {
        long a = r.nextInt(upper - lower) + lower;
        long b = r.nextInt(upper - lower) + lower;
        this.left = Math.min(a, b);
        this.right = Math.max(a, b);
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
        return this.getRight() - this.getLeft() + 1;
    }

    @Override
    public long getLeft() {
        return left;
    }

    @Override
    public long getRight() {
        return right;
    }

    /**
     * Test if an interval intersects this one
     *
     * @param that Interval to test
     * @return True if it intersects
     */
    public boolean intersects(Interval1D that) {
        return intersects(that, 0, 0);
    }

    /**
     * does this interval reciprocal overlap intersect that one, with a particular
     * ratio. The ratio is rounded up. Require at least one base overlap.
     * <p/>
     * Can suffer from floating point precision issues for large numbers
     *
     * @param that            Interval to compare to
     * @param reciprocalRatio Ratio between 0 and 1, values outside this range will have undefined behaviour
     * @return True if it intersects, with the required reciprocal overlap
     */
    public boolean intersects(Interval1D that, double reciprocalRatio) {
        return intersects(that, reciprocalRatio, 0);
    }

    /**
     * Checks reciprocal overlap with wiggle, each possible shift within the wiggle is checked
     * <p/>
     * TODO this is really slow, better to do simple intersect then further filter in most cases
     *
     * @param that            Interval to compare with
     * @param reciprocalRatio Ratio between 0 and 1, values outside this range will have undefined behaviour
     * @param wiggle          Shift the interval between this range and for each shift check the overlap
     * @return True if overlaps with given criterion
     */
    public boolean intersects(Interval1D that, double reciprocalRatio, int wiggle) {
        /*
        //assume wiggle >= 0
        when that.left > this.left
                ------------------that
        --------------this
        min(wiggle, that.left - this.left) >= 0, this can be shifted to the right

        when that.left < this.left
        -----------that
                --------------this
        min(wiggle, that.left - this.left) == that.left - this.left
        0 >= max(-wiggle, that.left - this.left) >= -wiggle, this can be shifted to the left

        when wiggle == 0,
        no shift is allowed.
         */
        //guarantee that reciprocalRatio is not zero, which does not make sense for overlapping
        reciprocalRatio = Math.max(reciprocalRatio, ((double)1/Long.MAX_VALUE));
        long maxAllowedShift = Math.max(-wiggle, Math.min(wiggle, that.getLeft() - this.getLeft()));
        long maxOverlap = Math.min(this.getRight() + maxAllowedShift, that.getRight()) - Math.max(this.getLeft() + maxAllowedShift, that.getLeft()) + 1l;
        /*assumption: 0-length interval vs 0-length interval or non-zero-length vs non-zero-length
        not good for 1-vs-0 intervals
         */
        return maxOverlap >= Math.max(this.length(), that.length()) * reciprocalRatio;

    }

    /**
     * @param point point to test
     * @return True if interval contains the point
     */
    public boolean contains(final long point) {
        return (getLeft() <= point) && (point <= getRight());
    }

    /**
     * @param that Interval to compare to
     * @return zero if equal
     */
    public int compareTo(final Interval1D that) {
        if (this.getLeft() < that.getLeft())
            return -1;
        else if (this.getLeft() > that.getLeft())
            return +1;
        else if (this.getRight() < that.getRight())
            return -1;
        else if (this.getRight() > that.getRight())
            return +1;
        else
            return 0;
    }

    /**
     * Note: Not really the union, maybe this method name should be changed
     *
     * @param that Another interval
     * @return interval that encloses both current interval and the provided one
     */
    public SimpleInterval1D union(final Interval1D that) {
        return new SimpleInterval1D(Math.min(this.getLeft(), that.getLeft()), Math.max(this.getRight(), that.getRight()));
    }

    /**
     * @return The midpoint of the interval floor((right+left)/2)
     */
    public long getCenter() {
        return Math.round(Math.floor(((getRight() / 2.0) + (getLeft() / 2.0))));
    }

    public String toString() {
        return "[" + getLeft() + ", " + getRight() + "]";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SimpleInterval1D)) return false;

        SimpleInterval1D that = (SimpleInterval1D) o;

        if (left != that.left) return false;
        if (right != that.right) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = (int) (left ^ (left >>> 32));
        result = 31 * result + (int) (right ^ (right >>> 32));
        return result;
    }
}
