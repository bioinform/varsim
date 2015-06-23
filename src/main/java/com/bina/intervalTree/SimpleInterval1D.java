package com.bina.intervaltree;

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
        if ((left <= right) && ((((right / 2) + 1) - ((left / 2) - 1)) < (Long.MAX_VALUE / 2))) {
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
        if (that.getRight() < this.getLeft())
            return false;
        if (this.getRight() < that.getLeft())
            return false;
        return true;
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
        if (reciprocalRatio == 0) {
            return intersects(that);
        }
        // Note: The max may be able to be removed in this case, left to be safe
        long thisLen = reciprocalRatio == 1.0 ? this.length() : (long) Math.max(Math.ceil(this.length() * reciprocalRatio), 1l);
        long thatLen = reciprocalRatio == 1.0 ? that.length() : (long) Math.max(Math.ceil(that.length() * reciprocalRatio), 1l);
        long overlap = Math.min(this.getRight(), that.getRight())
                - Math.max(this.getLeft(), that.getLeft()) + 1l;
        if (overlap >= Math.max(thisLen, thatLen)) {
            return true;
        }
        return false;
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
        if (wiggle == 0) {
            return intersects(that, reciprocalRatio);
        }
        long len_this = reciprocalRatio == 1.0 ? this.length() : (long) Math.max(Math.ceil(this.length() * reciprocalRatio), 1l);
        long len_that = reciprocalRatio == 1.0 ? that.length() : (long) Math.max(Math.ceil(that.length() * reciprocalRatio), 1l);

        long maxOverlap = 0;
        long rightLim;
        long leftLim;

        // right limit
        if (getRight() < that.getRight()) {
            rightLim = Math.min(that.getRight(), getRight() + wiggle);
        } else {
            rightLim = Math.max(that.getRight(), getRight() - wiggle);
        }
        leftLim = rightLim - length() + 1l;

        long overlap = Math.min(rightLim, that.getRight())
                - Math.max(leftLim, that.getLeft()) + 1l;

        maxOverlap = Math.max(maxOverlap, overlap);

        // left limit
        if (getLeft() < that.getLeft()) {
            leftLim = Math.min(that.getLeft(), getLeft() + wiggle);
        } else {
            leftLim = Math.max(that.getLeft(), getLeft() - wiggle);
        }
        rightLim = leftLim + length() - 1l;

        overlap = Math.min(rightLim, that.getRight())
                - Math.max(leftLim, that.getLeft()) + 1l;

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
        return (((getRight() / 2) + (getLeft() / 2)));
    }

    public String toString() {
        return "[" + getLeft() + ", " + getRight() + "]";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

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
