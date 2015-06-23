package com.bina.intervaltree;

import java.util.ArrayList;

/**
 * The current design only allows searching for intervals, searching for a point will
 * require creating an interval of length 1.
 * <p/>
 * Created by johnmu on 12/6/14.
 *
 * @author johnmu
 */
public class IntervalTreeNode<Key extends Interval1D> {
    protected final long centerCut; // intervals overlapping this point are placed in center
    protected int balanceFactor = 0; // positive for heavier to the right
    private IntervalTreeNode<Key> left;
    private IntervalTreeNode<Key> right;
    private ArrayList<Key> center = new ArrayList<>();

    public IntervalTreeNode(long centerCut) {
        this.centerCut = centerCut;
    }

    /**
     * Create new node from a key, the centerCut will be set to the middle of the Interval
     *
     * @param k
     */
    public IntervalTreeNode(Key k) {
        this(k.getCenter());
        center.add(k);
    }

    public IntervalTreeNode(final IntervalTreeNode<Key> left, final IntervalTreeNode<Key> right,
                            final ArrayList<Key> center, final long centerCut) {
        this.left = left;
        this.right = right;
        for (Key k : center) {
            this.center.add(k);
        }
        this.centerCut = centerCut;
    }

    public IntervalTreeNode(final IntervalTreeNode<Key> that) {
        this.left = that.left;
        this.right = that.right;
        for (Key k : that.center) {
            this.center.add(k);
        }
        this.centerCut = that.centerCut;
    }

    /**
     * @return positive values indicate more depth on the right
     */
    public int getBalanceFactor() {
        return balanceFactor;
    }

    public void setBalanceFactor(int balanceFactor) {
        this.balanceFactor = balanceFactor;
    }

    /**
     * This is for adding to the right
     */
    public void incBalanceFactor() {
        balanceFactor++;
    }

    public void addBalanceFactor(int val) {
        balanceFactor += val;
    }

    /**
     * This is for adding to the left
     */
    public void decBalanceFactor() {
        balanceFactor--;
    }


    public ArrayList<Key> getCenter() {
        return center;
    }

    void setCenter(ArrayList<Key> center) {
        this.center = center;
    }

    /**
     * @param k key to try and add
     * @return 0 if it overlaps center, -1 for left, 1 for right
     */
    public int addKey(final Key k) {
        int checkVal = checkKey(k);
        if (checkVal == 0) {
            center.add(k);
        }
        return checkVal;
    }

    /**
     * @param k key to check
     * @return 0 if it overlaps center, -1 for completely on left, 1 for completely on right
     */
    public int checkKey(final Interval1D k) {
        // check if the key overlaps the center cut
        if (overlapCenter(k)) {
            return 0;
        }
        if (k.getRight() < centerCut) {
            return -1;
        } else {
            return 1;
        }
    }

    public final boolean overlapCenter(final Interval1D k) {
        return k.contains(centerCut);
    }

    /**
     * Goes through the center list and finds all the overlaps with the key
     * <p/>
     * TODO This can be made more efficient, but maybe not necessary for most sane datasets
     *
     * @param k
     * @return All intervals that overlap at least one base with the provided one
     */
    public final ArrayList<Key> getOverlaps(final Interval1D k) {
        return getOverlaps(k, 0, 0);
    }

    /**
     * @param k key to search for
     * @return true if the node contains the key
     */
    public final boolean contains(final Interval1D k) {
        return contains(k, 0, 0);
    }

    /**
     * @param k
     * @param reciprocalRatio
     * @return
     */
    public final ArrayList<Key> getOverlaps(final Interval1D k, double reciprocalRatio) {
        return getOverlaps(k, reciprocalRatio, 0);
    }

    /**
     * @param k
     * @param reciprocalRatio
     * @return
     */
    public final boolean contains(final Interval1D k, double reciprocalRatio) {
        return contains(k, reciprocalRatio, 0);
    }

    /**
     * Goes through the center list and finds all the overlaps (reciprocal with wiggle) with the key
     *
     * @param k               Key to search for
     * @param reciprocalRatio Reciprocal overlap ratio
     * @param wiggle          Try to shift the interval within this wiggle
     * @return All intervals that overlap with the provided one matching the criteria
     */
    public final ArrayList<Key> getOverlaps(final Interval1D k, double reciprocalRatio, int wiggle) {
        ArrayList<Key> retVal = new ArrayList<>();
        for (Key centerKey : center) {
            if (centerKey.intersects(k, reciprocalRatio, wiggle)) {
                retVal.add(centerKey);
            }
        }
        return retVal;
    }

    /**
     * @param k               Key to search for
     * @param reciprocalRatio Reciprocal overlap ratio
     * @param wiggle          Try to shift the interval within this wiggle
     * @return True if node contains the key with specified criteria
     */
    public final boolean contains(final Interval1D k, double reciprocalRatio, int wiggle) {
        for (Interval1D centerKey : center) {
            if (centerKey.intersects(k, reciprocalRatio, wiggle)) {
                return true;
            }
        }
        return false;
    }

    public IntervalTreeNode<Key> getLeft() {
        return left;
    }

    public void setLeft(IntervalTreeNode<Key> left) {
        this.left = left;
    }

    public IntervalTreeNode<Key> getRight() {
        return right;
    }

    public void setRight(IntervalTreeNode<Key> right) {
        this.right = right;
    }

    /**
     * @param child sets child that matches this, either left or right
     * @param val   sets to this val
     */
    public void setChild(IntervalTreeNode<Key> child, IntervalTreeNode<Key> val) {
        if (getRight() == child) {
            setRight(val);
        } else if (getLeft() == child) {
            setLeft(val);
        } else {
            throw new RuntimeException("Tried to set a child that did to match");
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("cut: ");
        sb.append(centerCut);
        sb.append(" bf: ");
        sb.append(balanceFactor);
        sb.append(" -- ");
        for (Interval1D val : center) {
            sb.append(val).append(",");
        }
        return sb.toString();
    }
}
