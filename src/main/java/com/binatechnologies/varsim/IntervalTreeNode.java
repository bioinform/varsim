package com.binatechnologies.varsim;

import java.util.ArrayList;

/**
 *
 * The current design only allows searching for intervals, searching for a point will
 * require creating an interval of length 1.
 *
 * Created by johnmu on 12/6/14.
 */
public class IntervalTreeNode<Key extends Interval1D> {
    IntervalTreeNode<Key> left = null;
    IntervalTreeNode<Key> right = null;
    protected ArrayList<Key> center = new ArrayList<Key>();
    final long centerCut; // intervals overlapping this point are placed in center

    IntervalTreeNode(long centerCut) {
        this.centerCut = centerCut;
    }

    IntervalTreeNode(final IntervalTreeNode<Key> left, final IntervalTreeNode<Key> right, final ArrayList<Key> center, final long centerCut) {
        this.left = left;
        this.right = right;
        for (Key k : center) {
            this.center.add(k);
        }
        this.centerCut = centerCut;
    }

    IntervalTreeNode(final IntervalTreeNode<Key> that) {
        this.left = that.left;
        this.right = that.right;
        for (Key k : that.center) {
            this.center.add(k);
        }
        this.centerCut = that.centerCut;
    }

    /**
     * @param k key to try and add
     * @return 0 if it overlaps center, -1 for left, 1 for right
     */
    int addKey(final Key k) {
        int checkVal = checkKey(k);
        if(checkVal == 0){
            center.add(k);
        }
        return checkVal;
    }

    /**
     * @param k key to check
     * @return 0 if it overlaps center, -1 for completely on left, 1 for completely on right
     */
    int checkKey(final Key k) {
        // check if the key overlaps the center cut
        if (overlapCenter(k)) {
            return 0;
        }
        if (k.right < centerCut) {
            return -1;
        } else {
            return 1;
        }
    }

    final boolean overlapCenter(final Key k){
        return k.contains(centerCut);
    }

    /**
     * Goes through the center list and finds all the overlaps with the key
     *
     * TODO This can be made more efficient, but maybe not necessary for most sane datasets
     *
     * @param k
     * @return
     */
    final ArrayList<Key> getOverlaps(final Key k){
        ArrayList<Key> retVal = new ArrayList<Key>();
        for(Key centerKey : center){
            if(centerKey.intersects(k)){
                retVal.add(centerKey);
            }
        }
        return retVal;
    }

    final ArrayList<Key> getOverlaps(final Key k, double ratio, int wiggle){
        ArrayList<Key> retVal = new ArrayList<Key>();
        for(Key centerKey : center){
            if(centerKey.intersects(k,ratio, wiggle)){
                retVal.add(centerKey);
            }
        }
        return retVal;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        IntervalTreeNode that = (IntervalTreeNode) o;

        if (centerCut != that.centerCut) return false;
        if (center != null ? !center.equals(that.center) : that.center != null) return false;
        if (left != null ? !left.equals(that.left) : that.left != null) return false;
        if (right != null ? !right.equals(that.right) : that.right != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = left != null ? left.hashCode() : 0;
        result = 31 * result + (right != null ? right.hashCode() : 0);
        result = 31 * result + (center != null ? center.hashCode() : 0);
        result = 31 * result + (int) (centerCut ^ (centerCut >>> 32));
        return result;
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
}
