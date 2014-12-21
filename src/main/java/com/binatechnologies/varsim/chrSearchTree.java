package com.binatechnologies.varsim;


/**
 * Chromosome search tree, essentially a map from chromosome name to an interval tree
 * @author johnmu
 */

import com.binatechnologies.varsim.intervalTree.Interval1D;
import com.binatechnologies.varsim.intervalTree.IntervalTree;
import com.binatechnologies.varsim.intervalTree.SimpleInterval1D;

import java.util.HashMap;


public class chrSearchTree<Key extends Interval1D> {

    HashMap<String, IntervalTree<Key>> data;
    boolean allow_duplicates = false; // this will allow duplicate intervals

    /**
     * Constructors
     */

    public chrSearchTree() {
        data = new HashMap<String, IntervalTree<Key>>();
    }


    /**
     * @param allow_duplicates Allow duplicate intervals
     */
    public chrSearchTree(boolean allow_duplicates) {
        this();
        this.allow_duplicates = allow_duplicates;
    }

    /**
     * @return Total number of intervals stored in all chromosomes
     */
    public long size() {
        long total = 0;
        for (IntervalTree<Key> val : data.values()) {
            total += val.size();
        }
        return total;
    }

    /**
     * Adding functions
     */

    /**
     * Add an interval to the appropriate search tree
     *
     * @param chrname  Chromosome name as a string
     * @param key      Entry to be inserted
     */
    public void put(String chrname, Key key) {
        IntervalTree<Key> out = data.get(chrname);
        if (out == null) {
            IntervalTree<Key> contents = new IntervalTree<Key>();
            contents.add(key);
            System.err.println("Put: " + chrname);
            data.put(chrname, contents);
        } else {
            out.add(key);
        }
    }

    /**
     * Search functions
     */


    public Iterable<Key> getOverlaps(String chrname, Interval1D key) {
        return getOverlaps(chrname,key,0,0);
    }

    public Iterable<Key> getOverlaps(String chrname, Interval1D key, double reciprocalRatio) {
        return getOverlaps(chrname,key,reciprocalRatio,0);
    }

    /**
     * @param chrname  Chromosome name as a string
     * @param key Interval to be searched for
     * @param reciprocalRatio    minimum reciprocal overlap required, 0 means minimum one position overlap
     * @param wiggle Amount the interval can be shifted
     * @return All the values corresponding to the intervals overlapping the specified interval
     */
    public Iterable<Key> getOverlaps(String chrname, Interval1D key, double reciprocalRatio, int wiggle) {
        IntervalTree<Key> out = data.get(chrname);
        if (out == null) {
            return null;
        } else {
            return out.getOverlaps(key, reciprocalRatio, wiggle);
        }
    }

    public boolean contains(String chrname, Interval1D key) {
        return contains(chrname,key,0,0);
    }

    public boolean contains(String chrname, Interval1D key, double reciprocalRatio) {
        return contains(chrname,key,reciprocalRatio,0);
    }

    /**
     * @param chrname  Chromosome name as a string
     * @param key Interval to be searched for
     * @param reciprocalRatio    minimum reciprocal overlap required, 0 means minimum one position overlap
     * @param wiggle Amount the interval can be shifted
     * @return Whether the specified interval is overlapped at all
     */
    public boolean contains(String chrname, Interval1D key, double reciprocalRatio, int wiggle) {
        IntervalTree<Key> out = data.get(chrname);
        if (out == null) {
            return false;
        } else {
            return out.contains(key, reciprocalRatio, wiggle);
        }
    }

    /**
     * @param chrname Chromosome name as a string
     * @param start   start location of interval (inclusive)
     * @param end     end location of interval (inclusive)
     * @param reciprocalRatio   minimum reciprocal overlap required, 0 means minimum one position overlap
     * @return Whether the specified interval is overlapped at all
     */
    public boolean contains(String chrname, int start, int end, double reciprocalRatio) {
        return contains(chrname, new SimpleInterval1D(start, end), reciprocalRatio);
    }

}
