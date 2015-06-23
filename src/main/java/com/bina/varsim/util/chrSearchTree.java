package com.bina.varsim.util;


/**
 * Chromosome search tree, essentially a map from chromosome name to an interval tree
 *
 * @author johnmu
 */

import com.bina.intervaltree.Interval1D;
import com.bina.intervaltree.IntervalTree;
import com.bina.intervaltree.SimpleInterval1D;
import com.bina.varsim.types.ChrString;
import org.apache.log4j.Logger;

import java.util.HashMap;


public class chrSearchTree<K extends Interval1D> {
    private final static Logger log = Logger.getLogger(chrSearchTree.class.getName());

    HashMap<ChrString, IntervalTree<K>> data;
    boolean allow_duplicates = false; // this will allow duplicate intervals

    /**
     * Constructors
     */

    public chrSearchTree() {
        data = new HashMap<>();
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
        for (IntervalTree<K> val : data.values()) {
            total += val.size();
        }
        return total;
    }

    public long maxDepth() {
        long maxval = 0;
        for (IntervalTree<K> val : data.values()) {
            maxval = Math.max(maxval, val.maxDepth());
        }
        return maxval;
    }

    /**
     * Adding functions
     */

    /**
     * Add an interval to the appropriate search tree
     *
     * @param chrname Chromosome name as a string
     * @param key     Entry to be inserted
     */
    public void put(ChrString chrname, K key) {
        IntervalTree<K> out = data.get(chrname);
        if (out == null) {
            IntervalTree<K> contents = new IntervalTree<>();
            contents.add(key);
            log.info("Added chromosome: " + chrname);
            data.put(chrname, contents);
        } else {
            out.add(key);
        }
    }

    /**
     * Search functions
     */


    public Iterable<K> getOverlaps(ChrString chrname, Interval1D key) {
        return getOverlaps(chrname, key, 0, 0);
    }

    public Iterable<K> getOverlaps(ChrString chrname, Interval1D key, double reciprocalRatio) {
        return getOverlaps(chrname, key, reciprocalRatio, 0);
    }

    /**
     * @param chrname         Chromosome name as a string
     * @param key             Interval to be searched for
     * @param reciprocalRatio minimum reciprocal overlap required, 0 means minimum one position overlap
     * @param wiggle          Amount the interval can be shifted
     * @return All the values corresponding to the intervals overlapping the specified interval
     */
    public Iterable<K> getOverlaps(ChrString chrname, Interval1D key, double reciprocalRatio, int wiggle) {
        IntervalTree<K> out = data.get(chrname);
        if (out == null) {
            return null;
        } else {
            return out.getOverlaps(key, reciprocalRatio, wiggle);
        }
    }

    public boolean contains(ChrString chrname, Interval1D key) {
        return contains(chrname, key, 0, 0);
    }

    public boolean contains(ChrString chrname, Interval1D key, double reciprocalRatio) {
        return contains(chrname, key, reciprocalRatio, 0);
    }

    /**
     * @param chrname         Chromosome name as a string
     * @param key             Interval to be searched for
     * @param reciprocalRatio minimum reciprocal overlap required, 0 means minimum one position overlap
     * @param wiggle          Amount the interval can be shifted
     * @return Whether the specified interval is overlapped at all
     */
    public boolean contains(ChrString chrname, Interval1D key, double reciprocalRatio, int wiggle) {
        IntervalTree<K> out = data.get(chrname);
        if (out == null) {
            return false;
        } else {
            return out.contains(key, reciprocalRatio, wiggle);
        }
    }

    /**
     * @param chrname         Chromosome name as a string
     * @param start           start location of interval (inclusive)
     * @param end             end location of interval (inclusive)
     * @param reciprocalRatio minimum reciprocal overlap required, 0 means minimum one position overlap
     * @return Whether the specified interval is overlapped at all
     */
    public boolean contains(ChrString chrname, int start, int end, double reciprocalRatio) {
        return contains(chrname, new SimpleInterval1D(start, end), reciprocalRatio);
    }

}
