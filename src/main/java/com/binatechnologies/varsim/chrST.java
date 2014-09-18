package com.binatechnologies.varsim;


/**
 * Chromosome search tree, essentially a map from chromosome name to an interval tree
 * @author johnmu
 */

import java.util.HashMap;


public class chrST<Value> {

    HashMap<String, IntervalST<Value>> data;
    boolean allow_duplicates = false; // this will allow duplicate intervals

    /**
     * Constructors
     */

    public chrST() {
        data = new HashMap<String, IntervalST<Value>>();
    }


    /**
     * @param allow_duplicates Allow duplicate intervals
     */
    public chrST(boolean allow_duplicates) {
        this();
        this.allow_duplicates = allow_duplicates;
    }

    /**
     * @return Total number of intervals stored in all chromosomes
     */
    public int size() {
        int total = 0;
        for (IntervalST<Value> val : data.values()) {
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
     * @param interval interval to add (inclusive)
     * @param value    Value to be associated with the interval
     */
    public void put(String chrname, Interval1D interval, Value value) {
        IntervalST<Value> out = data.get(chrname);
        if (out == null) {
            IntervalST<Value> contents = new IntervalST<Value>(allow_duplicates);
            contents.put(interval, value);
            System.err.println("Put: " + chrname);
            data.put(chrname, contents);
        } else {
            out.put(interval, value);
        }
    }

    /**
     * Search functions
     */

    /**
     * @param chrname  Chromosome name as a string
     * @param interval Interval to be searched for
     * @param ratio    minimum reciprocal overlap required, 0 means minimum one position overlap
     * @return All intervals overlapping the specified interval
     */
    public Iterable<Interval1D> searchAll(String chrname, Interval1D interval,
                                          double ratio) {
        IntervalST<Value> out = data.get(chrname);
        if (out == null) {
            // System.err.println("chrname not_found: " + chrname);
            return null;
        } else {
            // System.err.println("Found chrname: " + chrname);
            return out.searchAll(interval, ratio);
        }
    }

    /**
     * @param chrname  Chromosome name as a string
     * @param interval Interval to be searched for
     * @param ratio    minimum reciprocal overlap required, 0 means minimum one position overlap
     * @return All the values corresponding to the intervals overlapping the specified interval
     */
    public Iterable<Value> getAll(String chrname, Interval1D interval,
                                  double ratio) {
        IntervalST<Value> out = data.get(chrname);
        if (out == null) {
            // System.err.println("chrname not_found: " + chrname);
            return null;
        } else {
            // System.err.println("Found chrname: " + chrname);
            return out.getAll(interval, ratio);
        }
    }

    /**
     * @param chrname  Chromosome name as a string
     * @param interval Interval to be searched for
     * @param ratio    minimum reciprocal overlap required, 0 means minimum one position overlap
     * @return Whether the specified interval is overlapped at all
     */
    public boolean contains(String chrname, Interval1D interval, double ratio) {
        IntervalST<Value> out = data.get(chrname);
        if (out == null) {
            // System.err.println("chrname not_found: " + chrname);
            return false;
        } else {
            // System.err.println("Found chrname: " + chrname);
            return (out.search(interval, ratio) != null);
        }
    }

    /**
     * @param chrname Chromosome name as a string
     * @param start   start location of interval (inclusive)
     * @param end     end location of interval (inclusive)
     * @param ratio   minimum reciprocal overlap required, 0 means minimum one position overlap
     * @return Whether the specified interval is overlapped at all
     */
    public boolean contains(String chrname, int start, int end, double ratio) {
        return contains(chrname, new Interval1D(start, end), ratio);
    }

}
