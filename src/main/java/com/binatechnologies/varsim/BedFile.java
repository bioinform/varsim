package com.binatechnologies.varsim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Reads in a BED file and allows testing of regions
 * Remember that a BED file is 0-based
 *
 * TODO this current implementation ignores the other bed columns other than chr,start,end
 * TODO implement the other types of overlap
 */

/**
 * @author johnmu
 */
public class BedFile {
    String _filename; // file name of the BED file
    chrSearchTree<Interval1D> bedST; // the interval search tree for bed file

    /**
     * Reads the BED file into a search tree
     *
     * @param filename BED file
     */
    public BedFile(String filename) {
        bedST = new chrSearchTree<Interval1D>();
        _filename = filename;
        try {
            readBedFile(new File(_filename));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Reads the bed file into a search tree
     *
     * @param f BED file
     * @throws IOException
     */
    private void readBedFile(File f) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.charAt(0) == '#') {
                // comment line
                continue;
            }

            String[] ll = line.split("\t");

            String chr_name = ll[0];
            int start = Integer.parseInt(ll[1]);
            int end = Integer.parseInt(ll[2]);

            bedST.put(chr_name, new Interval1D(start, end));
        }
    }

    /**
     * Check whether an interval in contained within the interval tree, don't allow for wiggle
     *
     * @param chrname chromosome name as a string
     * @param start   start location (inclusive)
     * @param end     end location (inclusive)
     * @return
     */
    public boolean contains(String chrname, int start, int end) {
        return bedST.contains(chrname, new Interval1D(start, end));
    }

    /**
     * Check whether an interval in contained within the interval tree, don't allow for wiggle
     *
     * @param chrname  chromosome name as a string
     * @param interval Interval to search (inclusive)
     * @return
     */
    public boolean contains(String chrname, Interval1D interval) {
        return bedST.contains(chrname, interval);
    }

}
