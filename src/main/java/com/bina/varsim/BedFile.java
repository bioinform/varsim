package com.bina.varsim;

import com.bina.varsim.intervalTree.SimpleInterval1D;

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
    chrSearchTree<SimpleInterval1D> bedST; // the interval search tree for bed file

    /**
     * Reads the BED file into a search tree
     *
     * @param filename BED file
     */
    public BedFile(String filename) {
        bedST = new chrSearchTree<>();
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
            if (line.length() == 0 || line.charAt(0) == '#') {
                // comment line
                continue;
            }

            // TODO replace this with apache-commons for speed
            String[] ll = line.split("\t");

            String chr_name = ll[0];
            int start;
            int end;
            try {
                start = Integer.parseInt(ll[1]);
                end = Integer.parseInt(ll[2]);
            } catch (NumberFormatException e) {
                throw new RuntimeException("Malformed BED line (nfe): " + line);
            }

            bedST.put(new ChrString(chr_name), new SimpleInterval1D(start, end - 1));
        }
    }

    /**
     * Check whether an interval in contained within the interval tree, don't allow for wiggle
     *
     * @param chr chromosome name as a string
     * @param start   start location (inclusive)
     * @param end     end location (inclusive)
     * @return
     */
    public boolean contains(ChrString chr, int start, int end) {
        return bedST.contains(chr, new SimpleInterval1D(start, end));
    }

    /**
     * Check whether an interval in contained within the interval tree, don't allow for wiggle
     *
     * @param chr  chromosome name as a string
     * @param interval Interval to search (inclusive)
     * @return
     */
    public boolean contains(ChrString chr, SimpleInterval1D interval) {
        return bedST.contains(chr, interval);
    }

    /**
     * Checks whether either endpoint of the interval is in the BED file
     *
     * @param chr
     * @param interval
     * @return
     */
    public boolean containsEitherEndpoint(ChrString chr, SimpleInterval1D interval) {
        if (interval.getLeft() == interval.getRight()) {
            return bedST.contains(chr, interval);
        } else {
            return bedST.contains(chr, new SimpleInterval1D(interval.getLeft(), interval.getLeft())) |
                    bedST.contains(chr, new SimpleInterval1D(interval.getRight(), interval.getRight()));
        }
    }

}
