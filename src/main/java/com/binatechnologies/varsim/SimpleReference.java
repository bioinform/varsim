package com.binatechnologies.varsim;

/**
 *  Reads in a FASTA reference to memory and allows fast queries
 *  @author johnmu
 */

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;


public class SimpleReference {
    private final static Logger log = Logger.getLogger(SimpleReference.class.getName());

    // chr_idx -> reference_string
    HashMap<Integer, Sequence> data;

    // list of chromosomes already loaded
    ArrayList<Integer> chr_list;

    SimpleReference() {
        data = new HashMap<Integer, Sequence>();
        chr_list = new ArrayList<Integer>();
    }

    /**
     * @param filename FASTA file, should include the fai index
     */
    SimpleReference(String filename) {
        this();
        addReference(filename);
    }

    /**
     * This loads a FASTA file into the memory, can load multi-FASTA files, can call multiple times
     *
     * @param filename FASTA file with fai index
     */
    void addReference(String filename) {
        File f = new File(filename);
        FastaSequenceFile fa = new FastaSequenceFile(f, true);
        ReferenceSequence seq = fa.nextSequence();
        while (seq != null) {
            int name = variantFileParser.getChromIndex(seq.getName());
            byte[] seq_bytes = seq.getBases();

            if (seq_bytes == null) {
                System.err.println("Contig error: " + seq.getName());
            } else {
                if (name > 0) {
                    System.err.println("Read ref: " + seq.getName() + ":"
                            + name);
                    if (!data.containsKey(name)) {
                        Sequence contig = new Sequence(seq.getName(),
                                seq_bytes, seq_bytes.length);
                        data.put(name, contig);
                        chr_list.add(name);
                    } else {
                        System.err.println("Duplicate Key!");
                    }
                }
            }
            seq = fa.nextSequence();
        }
    }

    /**
     * @param chr_name chromosome name as a string
     * @return Entire sequence of the chromosome
     */
    public Sequence getSequence(String chr_name) {
        return getSequence(variantFileParser.getChromIndex(chr_name));
    }

    /**
     * @param chr_idx chromosome name as a number
     * @return Entire sequence of the chromosome
     */
    public Sequence getSequence(int chr_idx) {
        return data.get(chr_idx);
    }

    /**
     * 1-based
     *
     * @param chr_name chromosome name as a string
     * @param loc      location in chromosome 1-based
     * @return base at the specified position
     */
    byte byteAt(String chr_name, int loc) {
        return byteAt(variantFileParser.getChromIndex(chr_name), loc);
    }

    /**
     * 1-based
     *
     * @param chr_idx chromosome name as a number
     * @param loc     location in chromosome 1-based
     * @return base at the specified position
     */
    byte byteAt(int chr_idx, int loc) {
        if (loc < 1) {
            return 0;
        }

        Sequence contig = data.get(chr_idx);
        if (contig == null) {
            return 0;
        } else {
            return contig.byteAt(loc);
        }
    }

    char charAt(int chr_idx, int loc){
        return (char)byteAt(chr_idx,loc);
    }

    char charAt(String chr_name, int loc){
        return (char)byteAt(chr_name,loc);
    }

    // 1-based, inclusive start, exclusive end

    /**
     * @param chr_name  chromosome name as a string
     * @param start_loc inclusive
     * @param end_loc   exclusive
     * @return an array of bytes that is the sequence of bases
     */
    byte[] byteRange(String chr_name, int start_loc, int end_loc) {
        return byteRange(variantFileParser.getChromIndex(chr_name), start_loc, end_loc);
    }

    /**
     * @param chr_idx   chromosome name as a number
     * @param start_loc inclusive
     * @param end_loc   exclusive
     * @return an array of bytes that is the sequence of bases
     */
    byte[] byteRange(int chr_idx, int start_loc, int end_loc) {
        if (start_loc < 1 || end_loc < 1 || end_loc < start_loc) {
            System.err.println("byteRange: Invalide range");
            return null;
        }

        Sequence contig = data.get(chr_idx);
        if (contig == null) {
            log.error("Contig not found: " + chr_idx);
            return null;
        } else {
            return contig.subSeq(start_loc, end_loc);
        }
    }

    /**
     * @param chr_idx Chromosome name as a number
     * @return length of the specified chromosome
     */
    int getRefLen(int chr_idx) {
        Sequence contig = data.get(chr_idx);
        if (contig == null) {
            return 0;
        } else {
            return contig.length();
        }
    }

    /**
     * @return Number of chromosomes loaded
     */
    int getNumContigs() {
        return chr_list.size();
    }

    /**
     * This is used to iterate through the loaded chromosomes... there is probably a better solution
     *
     * @param inner_idx idx in internal array
     * @return actual chromosome idx as a number
     */
    int nameAt(int inner_idx) {
        return chr_list.get(inner_idx);
    }
}
