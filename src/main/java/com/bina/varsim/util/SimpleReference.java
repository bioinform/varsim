package com.bina.varsim.util;

/**
 * Reads in a FASTA reference to memory and allows fast queries
 *
 * @author johnmu
 */

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.Sequence;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.HashMap;
import java.util.Set;


public class SimpleReference {
    private final static Logger log = Logger.getLogger(SimpleReference.class.getName());

    // chr_idx -> reference_string
    private final HashMap<ChrString, Sequence> data;

    public SimpleReference() {
        data = new HashMap<>();
    }

    /**
     * @param filename FASTA file, should include the fai index
     */
    public SimpleReference(String filename) {
        this();
        addReference(filename);
    }

    /**
     * This loads a FASTA file into the memory, can load multi-FASTA files, can call multiple times
     *
     * @param filename FASTA file with fai index
     */
    public void addReference(String filename) {
        File f = new File(filename);
        FastaSequenceFile fa = new FastaSequenceFile(f, true);
        ReferenceSequence seq = fa.nextSequence();
        while (seq != null) {
            ChrString name = new ChrString(seq.getName());
            byte[] seq_bytes = seq.getBases();

            if (seq_bytes == null) {
                log.error("Contig error: " + seq.getName());
            } else {
                log.info("Read ref: " + seq.getName() + ":" + name);
                if (!data.containsKey(name)) {
                    Sequence contig = new Sequence(seq.getName(),
                            seq_bytes, seq_bytes.length);
                    data.put(name, contig);
                } else {
                    log.warn("Duplicate Key!");
                }
            }
            seq = fa.nextSequence();
        }
    }

    /**
     * @param chr_name chromosome name as a class
     * @return Entire sequence of the chromosome
     */
    public Sequence getSequence(ChrString chr_name) {
        return data.get(chr_name);
    }

    /**
     * 1-based
     *
     * @param chr_name chromosome name as a class
     * @param loc      location in chromosome 1-based
     * @return base at the specified position
     */
    public byte byteAt(ChrString chr_name, int loc) {
        if (loc < 1) {
            return 0;
        }

        Sequence contig = data.get(chr_name);
        if (contig == null) {
            return 0;
        } else {
            return contig.byteAt(loc);
        }
    }

    public char charAt(ChrString chr_name, int loc) {
        return (char) byteAt(chr_name, loc);
    }

    // 1-based, inclusive start, exclusive end

    /**
     * @param chr_name  chromosome name as a number
     * @param start_loc inclusive
     * @param end_loc   exclusive
     * @return an array of bytes that is the sequence of bases
     */
    public byte[] byteRange(ChrString chr_name, int start_loc, int end_loc) {
        if (start_loc < 1 || end_loc < 1 || end_loc < start_loc) {
            log.error("byteRange: Invalid range");
            return null;
        }

        Sequence contig = data.get(chr_name);
        if (contig == null) {
            log.error("Contig not found: " + chr_name);
            return null;
        } else {
            return contig.subSeq(start_loc, end_loc);
        }
    }

    /**
     * @param chr_name Chromosome name as a number
     * @return length of the specified chromosome
     */
    public int getRefLen(ChrString chr_name) {
        Sequence contig = data.get(chr_name);
        if (contig == null) {
            return 0;
        } else {
            return contig.length();
        }
    }

    /**
     * @return Number of chromosomes loaded
     */
    public int getNumContigs() {
        return data.keySet().size();
    }

    /**
     * @return Chromosomes currently loaded
     */
    public Set<ChrString> keySet() {
        return data.keySet();
    }
}
