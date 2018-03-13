package com.bina.varsim.util;

/**
 * Reads in a FASTA reference to memory and allows fast queries
 *
 * @author johnmu
 */

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.Sequence;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.*;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.LineIterator;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;


public class SimpleReference {
    private final static Logger log = Logger.getLogger(SimpleReference.class.getName());

    // chr_idx -> reference_string
    private final Map<ChrString, Sequence> data = new HashMap<>();
    private final Map<ChrString, IndexedFastaSequenceFile> dataSources = new HashMap<>();
    private long nonNCount = -1;
    private String referenceFileName = null;
    private boolean sequencesFullyLoaded = false;

    public SimpleReference() {
    }

    /**
     * @param filename FASTA file, should include the fai index
     */
    public SimpleReference(String filename) {
        addReference(filename);
        referenceFileName = filename;
    }

    public SimpleReference(final Collection<String> filenames) {
        if (filenames != null) {
            StringBuilder stringBuilder = new StringBuilder();
            for (final String filename : filenames) {
                addReference(filename);
                stringBuilder.append(filename);
            }
            referenceFileName = stringBuilder.toString();
        }
    }

    /**
     * This loads a FASTA file into the memory, can load multi-FASTA files, can call multiple times
     *
     * @param filename FASTA file with fai index
     */
    private void addReference(String filename) {
        File f = new File(filename);
        Path indexFile = ReferenceSequenceFileFactory.getFastaIndexFileName(Paths.get(filename));
        Path dictionaryFile = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(Paths.get(filename));
        try {
            if (!Files.exists(indexFile, new LinkOption[0])) {
                FastaSequenceIndexCreator.create(Paths.get(filename), false);
            }
            if (!Files.exists(dictionaryFile, new LinkOption[0])) {
                CreateSequenceDictionary sequenceDictionaryCreator = new CreateSequenceDictionary(filename);
            }
            IndexedFastaSequenceFile fa = new IndexedFastaSequenceFile(f);
            for (SAMSequenceRecord s : fa.getSequenceDictionary().getSequences()) {
                ChrString name = new ChrString(s.getSequenceName());
                if (!name.toString().equals(s.getSequenceName())) {
                    throw new IllegalArgumentException("Internal name " + name + " is different from name in file (" + filename + "): " + s.getSequenceName());
                }
                if (!data.containsKey(name)) {
                    data.put(name, null); //be lazy here
                    dataSources.put(name, fa);
                } else {
                    log.warn("Duplicate Key!");
                }
            }
        } catch (IOException e) {
            log.error(filename + " not found.");
            e.printStackTrace();
            System.exit(1);
        }
        sequencesFullyLoaded = false; //everytime we add new sequences, reset this flag.
    }

    /**
     * @param chr_name chromosome name as a class
     * @return Entire sequence of the chromosome
     */
    public Sequence getSequence(ChrString chr_name) {
        if (data.containsKey(chr_name)) {
            if (data.get(chr_name) == null) {
                ReferenceSequence sequence = dataSources.get(chr_name).getSequence(chr_name.toString());
                data.put(chr_name, new Sequence(chr_name.toString(), sequence.getBases(), sequence.length()));
            }
        } else {
            return null;
        }
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

        if (data.containsKey(chr_name)) {
            if (data.get(chr_name) == null) {
                return dataSources.get(chr_name).getSubsequenceAt(chr_name.toString(), loc, loc).getBases()[0];
            } else {
                return data.get(chr_name).byteAt(loc);
            }
        } else {
            return 0;
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

        if (data.containsKey(chr_name)) {
            if (data.get(chr_name) == null) {
                return dataSources.get(chr_name).getSubsequenceAt(chr_name.toString(), start_loc + 1, end_loc).getBases();
            } else {
                Sequence contig = data.get(chr_name);
                return contig.subSeq(start_loc, end_loc);
            }
        } else {
            log.error("Contig not found: " + chr_name);
            return null;
        }
    }

    /**
     * @param chr_name Chromosome name as a number
     * @return length of the specified chromosome
     */
    public int getRefLen(ChrString chr_name) {
        if (data.containsKey(chr_name)) {
            return dataSources.get(chr_name).getSequenceDictionary().getSequence(chr_name.toString()).getSequenceLength();
        } else {
            return 0;
        }
    }

    /**
     * @return Number of chromosomes loaded
     */
    public int getNumContigs() {
        return data.keySet().size();
    }

    private void loadAllSequences() {
        if (!sequencesFullyLoaded) {
            for (ChrString contig : data.keySet()) {
                if (data.get(contig) == null) {
                    ReferenceSequence contigSequence = dataSources.get(contig).getSequence(contig.toString());
                    data.put(contig, new Sequence(contig.toString(), contigSequence.getBases(), contigSequence.length()));
                }
            }
            sequencesFullyLoaded = true;
        }
    }
    /**
     * This is computed lazily, so the first call to this will be slow
     * @return number of non-N bases in the sequence
     */
    public long getNumNonNBases(){
        if (nonNCount == -1) {
            loadAllSequences();
            for (Sequence sequence : data.values()) {
                nonNCount += sequence.getNumNonNBases();
            }
        }
        return nonNCount;
    }

    public long getNumNonNBases(final File regions) throws IOException {
        loadAllSequences();
        long count = 0;

        final FeatureCodec<BEDFeature, LineIterator> bedCodec = new BEDCodec(BEDCodec.StartOffset.ONE);
        final LineIterator lineIterator = new AsciiLineReaderIterator(new AsciiLineReader(new FileInputStream(regions)));

        while (lineIterator.hasNext()) {
            final BEDFeature bedFeature = bedCodec.decode(lineIterator);
            count += data.get(new ChrString(bedFeature.getContig())).getNumNonNBases(bedFeature.getStart(), bedFeature.getEnd());
        }
        return count;
    }

    /**
     * @return Chromosomes currently loaded
     */
    public Set<ChrString> keySet() {
        return data.keySet();
    }

    /**
     * @return referenceFileName name(s) read
     */
    public String getReferenceFileName() {
        return referenceFileName;
    }
}
