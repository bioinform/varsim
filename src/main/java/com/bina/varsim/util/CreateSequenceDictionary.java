package com.bina.varsim.util;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Md5CalculatingOutputStream;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * Modified by guoy28 on 3/12/18.
 */
public class CreateSequenceDictionary {
  /*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
  /**
   * Create a SAM/BAM file from a fasta containing reference sequence. The output SAM file contains a header but no
   * SAMRecords, and the header contains only sequence records.
   */
    private static final Log logger = Log.getInstance(CreateSequenceDictionary.class);

    public File OUTPUT;
    public String GENOME_ASSEMBLY;
    public String URI;
    public String SPECIES;
    public boolean TRUNCATE_NAMES_AT_WHITESPACE = true;
    public int NUM_SEQUENCES = Integer.MAX_VALUE;
    private final MessageDigest md5;
    public File REFERENCE;
    public File REFERENCE_SEQUENCE;
    public boolean CREATE_MD5_FILE = false;

  public CreateSequenceDictionary(String FastaSequenceFileName) {
    try {
      md5 = MessageDigest.getInstance("MD5");
    } catch (NoSuchAlgorithmException e) {
      throw new RuntimeException("MD5 algorithm not found", e);
    }
    REFERENCE_SEQUENCE = new File(FastaSequenceFileName);
    URI = "file:" + REFERENCE_SEQUENCE.getAbsolutePath();
    OUTPUT = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(REFERENCE_SEQUENCE);
    logger.info("Output dictionary will be written in ", OUTPUT);
    doWork();
  }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    protected int doWork() {
      if (OUTPUT.exists()) {
        throw new RuntimeException(OUTPUT.getAbsolutePath() +
                " already exists.  Delete this file and try again, or specify a different output file.");
      }

      // SortingCollection is used to check uniqueness of sequence names
      final SortingCollection<String> sequenceNames = makeSortingCollection();
      try (BufferedWriter writer = makeWriter()) {
        final ReferenceSequenceFile refSeqFile = ReferenceSequenceFileFactory.
                getReferenceSequenceFile(REFERENCE_SEQUENCE, TRUNCATE_NAMES_AT_WHITESPACE);
        SAMSequenceDictionaryCodec samDictCodec = new SAMSequenceDictionaryCodec(writer);

        samDictCodec.encodeHeaderLine(false);
        // read reference sequence one by one and write its metadata
        for (ReferenceSequence refSeq = refSeqFile.nextSequence(); refSeq != null; refSeq = refSeqFile.nextSequence()) {
          final SAMSequenceRecord samSequenceRecord = makeSequenceRecord(refSeq);
          samDictCodec.encodeSequenceRecord(samSequenceRecord);
          sequenceNames.add(refSeq.getName());
        }
      } catch (FileNotFoundException e) {
        throw new RuntimeException("File " + OUTPUT.getAbsolutePath() + " not found");
      } catch (IOException e) {
        throw new RuntimeException("Can't write to or close output file " + OUTPUT.getAbsolutePath());
      }

      // check uniqueness of sequences names
      final CloseableIterator<String> iterator = sequenceNames.iterator();

      if(!iterator.hasNext()) return 0;

      String current = iterator.next();
      while (iterator.hasNext()) {
        final String next = iterator.next();
        if (current.equals(next)) {
          OUTPUT.delete();
          throw new RuntimeException("Sequence name " + current +
                  " appears more than once in reference file");
        }
        current = next;
      }
      return 0;
    }

    private BufferedWriter makeWriter() throws FileNotFoundException {
      return new BufferedWriter(
              new AsciiWriter(this.CREATE_MD5_FILE ?
                      new Md5CalculatingOutputStream(
                              new FileOutputStream(OUTPUT, false),
                              new File(OUTPUT.getAbsolutePath() + ".md5")
                      )
                      : new FileOutputStream(OUTPUT)
              )
      );
    }

    /**
     * Create one SAMSequenceRecord from a single fasta sequence
     */
    private SAMSequenceRecord makeSequenceRecord(final ReferenceSequence refSeq) {
      final SAMSequenceRecord ret = new SAMSequenceRecord(refSeq.getName(), refSeq.length());

      // Compute MD5 of upcased bases
      final byte[] bases = refSeq.getBases();
      for (int i = 0; i < bases.length; ++i) {
        bases[i] = StringUtil.toUpperCase(bases[i]);
      }

      ret.setAttribute(SAMSequenceRecord.MD5_TAG, md5Hash(bases));
      if (GENOME_ASSEMBLY != null) {
        ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, GENOME_ASSEMBLY);
      }
      ret.setAttribute(SAMSequenceRecord.URI_TAG, URI);
      if (SPECIES != null) {
        ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, SPECIES);
      }
      return ret;
    }

    private String md5Hash(final byte[] bytes) {
      md5.reset();
      md5.update(bytes);
      String s = new BigInteger(1, md5.digest()).toString(16);
      if (s.length() != 32) {
        final String zeros = "00000000000000000000000000000000";
        s = zeros.substring(0, 32 - s.length()) + s;
      }
      return s;
    }

    private SortingCollection<String> makeSortingCollection() {
      final String name = getClass().getSimpleName();
      final File tmpDir = IOUtil.createTempDir(name, null);
      tmpDir.deleteOnExit();
      // 256 byte for one name, and 1/10 part of all memory for this, rough estimate
      long maxNamesInRam = Runtime.getRuntime().maxMemory() / 256 / 10;
      return SortingCollection.newInstance(
              String.class,
              new StringCodec(),
              String::compareTo,
              (int) Math.min(maxNamesInRam, Integer.MAX_VALUE),
              tmpDir
      );
    }

    private static class StringCodec implements SortingCollection.Codec<String> {
      private DataInputStream dis;
      private DataOutputStream dos;

      public StringCodec clone() {
        return new StringCodec();
      }

      public void setOutputStream(final OutputStream os) {
        dos = new DataOutputStream(os);
      }

      public void setInputStream(final InputStream is) {
        dis = new DataInputStream(is);
      }

      public void encode(final String str) {
        try {
          dos.writeUTF(str);
        } catch (IOException e) {
          throw new RuntimeIOException(e);
        }
      }

      public String decode() {
        try {
          return dis.readUTF();
        } catch (EOFException e) {
          return null;
        } catch (IOException e) {
          throw new RuntimeException("Exception reading sequence name from temporary file.", e);
        }
      }
    }
}
