package com.bina.varsim;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by johnmu on 8/21/14.
 */
public abstract class variantFileParser {
    public static final int X = 23;
    public static final int Y = 24;
    public static final int MT = 25;
    protected BufferedReader _br = null;
    protected String _line = "";

    /**
     * Converts hg19 or b37 format chromosomes to b37 format
     *
     * @param chr Chromosome name as a string
     * @return chromosome name as a string in b37 format
     */
    public static String stripChr(String chr) {
        if (chr.length() > 3 && chr.substring(0, 3).equalsIgnoreCase("chr")) {
            return chr.substring(3);
        }
        if (chr.equals("M")) {
            chr = "MT";
        }
        return chr;
    }

    /**
     * Converts the b37 chromosome name to an integer representation, only works for human
     *
     * @param chr chromosome name as a string in b37 format
     * @return chromosome name as a integer
     */
    public static int getChromIndex(String chr) {
        int ret = -1;
        chr = stripChr(chr);

        if (chr.equalsIgnoreCase("X"))
            ret = X;
        else if (chr.equalsIgnoreCase("Y"))
            ret = Y;
        else if (chr.equalsIgnoreCase("M") || chr.equalsIgnoreCase("MT"))
            ret = MT;
        else
            try {
                ret = Integer.parseInt(chr);
            } catch (NumberFormatException e) {
                System.err.println("Unknown human chromosome " + chr + ".");
            }

        return ret;
    }

    protected InputStream decompressStream(final String fileName) throws IOException {
        PushbackInputStream pb = new PushbackInputStream(new FileInputStream(fileName), 1024 * 1024); //we need a pushbackstream to look ahead
        byte[] signature = new byte[2];
        pb.read(signature); //read the signature
        pb.unread(signature); //push back the signature to the stream
        if (signature[0] == (byte) (GZIPInputStream.GZIP_MAGIC & 0xff) && signature[1] == (byte) ((GZIPInputStream.GZIP_MAGIC >> 8) & 0xff)) {
            return new GZIPInputStream(pb, 1024 * 1024);
        } else {
            return pb;
        }
    }

    /**
     * @return true if the file has more to read
     */
    public boolean hasMoreInput() {
        return (_line != null);
    }


    /**
     * @return the variant read, null if the line is not a variant line or there was an error parsing
     */
    public abstract Variant parseLine();


    protected void readLine() {
        _line = null;
        try {
            _line = _br.readLine();
        } catch (IOException ex) {
            System.err.println(ex);
        }

        if (_line == null) {
            try {
                _br.close();
            } catch (IOException ex) {
                System.err.println(ex);
            }
            _br = null;
        }
    }

}






































