package com.bina.varsim.util;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by johnmu on 8/21/14.
 */
public abstract class GzFileParser<T> {
    protected BufferedReader bufferedReader = null;
    protected String line = "";

    protected InputStream decompressStream(final String fileName) throws IOException {
        return decompressStream(new File(fileName));
    }

    protected InputStream decompressStream(final File file) throws IOException {
        PushbackInputStream pb = new PushbackInputStream(new FileInputStream(file), 1024 * 1024); //we need a pushbackstream to look ahead
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
        return (line != null);
    }


    /**
     * @return the variant read, null if the line is not a variant line or there was an error parsing
     */
    public abstract T parseLine();


    protected void readLine() {
        line = null;
        try {
            line = bufferedReader.readLine();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        if (line == null) {
            try {
                bufferedReader.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            bufferedReader = null;
        }
    }

}






































