package com.bina.varsim;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by johnmu on 8/21/14.
 */
public abstract class GzFileParser<T> {
    protected BufferedReader _br = null;
    protected String _line = "";

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
    public abstract T parseLine();


    protected void readLine() {
        _line = null;
        try {
            _line = _br.readLine();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        if (_line == null) {
            try {
                _br.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            _br = null;
        }
    }

}






































