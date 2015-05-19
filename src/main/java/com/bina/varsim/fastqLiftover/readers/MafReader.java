package com.bina.varsim.fastqLiftover.readers;

import com.bina.varsim.fastqLiftover.types.MafRecord;

import java.util.Iterator;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.IOException;

public class MafReader implements Iterator<MafRecord> {
    private String lastLine_ = null;
    private MafRecord nextRecord_ = null;
    private BufferedReader brMaf_ = null;

    public MafReader(final BufferedReader brMaf) {
        brMaf_ = brMaf;
        nextRecord_ = readNextEntry();
    }

    /**
     * return the next complete MafRecord record, null otherwise
     *
     * @return  a new instance of MafRecord if a complete entry is retrieved, null otherwise
     */
    @Override
    public MafRecord next() {
        MafRecord ret = nextRecord_;
        nextRecord_ = readNextEntry();
        return ret;
    }

    @Override
    public boolean hasNext() {
        return nextRecord_ != null;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * attempt to read a MAF record from brMaf_
     *
     * @return  a new instance of MafRecord if a complete entry is retrieved, null otherwise
     */
    private MafRecord readNextEntry() {
        try {
            String buffer = lastLine_ != null ? lastLine_ : brMaf_.readLine();
            lastLine_ = null;
            for( ; buffer != null ; buffer = brMaf_.readLine() ){
                buffer = buffer.trim();
                if( buffer.length() > 0 && buffer.charAt(0) != '#') break;
            }
            if (buffer == null) return null;
            if (buffer.charAt(0) != 'a') throw new RuntimeException("unexpected MAF line: "+buffer);
            String header = buffer;
            ArrayList<String> entries = new ArrayList<String>(2);
            for (buffer=brMaf_.readLine() ; buffer != null ; buffer=brMaf_.readLine()) {
                buffer = buffer.trim();
                if (buffer.length() > 0) {
                    if (buffer.charAt(0) == 's' ) {
                      entries.add(buffer);
                    }
                    else if(buffer.charAt(0) == 'a') {
                      break;
                    }
                    else if(buffer.charAt(0) != '#') {
                        throw new RuntimeException("unexpected MAF line: "+buffer);
                    }
                }
            }
            lastLine_ = buffer;
            return new MafRecord(header,entries);
        }
        catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }
}
