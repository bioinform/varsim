package com.bina.varsim.fastqLiftover.readers;

import java.io.IOException;
import java.io.LineNumberReader;

public class FastqReader {
    private LineNumberReader fastqBr;

    public FastqReader(LineNumberReader fastqBr) {
        this.fastqBr = fastqBr;
    }

    public String[] getNextFastqEntry() throws IOException {
        String[] fastqEntry = new String[4];
        for (int i = 0; i < 4; i++) {
            fastqEntry[i] = fastqBr.readLine();
            if (fastqEntry[i] == null) {
                return null;
            }
        }
        return fastqEntry;
    }

    public int getLineNumber() {
        return fastqBr.getLineNumber();
    }
}