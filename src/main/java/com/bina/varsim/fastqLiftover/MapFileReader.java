package com.bina.varsim.fastqLiftover;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class MapFileReader {
    private Scanner scanner;

    public MapFileReader(final File file) throws IOException {
        this.scanner = new Scanner(file);
    }

    public MapBlock getNext() throws IOException, IllegalArgumentException {
        MapBlock mapBlock = null;
        if (scanner.hasNext()) {
            final int size = scanner.nextInt();
            final String srcChr = scanner.next();
            final int srcLocation = scanner.nextInt();
            final String dstChr = scanner.next();
            final int dstLocation = scanner.nextInt();
            final String direction = scanner.next();
            final String featureType = scanner.next();
            final String featureName = scanner.next();

            if (size <= 0) {
                throw new IllegalArgumentException("Encountered invalid size " + size);
            }
            mapBlock = new MapBlock(size, srcChr, srcLocation, dstChr, dstLocation, direction, featureType, featureName);
        }
        return mapBlock;
    }
}
