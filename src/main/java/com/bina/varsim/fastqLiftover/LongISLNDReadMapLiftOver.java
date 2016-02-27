package com.bina.varsim.fastqLiftover;

import com.bina.varsim.fastqLiftover.readers.ARTPairedFastqAlnReader;
import com.bina.varsim.fastqLiftover.readers.DWGSIMPairedFastqReader;
import com.bina.varsim.fastqLiftover.readers.PBSIMFastqReader;
import com.bina.varsim.fastqLiftover.readers.PairedFastqReader;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MapBlocks;
import com.bina.varsim.fastqLiftover.types.SimulatedReadPair;
import com.bina.varsim.readers.longislnd.LongISLNDReadAlignmentMap;
import com.bina.varsim.types.ReadMapRecord;
import org.apache.commons.codec.binary.Base64;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class LongISLNDReadMapLiftOver {
    private final static Logger log = Logger.getLogger(LongISLNDReadMapLiftOver.class.getName());
    String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
    @Option(name = "-map", usage = "Map file", metaVar = "file")
    private File mapFile;
    @Option(name = "-longislnd", usage = "Read map file from LongISLND", metaVar = "file")
    private File longislnd;
    @Option(name = "-out", usage = "Output file", metaVar = "file")
    private File outFile;
    @Option(name = "-compress", usage = "Use GZIP compression")
    private boolean compress;

    public static void main(String[] args) throws IOException {
        new LongISLNDReadMapLiftOver().run(args);
    }

    public PrintStream getOutStream(final File outFile, boolean compress) throws IOException {
        PrintStream ps = System.out;
        if (outFile != null) {
            if (compress) {
                ps = new PrintStream(new GZIPOutputStream(new FileOutputStream(outFile)) {{
                    def.setLevel(Deflater.BEST_SPEED);
                }});
            } else {
                ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFile)));
            }
        }
        return ps;
    }

    public void run(String[] args) throws IOException {
        CmdLineParser parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java LongISLNDReadMapLiftOver [options...] arguments...");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println();
            return;
        }

        long tStart = System.currentTimeMillis();
        MapBlocks mapBlocks = new MapBlocks(mapFile);

        doLiftOver(mapBlocks, longislnd, getOutStream(outFile, compress));

        log.info("Conversion took " + (System.currentTimeMillis() - tStart) / 1e3 + " seconds.");
    }

    public void doLiftOver(final MapBlocks mapBlocks, final File longislnd, final PrintStream ps) throws IOException {
        log.info("Reading read map file " + longislnd.getName());
        final Collection<ReadMapRecord> readMapRecords = new LongISLNDReadAlignmentMap(longislnd).getReadAlignmentMap().values();
        int readCount = 0;
        for (final ReadMapRecord readMapRecord : readMapRecords) {
            ps.println(mapBlocks.liftOverReadMapRecord(readMapRecord));
            if ((readCount % 100000) == 0) {
                log.info(readCount + " reads processed");
            }
            readCount++;
        }
        ps.close();
        log.info(readCount + "reads processed");
    }
}
