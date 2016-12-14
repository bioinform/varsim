package com.bina.varsim.fastqLiftover;

import com.bina.varsim.VarSimTool;
import com.bina.varsim.fastqLiftover.types.MapBlocks;
import com.bina.varsim.readers.longislnd.LongISLNDReadAlignmentMap;
import com.bina.varsim.types.ReadMapRecord;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.Collection;
import java.util.List;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

public class LongISLNDReadMapLiftOver extends VarSimTool {
    private final static Logger log = Logger.getLogger(LongISLNDReadMapLiftOver.class.getName());
    @Option(name = "-map", usage = "Map file", metaVar = "file", required = true)
    private File mapFile;
    @Option(name = "-longislnd", usage = "Read map files from LongISLND", metaVar = "file", required = true)
    private List<String> longislnd;
    @Option(name = "-out", usage = "Output file", metaVar = "file", required = true)
    private File outFile;
    @Option(name = "-compress", usage = "Use GZIP compression")
    private boolean compress;
    @Option(name = "-minIntervalLength", usage = "Minimum interval length in liftover")
    private int minIntervalLength = MapBlocks.MIN_LENGTH_INTERVAL;

    public LongISLNDReadMapLiftOver(final String command, final String description) {
        super(command, description);
    }

    public static void main(String[] args) throws IOException {
        new LongISLNDReadMapLiftOver("", "").run(args);
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
        if (!parseArguments(args)) {
            return;
        }

        long tStart = System.currentTimeMillis();
        MapBlocks mapBlocks = new MapBlocks(mapFile);

        doLiftOver(mapBlocks, longislnd, getOutStream(outFile, compress));

        log.info("Conversion took " + (System.currentTimeMillis() - tStart) / 1e3 + " seconds.");
    }

    public void doLiftOver(final MapBlocks mapBlocks, final Collection<String> longislnd, final PrintStream ps) throws IOException {
        final Collection<ReadMapRecord> readMapRecords = new LongISLNDReadAlignmentMap(longislnd).getReadAlignmentMap().values();
        int readCount = 0;
        for (final ReadMapRecord readMapRecord : readMapRecords) {
            ps.println(mapBlocks.liftOverReadMapRecord(readMapRecord, minIntervalLength));
            if ((readCount % 100000) == 0) {
                log.info(readCount + " reads processed");
            }
            readCount++;
        }
        ps.close();
        log.info(readCount + " reads processed");
    }
}
