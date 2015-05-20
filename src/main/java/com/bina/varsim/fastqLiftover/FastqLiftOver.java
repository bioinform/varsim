package com.bina.varsim.fastqLiftover;

import com.bina.varsim.fastqLiftover.readers.ARTPairedFastqAlnReader;
import com.bina.varsim.fastqLiftover.readers.DWGSIMPairedFastqReader;
import com.bina.varsim.fastqLiftover.readers.PBSIMFastqReader;
import com.bina.varsim.fastqLiftover.readers.PairedFastqReader;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MapBlocks;
import com.bina.varsim.fastqLiftover.types.SimulatedReadPair;
import org.apache.commons.codec.binary.Base64;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class FastqLiftOver {
    String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();

    private final static Logger log = Logger.getLogger(FastqLiftOver.class.getName());
    @Option(name = "-map", usage = "Map file", metaVar = "file")
    private File mapFile;
    @Option(name = "-type", usage = "Type of FASTQ (art/dwgsim)", metaVar = "fastqType")
    private String fastqType = "dwgsim";
    @Option(name = "-fastq", usage = "FASTQ input", metaVar = "file", required = true)
    private List<File> fastqFiles;
    @Option(name = "-aln", usage = "ART aln file", metaVar = "file")
    private List<File> alnFiles;
    @Option(name = "-maf", usage = "PBSIM maf file", metaVar = "file")
    private List<File> mafFiles;
    @Option(name = "-ref", usage = "PBSIM fasta headers of references", metaVar = "file")
    private List<File> refFiles;
    @Option(name = "-out", usage = "Output file", metaVar = "file")
    private List<File> outFiles;
    @Option(name = "-compress", usage = "Use GZIP compression")
    private boolean compress;
    @Option(name = "-id", usage = "Id to append to each read", required = true)
    private int laneId;
    @Option(name = "-force_five_base_encoding", usage = "For ART, force bases to be ACTGN")
    private boolean forceFiveBaseEncoding = false;

    public static void main(String[] args) throws IOException {
        new FastqLiftOver().run(args);
    }

    public InputStream decompressStream(final File inputFile) throws IOException {
        PushbackInputStream pb = new PushbackInputStream(new FileInputStream(inputFile), 1024 * 1024); //we need a pushbackstream to look ahead
        byte[] signature = new byte[2];
        pb.read(signature); //read the signature
        pb.unread(signature); //push back the signature to the stream
        if (signature[0] == (byte) (GZIPInputStream.GZIP_MAGIC & 0xff) && signature[1] == (byte) ((GZIPInputStream.GZIP_MAGIC >> 8) & 0xff)) {
            log.info(inputFile.getName() + " is gzip compressed");
            return new GZIPInputStream(pb, 1024 * 1024);
        } else {
            log.info(inputFile.getName() + " is not compressed");
            return pb;
        }
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
            System.err.println("java Fastq_liftOver [options...] arguments...");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println();
            return;
        }

        long tStart = System.currentTimeMillis();
        MapBlocks mapBlocks = new MapBlocks(mapFile);

        if (fastqType.equals("dwgsim")) {
            log.info("fastqFiles " + fastqFiles.get(0).getName() + " " + fastqFiles.get(1).getName());
            doLiftOverDwgsimFastqMap(mapBlocks, getOutStream(outFiles.get(0), compress), getOutStream(outFiles.get(1), compress));
        }

        if (fastqType.equals("art")) {
            log.info("fastqFiles " + fastqFiles.get(0).getName() + " " + fastqFiles.get(1).getName());
            log.info("alnFiles   " + alnFiles.get(0).getName() + " " + alnFiles.get(1).getName());
            log.info("outFiles   " + outFiles.get(0).getName() + " " + outFiles.get(1).getName());
            doLiftOverArtFastqMap(mapBlocks, getOutStream(outFiles.get(0), compress), getOutStream(outFiles.get(1), compress));
        }

        if (fastqType.equals("pbsim")) {
            log.info("fastqFiles " + fastqFiles.get(0).getName() + "  ");
            log.info("mafFiles   " + mafFiles.get(0).getName() + "  ");
            log.info("refFiles   " + refFiles.get(0).getName() + "  ");
            log.info("outFiles   " + outFiles.get(0).getName() + "  ");
            doLiftOverPbsimFastqMap(mapBlocks, getOutStream(outFiles.get(0), compress));
        }

        log.info("Conversion took " + (System.currentTimeMillis() - tStart) / 1e3 + " seconds.");
    }

    public void doLiftOverArtFastqMap(final MapBlocks mapBlocks, final PrintStream ps1, final PrintStream ps2) throws IOException {
        BufferedReader brFastq1 = new BufferedReader(new InputStreamReader(decompressStream(fastqFiles.get(0))));
        BufferedReader brAln1 = new BufferedReader(new InputStreamReader(decompressStream(alnFiles.get(0))));
        BufferedReader brFastq2 = new BufferedReader(new InputStreamReader(decompressStream(fastqFiles.get(1))));
        BufferedReader brAln2 = new BufferedReader(new InputStreamReader(decompressStream(alnFiles.get(1))));

        doLiftOverPairedFastq(mapBlocks, new ARTPairedFastqAlnReader(brAln1, brFastq1, brAln2, brFastq2, forceFiveBaseEncoding), ps1, ps2);
    }

    public void doLiftOverDwgsimFastqMap(final MapBlocks mapBlocks, final PrintStream ps1, final PrintStream ps2) throws IOException {
        BufferedReader br1 = new BufferedReader(new InputStreamReader(decompressStream(fastqFiles.get(0))));
        BufferedReader br2 = new BufferedReader(new InputStreamReader(decompressStream(fastqFiles.get(1))));

        doLiftOverPairedFastq(mapBlocks, new DWGSIMPairedFastqReader(br1, br2, forceFiveBaseEncoding), ps1, ps2);
    }

    public void doLiftOverPbsimFastqMap(final MapBlocks mapBlocks, final PrintStream ps) throws IOException {
        BufferedReader brFastq = new BufferedReader(new InputStreamReader(decompressStream(fastqFiles.get(0))));
        BufferedReader brMAF = new BufferedReader(new InputStreamReader(decompressStream(mafFiles.get(0))));
        BufferedReader brREF = new BufferedReader(new InputStreamReader(decompressStream(refFiles.get(0))));

        doLiftOverPairedFastq(mapBlocks, new PBSIMFastqReader(brREF, brMAF, brFastq, forceFiveBaseEncoding), ps, null);
    }

    public void doLiftOverPairedFastq(final MapBlocks mapBlocks, final PairedFastqReader pfr, final PrintStream ps1, final PrintStream ps2) throws IOException {
        SimulatedReadPair readPair;
        int readCount = 0;
        while ((readPair = pfr.getNextReadPair()) != null) {

            final boolean hasRead2 = readPair.read2 != null;
            if (hasRead2 && (ps2 == null)) throw new RuntimeException("found read2 but read2 output not provided");

            readPair.read1.laneId = laneId;
            if (hasRead2) {
                readPair.read2.laneId = laneId;
            }

            // if read2 not present, set locs2 to output the same as locs1 to preserve pair-ended behavior
            final GenomeLocation loc1 = readPair.read1.locs1.get(0);
            final GenomeLocation loc2 = hasRead2 ? readPair.read2.locs2.get(0) : loc1;
            final List<GenomeLocation> newLocs1 = mapBlocks.liftOverInterval(loc1.chromosome, loc1.location, loc1.location + readPair.read1.alignedBases1, loc1.direction);
            final List<GenomeLocation> newLocs2 = hasRead2 ? mapBlocks.liftOverInterval(loc2.chromosome, loc2.location, loc2.location + readPair.read2.alignedBases2, loc2.direction) : newLocs1;

            String readId = Base64.encodeBase64String(String.valueOf(readCount).getBytes());

            readPair.read1.locs1 = newLocs1;
            readPair.read1.locs2 = newLocs2;
            readPair.read1.setReadId(readPair.read1.getReadId() + readId);
            readPair.read1.origLocs1 = new ArrayList<>();
            readPair.read1.origLocs2 = new ArrayList<>();
            ++readCount;
            ps1.println(readPair.read1);

            if (hasRead2) {
                readPair.read2.locs1 = newLocs1;
                readPair.read2.locs2 = newLocs2;
                readPair.read2.setReadId(readPair.read2.getReadId() + readId);
                readPair.read2.origLocs1 = new ArrayList<>();
                readPair.read2.origLocs2 = new ArrayList<>();
                ++readCount;
                ps2.println(readPair.read2);
            }

            if ((readCount % 1000000) == 0) {
                log.info(readCount + " reads processed");
            }
        }
        ps1.close();
        if (ps2 != null) ps2.close();
        log.info("Processed " + readCount + " reads.");
    }
}
