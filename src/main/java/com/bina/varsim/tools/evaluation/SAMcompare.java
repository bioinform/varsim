package com.bina.varsim.tools.evaluation;

/**
 * Read in SAM file with appropriately formatted read names and output accuracy statistics
 *
 * @author johnmu
 */


import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MapBlock.BlockType;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;
import com.bina.varsim.types.BedFile;
import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.stats.MapRatioRecordSum;
import com.bina.varsim.types.stats.StatsNamespace;
import com.bina.varsim.util.ReadMap;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.ObjectMapper;
import htsjdk.samtools.*;
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class SAMcompare {
    static final int WIGGLE_ARG = 20;
    static final int MAPQ_CUTOFF = 10;
    static final int MAPQ_UNMAPPED = 255;
    private final static Logger log = Logger.getLogger(SAMcompare.class.getName());
    @Option(name = "-wig", usage = "Wiggle allowance in validation [" + WIGGLE_ARG + "]")
    int wiggle = WIGGLE_ARG;
    @Option(name = "-mapq_cutoff", usage = "Mapping quality cutoff [" + MAPQ_CUTOFF + "]")
    int mapqCutoff = MAPQ_CUTOFF;
    @Option(name = "-prefix", usage = "Prefix for output file [Required]", metaVar = "file", required = true)
    String out_prefix;
    @Option(name = "-bed", usage = "BED file to restrict the analysis [Optional]", metaVar = "BED_file")
    String bed_filename = "";
    @Option(name = "-html", usage = "Insert JSON to HTML file [Optional, internal]", metaVar = "HTML_file", hidden = true)
    File html_file = null;
    @Option(name = "-read_map", usage = "Read MAP file [Optional]", metaVar = "ReadMap_file")
    File readMapFile = null;
    @Argument(usage = "One or more BAM files", metaVar = "bam_files ...", required = true)
    private ArrayList<String> bam_filename = new ArrayList<>();

    /**
     * @param args
     */
    public static void main(String[] args) {
        SAMcompare runner = new SAMcompare();
        runner.run(args);
    }

    /**
     * @param is_first is it the first read in a pair
     * @return return 0 for first, 1 for not first
     */
    private int getPairIdx(boolean is_first) {
        if (is_first) {
            return 0;
        } else {
            return 1;
        }
    }


    public void run(String[] args) {
        String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
        String usage = "Analyses the accuracy of the alignments in a SAM/BAM file\n" +
                "bed_file restricts the analysis to the bed regions\n";

        System.err.println(VERSION);

        CmdLineParser parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("java -jar VarSim.jar samcompare [options...] bam_files ...");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }

        final BedFile intersector = (bed_filename != null && new File(bed_filename).isFile()) ? new BedFile(bed_filename) : null;

        output_class output_blob = new output_class();

        // TODO think about better way to deal with multiple bams
        output_blob.setParams(new CompareParams(bam_filename.get(0).substring(0, Math.min(64, bam_filename.get(0).length())), wiggle, bed_filename));
        output_blob.setStats(new MapRatioRecordSum());

        // generate the output files
        PrintWriter jsonWriter = null;
        final Map<BlockType, SAMFileWriter> fpWriters = new HashMap<>();
        try {
            jsonWriter = new PrintWriter(out_prefix + "_report.json", "UTF-8");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        // read sam/bam file
        final SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();

        final SamReaderFactory factory =
                SamReaderFactory.makeDefault()
                        .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
                                SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
                        .validationStringency(ValidationStringency.LENIENT);

        final SAMFileWriter tumFPWriter = writerFactory.makeBAMWriter(factory.getFileHeader(new File(bam_filename.get(0))), false, new File(out_prefix + "_TUM_FP.bam"));
        final SAMFileWriter fpWriter = writerFactory.makeBAMWriter(factory.getFileHeader(new File(bam_filename.get(0))), false, new File(out_prefix + "_FP.bam"));
        for (final BlockType blockType : BlockType.values()) {
            if (blockType != BlockType.UNKNOWN) {
                fpWriters.put(blockType, writerFactory.makeBAMWriter(factory.getFileHeader(new File(bam_filename.get(0))), false, new File(out_prefix + "_" + blockType.getShortName() + "_FP.bam")));
            }
        }

        ReadMap readMap = null;
        try {
            readMap = readMapFile != null ? new ReadMap(readMapFile) : null;
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }

        int numReads = 0;
        for (String filename : bam_filename) {
            log.info("Reading file: " + filename);
            final SamReader reader = factory.open(new File(filename));
            try {
                for (SAMRecord rec : reader) {

                    // TODO this will need to change when we start considering supplementary alignements
                    if (rec.getNotPrimaryAlignmentFlag() || rec.getSupplementaryAlignmentFlag()) {
                        continue;
                    }

                    numReads++;

                    if (numReads % 100000 == 0) {
                        log.info("Read " + numReads + " records ...");
                    }

                    String name = rec.getReadName();

                    // parse the name
                    // TODO need to check for errors here
                    int pair_idx = rec.getReadPairedFlag() ? getPairIdx(rec.getFirstOfPairFlag()): 0;
                    log.trace("Getting true locations for read " + name);
                    final Collection<GenomeLocation> true_locs = readMap != null ? readMap.getReadMapRecord(name).getUnclippedStarts(pair_idx) : new SimulatedRead(name).getLocs(pair_idx);

                    if (!(intersector == null)) {
                        boolean contained_in_bed = false;
                        for (GenomeLocation loc : true_locs) {
                            if (intersector.contains(new ChrString(loc.chromosome), loc.location
                                    , loc.location + rec.getReadLength() - 1)) {
                                contained_in_bed = true;
                            }
                        }
                        if (!contained_in_bed) {
                            // skip this read
                            continue;
                        }
                    }

                    boolean true_unmapped;

                    Set<String> features = new HashSet<>(4);
                    features.add("All");

                    // determine if the read really should be unmapped
                    true_unmapped = true;
                    if (!true_locs.isEmpty()) {
                        // get types of genome features
                        // if only insertions and deletions then it is also unmapped
                        for (GenomeLocation loc : true_locs) {
                            final BlockType feat = loc.feature;
                            features.add(feat.getLongName());

                            // TODO check this with marghoob
                            true_unmapped &= !feat.isMappable();
                        }
                    }

                    if (true_unmapped) {
                        features.add("True_Unmapped");
                    } else {
                        output_blob.getStats().incStat(features, -1, StatsNamespace.T); // records the mappable reads
                    }


                    boolean unmapped = rec.getReadUnmappedFlag();
                    int mapping_quality = unmapped ? MAPQ_UNMAPPED : rec.getMappingQuality();
                    final StatsNamespace validationStatus;
                    if (unmapped) {
                        validationStatus = true_unmapped ? StatsNamespace.TN : StatsNamespace.FN;
                    } else {
                        // check if the it mapped to the correct location
                        boolean closeAln = false;
                        if (true_unmapped) {
                            if (mapping_quality > mapqCutoff) {
                                tumFPWriter.addAlignment(rec);
                            }
                        } else {
                            // Use unclipped location since the true locations are also unclipped
                            final GenomeLocation mappedLocation = new GenomeLocation(rec.getReferenceName(), rec.getUnclippedStart());

                            for (GenomeLocation loc : true_locs) {
                                closeAln |= loc.feature.isMappable() && loc.isClose(mappedLocation, wiggle);
                            }

                            if (!closeAln) {
                                if (mapping_quality > mapqCutoff) {
                                    fpWriter.addAlignment(rec);
                                    for (final BlockType blockType : BlockType.values()) {
                                        if (fpWriters.containsKey(blockType) && features.contains(blockType.getLongName())) {
                                            fpWriters.get(blockType).addAlignment(rec);
                                        }
                                    }
                                }
                            }
                        }
                        validationStatus = closeAln ? StatsNamespace.TP : StatsNamespace.FP;
                    }
                    output_blob.getStats().incStat(features, mapping_quality, validationStatus);
                }

            } finally {
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }

        log.info("Number of reads read: " + numReads);

        // output the statistics
        System.err.println("Statistics for all reads");
        System.err.println(output_blob.getStats());

        // output a JSON object
        ObjectMapper mapper = new ObjectMapper();
        mapper.configure(JsonGenerator.Feature.AUTO_CLOSE_TARGET, false);

        String jsonStr = "";
        try {
            jsonStr = mapper.writeValueAsString(output_blob);
            jsonWriter.print(jsonStr);
        } catch (Exception e) {
            e.printStackTrace();
        }
        jsonWriter.close();

        if (html_file != null) {
            try {
                FileUtils.writeStringToFile(new File(out_prefix + "_aligncomp.html"), JSONInserter.insertJSON(FileUtils.readFileToString(html_file), jsonStr));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        for (final SAMFileWriter fw : fpWriters.values()) {
            fw.close();
        }
        fpWriter.close();
        tumFPWriter.close();

        log.info("Done!"); // used to record the time
    }

    /**
     * This class is for outputting as a JSON
     */
    private static class output_class {
        CompareParams params;
        MapRatioRecordSum stats;

        output_class(CompareParams params, MapRatioRecordSum stats) {
            this.params = params;
            this.stats = stats;
        }

        output_class() {
        }

        public MapRatioRecordSum getStats() {
            return stats;
        }

        public void setStats(MapRatioRecordSum stats) {
            this.stats = stats;
        }

        public CompareParams getParams() {
            return params;
        }

        public void setParams(CompareParams params) {
            this.params = params;
        }
    }

    /**
     * Stores the parameters. this is mainly for outputting as a JSON.
     */
    class CompareParams {
        String bam_filename;
        int wiggle;
        String bed_filename;

        CompareParams(String bam_filename, int wiggle, String bed_filename) {
            this.bam_filename = bam_filename;
            this.wiggle = wiggle;
            this.bed_filename = bed_filename;
        }

        public String getBam_filename() {
            return bam_filename;
        }

        public void setBam_filename(String bam_filename) {
            this.bam_filename = bam_filename;
        }

        public int getWiggle() {
            return wiggle;
        }

        public void setWiggle(int wiggle) {
            this.wiggle = wiggle;
        }

        public String getBed_filename() {
            return bed_filename;
        }

        public void setBed_filename(String bed_filename) {
            this.bed_filename = bed_filename;
        }
    }

}
