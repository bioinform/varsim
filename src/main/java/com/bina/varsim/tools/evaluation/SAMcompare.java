package com.bina.varsim.tools.evaluation;

/**
 * Read in SAM file with appropriately formatted read names and output accuracy statistics
 *
 * truth alignment can also be supplied via a read map file
 *
 * @author johnmu
*/


import com.bina.varsim.VarSimTool;
import com.bina.varsim.VarSimToolNamespace;
import com.bina.varsim.fastqLiftover.types.GenomeLocation;
import com.bina.varsim.fastqLiftover.types.MapBlock.BlockType;
import com.bina.varsim.fastqLiftover.types.SimulatedRead;
import com.bina.varsim.types.BedFile;
import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.stats.MapRatioRecordSum;
import com.bina.varsim.types.stats.MapRatioRecordSum.EventTypesForStats;
import com.bina.varsim.types.stats.StatsNamespace;
import com.bina.varsim.util.ReadMap;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import htsjdk.samtools.*;
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class SAMcompare extends VarSimTool {
    static final int WIGGLE_ARG = 20;
    static final int MAPQ_CUTOFF = 10;
    static final int MAPQ_UNMAPPED = 255;
    static final double DEFAULT_IDENTITY_THRESHOLD = 1.0;
    //bwa outputs matches as M, whereas = is defined as match by SAM spec
    static final Set<CigarOperator> ALL_POSSIBLE_MATCHES = EnumSet.of(CigarOperator.EQ, CigarOperator.MATCH_OR_MISMATCH);
    private final static Logger log = Logger.getLogger(SAMcompare.class.getName());
    @Option(name = "-wig", usage = "Wiggle allowance in validation [" + WIGGLE_ARG + "]")
    int wiggle = WIGGLE_ARG;
    @Option(name = "-mapq_cutoff", usage = "Mapping quality cutoff [" + MAPQ_CUTOFF + "]")
    int mapqCutoff = MAPQ_CUTOFF;
    @Option(name = "-prefix", usage = "Prefix for output file [Required]", metaVar = "file", required = true)
    String outPrefix;
    @Option(name = "-bed", usage = "BED file to restrict the analysis [Optional]", metaVar = "BED_file")
    String bedFilename = "";
    @Option(name = "-html", usage = "Insert JSON to HTML file [Optional, internal]", metaVar = "HTML_file", hidden = true)
    File htmlFile = null;
    @Option(name = "-read_map", usage = "Read MAP file [Optional]", metaVar = "ReadMap_file")
    File readMapFile = null;
    @Option(name = "-use_nonprimary", usage = "Do not skip non-primary alignments")
    boolean useNonPrimary = false;
    @Option(name = "-identity_threshold", usage = "an alignment with higher identity (% of matches in a read) will be considered correct alignment regardless of mapping location. DEFAULT: 1.0")
    double identityThreshold = DEFAULT_IDENTITY_THRESHOLD;
    @Argument(usage = "One or more BAM files, header coming from the first file", metaVar = "bam_files ...", required = true)
    private List<String> bamFilenames = new ArrayList<>();

    /**
     * wrapper for actual execution code
     * @param command placeholder
     * @param description program information
     */
    public SAMcompare(final String command, final String description) {
        super(command, description);
    }

    /**
     * calls the execution code with command line arguments
     * @param args command line options
     */
    public static void main(String[] args) {
        new SAMcompare("", VarSimToolNamespace.SAMCompare.description).run(args);
    }

    /**
     * returns index of read in a pair
     * @param isFirst is it the first read in a pair
     * @return return 0 for first, 1 for not first
     */
    private int getPairIdx(boolean isFirst) {
        return isFirst? 0:1;
    }


    /**
     * main execution code block
     * expects alignments in SAM/BAM format
     * read in alignment line by line, filter by
     * alignment quality and mapping location
     * determine if an alignment is false positive
     * or not.
     * outputs false positive alignments by type of sequence changes in
     * the perturbed genome and a statistical report in json format
     *
     * @param args command line arguments
     */
    public void run(String[] args) {
        if (!parseArguments(args)) {
            return;
        }

        OutputClass outputBlob = new OutputClass();

        // TODO think about better way to deal with multiple bams
        outputBlob.setParams(new CompareParams(bamFilenames.get(0).substring(0, Math.min(64, bamFilenames.get(0).length())), wiggle, bedFilename, identityThreshold));
        outputBlob.setStats(new MapRatioRecordSum());

        // read sam/bam file
        final SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
        final SamReaderFactory factory =
                SamReaderFactory.makeDefault()
                        .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
                                SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
                        .validationStringency(ValidationStringency.LENIENT);
        // generate the output files
        final Map<BlockType, SAMFileWriter> fpWriters = new EnumMap<>(BlockType.class);
        //TODO: put file open into try()
        try (
                final PrintWriter jsonWriter = new PrintWriter(outPrefix + "_report.json", "UTF-8");
                final SAMFileWriter tumFPWriter = writerFactory.makeBAMWriter(
                        factory.getFileHeader(new File(bamFilenames.get(0))), false, new File(outPrefix + "_TUM_FP.bam"));
                final SAMFileWriter fpWriter = writerFactory.makeBAMWriter(
                        factory.getFileHeader(new File(bamFilenames.get(0))), false, new File(outPrefix + "_FP.bam"));
        ) {
            for (final BlockType blockType : BlockType.values()) {
                if (blockType != BlockType.UNKNOWN) {
                    fpWriters.put(blockType, writerFactory.makeBAMWriter(factory.getFileHeader(new File(bamFilenames.get(0))), false, new File(outPrefix + "_" + blockType.getShortName() + "_FP.bam")));
                }
            }

            ReadMap readMap = readMapFile != null ? new ReadMap(readMapFile) : null;
            final BedFile intersector = (bedFilename != null && new File(bedFilename).isFile()) ? new BedFile(bedFilename) : null;

            int numReads = 0;
            for (String filename : bamFilenames) {
                log.info("Reading file: " + filename);
                final SamReader reader = factory.open(new File(filename));
                for (SAMRecord rec : reader) {

                    // TODO this will need to change when we start considering supplementary alignements
                    if (!useNonPrimary && (rec.getNotPrimaryAlignmentFlag() || rec.getSupplementaryAlignmentFlag())) {
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
                    final Collection<GenomeLocation> trueLoci = readMap != null ? readMap.getReadMapRecord(name).getUnclippedStarts(pair_idx) : new SimulatedRead(name).getLocs(pair_idx);

                    if (!isContainedInBed(intersector, trueLoci, rec)) {
                        // skip this read
                        continue;
                    }
                    boolean isTrueUnmapped;

                    //TODO: create enum and replace HashSet with EnumSet
                    Set<EventTypesForStats> features = EnumSet.noneOf(EventTypesForStats.class);
                    features.add(EventTypesForStats.All);

                    // determine if the read really should be unmapped
                    isTrueUnmapped = true;
                    if (!trueLoci.isEmpty()) {
                        // get types of genome features
                        // if only insertions and deletions then it is also unmapped
                        for (GenomeLocation loc : trueLoci) {
                            final BlockType feat = loc.feature;
                            features.add(EventTypesForStats.valueOf(feat.getLongName()));

                            // TODO check this with marghoob
                            isTrueUnmapped &= !feat.isMappable();
                        }
                    }

                    if (isTrueUnmapped) {
                        features.add(EventTypesForStats.True_Unmapped);
                    } else {
                        outputBlob.getStats().incStat(features, -1, StatsNamespace.T); // records the mappable reads
                    }


                    boolean unmapped = rec.getReadUnmappedFlag();
                    //TODO: should mapq=0 for unmapped reads? otherwise all unmapped=true && true_unmapped=true reads will go to FP
                    int mappingQuality = unmapped ? MAPQ_UNMAPPED : rec.getMappingQuality();
                    StatsNamespace validationStatus = StatsNamespace.TP;
                    if (unmapped) {
                        validationStatus = isTrueUnmapped ? StatsNamespace.TN : StatsNamespace.FN;
                    } else {
                        /*
                         * some explanations about how an alignment is classified as TP or FP
                          *
                          * 1) if a read is aligned to correct location, then it is for sure TP
                          * 2) if a read is not aligned to correct location, but has very high identity (true multi-alignment), then TP
                          * 3) if a read is not aligned to correct location, doesn't have very high identity, but is assigned very informative mapq (low), then TP
                          * 4) if a read is supposed to be unmapped, but is mapped given low identity and is assigned high mapq, then FP
                          * 5) if a read is aligned, not FP, then it is TP
                          *
                          *
                         */
                        // check if the it mapped to the correct location
                        boolean closeAln = false;
                        if (isTrueUnmapped) {
                            if (getIdentity(rec) <= identityThreshold && mappingQuality > mapqCutoff) {
                                tumFPWriter.addAlignment(rec);
                            }
                        } else {
                            // Use unclipped location since the true locations are also unclipped
                            final GenomeLocation mappedLocation = new GenomeLocation(new ChrString(rec.getReferenceName()), rec.getUnclippedStart());

                            //if a read can be mapped to one of possible loci
                            //then we consider it as mapped
                            for (GenomeLocation loc : trueLoci) {
                                closeAln |= loc.feature.isMappable() && loc.isClose(mappedLocation, wiggle);
                            }

                            if (!closeAln) {
                                if (getIdentity(rec) <= identityThreshold && mappingQuality > mapqCutoff) {
                                    //got another false positive
                                    fpWriter.addAlignment(rec);
                                    validationStatus = StatsNamespace.FP;
                                    for (final BlockType blockType : BlockType.values()) {
                                        if (fpWriters.containsKey(blockType) && features.contains(EventTypesForStats.valueOf(blockType.getLongName()))) {
                                            //categorize false positive alignments
                                            fpWriters.get(blockType).addAlignment(rec);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    outputBlob.getStats().incStat(features, mappingQuality, validationStatus);
                }
                reader.close();
            }

            log.info("Number of reads read: " + numReads);

            // output the statistics
            log.info("Statistics for all reads");
            log.info(outputBlob.getStats());

            // output a JSON object
            ObjectMapper mapper = new ObjectMapper();
            mapper.configure(JsonGenerator.Feature.AUTO_CLOSE_TARGET, false);

            String jsonStr = "";
            try {
                jsonStr = mapper.writeValueAsString(outputBlob);
            } catch (JsonProcessingException e) {
                e.printStackTrace();
                System.exit(1);
            }
            jsonWriter.print(jsonStr);
            if (htmlFile != null) {
                FileUtils.writeStringToFile(new File(outPrefix + "_aligncomp.html"), JSONInserter.insertJSON(FileUtils.readFileToString(htmlFile), jsonStr));
            }
            for (final BlockType blockType : BlockType.values()) {
                if (blockType != BlockType.UNKNOWN) {
                    fpWriters.get(blockType).close();
                }
            }
            fpWriters.values().stream().forEach(SAMFileWriter::close);
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
        log.info("Done!"); // used to record the time
    }

    /**
     * This class is for outputting as a JSON
     */
    private static class OutputClass {
        CompareParams params;
        MapRatioRecordSum stats;

        /**
         * store parameters and statistics
         * @param params
         * @param stats
         */
        OutputClass(CompareParams params, MapRatioRecordSum stats) {
            this.params = params;
            this.stats = stats;
        }

        /**
         * initialize empty output object
         */
        OutputClass() {
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
        @JsonProperty(value = "bam_filename")
        String bamFilename;
        @JsonProperty(value = "wiggle")
        int wiggle;
        @JsonProperty(value = "bed_filename")
        String bedFilename;
        @JsonProperty(value = "identity_threshold")
        double identityThreshold;

        CompareParams(String bamFilename, int wiggle, String bedFilename, double identityThreshold) {
            this.bamFilename = bamFilename;
            this.wiggle = wiggle;
            this.bedFilename = bedFilename;
            this.identityThreshold = identityThreshold;
        }

        public String getBamFilename() {
            return bamFilename;
        }

        public void setBamFilename(String bamFilename) {
            this.bamFilename = bamFilename;
        }

        public int getWiggle() {
            return wiggle;
        }

        public void setWiggle(int wiggle) {
            this.wiggle = wiggle;
        }

        public String getBedFilename() {
            return bedFilename;
        }

        public void setBedFilename(String bedFilename) {
            this.bedFilename = bedFilename;
        }
    }

    /**
     * check if a read is supposed to be aligned in the user-supplied BED region
     * @param intersector BED-file object
     * @param trueLoci collection of true alignment loci
     * @param rec one alignment in SAM
     * @return
     */
    private static boolean isContainedInBed(BedFile intersector, Collection<GenomeLocation> trueLoci, SAMRecord rec) {
        if (intersector == null) {
            return true;
        }
        boolean isContainedInBed = false;
        for (GenomeLocation loc : trueLoci) {
            if (intersector.contains(loc.chromosome, loc.location
                    //here we should have used the true alignment length
                    //however that is difficult to know unless we have a perfect BAM
                    , loc.location + rec.getReadLength() - 1)) {
                isContainedInBed = true;
            }
        }
        return isContainedInBed;
    }

    /**
     * calculate identity ratio of the read
     * formula:
     *
     * sum of sequence matches (= or M in CIGAR) / total length of the read
     *
     * @param alignment
     * @return
     */
    private static double getIdentity(SAMRecord alignment) {
        int nMatches = 0;
        for (CigarElement c : alignment.getCigar().getCigarElements()) {
            if (ALL_POSSIBLE_MATCHES.contains(c.getOperator())) {
                nMatches += c.getLength();
            }
        }
        return (double) nMatches/alignment.getReadLength();
    }
}
