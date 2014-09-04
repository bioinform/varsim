package com.binatechnologies.seqalto.varsim;

/**
 * Read in SAM file with appropriately formatted read names and output accuracy statistics
 * @author johnmu
 */


import com.binatechnologies.varsim.fastq_liftover.GenomeLocation;
import com.binatechnologies.varsim.fastq_liftover.MapBlock;
import com.binatechnologies.varsim.fastq_liftover.SimulatedRead;
import com.fasterxml.jackson.databind.ObjectMapper;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class SAMcompare {
    private final static Logger log = Logger.getLogger(VCFcompare.class.getName());
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

    /**
     * @param orig_name original variant feature name as a string
     * @return convert the string name to a variable
     */
    private BlockTypeOut convert_feature_name(String orig_name) {
        BlockTypeOut out = BlockTypeOut.UNKNOWN;

        if (orig_name.equals(MapBlock.BlockType.DEL.toString())) {
            out = BlockTypeOut.DEL;
        } else if (orig_name.equals(MapBlock.BlockType.INS.toString())) {
            out = BlockTypeOut.INS;
        } else if (orig_name.equals(MapBlock.BlockType.DUP_TANDEM.toString())) {
            out = BlockTypeOut.DUP_TANDEM;
        } else if (orig_name.equals(MapBlock.BlockType.INV.toString())) {
            out = BlockTypeOut.INV;
        } else if (orig_name.equals(MapBlock.BlockType.SEQ.toString())) {
            out = BlockTypeOut.SEQ;
        } else if (orig_name.equals(MapBlock.BlockType.UNKNOWN.toString())) {
            out = BlockTypeOut.UNKNOWN;
        }

        return out;
    }

    public void run(String[] args) {
        String usage = "SAMcompare wiggle output_prefix [bed_file] bam_files ...\n" +
                "Analyses the accuracy of the alignments in a SAM/BAM file\n" +
                "bed_file restricts the analysis to the bed regions\n";
        // TODO args4j!!!
        if (args.length < 4) {
            System.err.println(usage);
            System.exit(1);
        }

        int wiggle = Integer.parseInt(args[0]);
        String out_prefix = args[1];

        String bed_filename = "";
        BedFile intersector = null;

        bed_filename = args[2];
        boolean bed_exists = false;
        // check if the file exists
        try{
            File f = new File(bed_filename);
            if(f.exists()){
                bed_exists = true;
            }
        }catch (Exception e){

        }

        if(bed_exists) {
            intersector = new BedFile(bed_filename);
        }

        // iterate to get all the bam files
        ArrayList<String> bam_filename = new ArrayList<String>(1);
        for(int i = 3;i<args.length;i++){
            bam_filename.add(args[i]);
        }



        /**
         * This class is for outputting as a JSON
         */
        class output_class {
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


        output_class output_blob = new output_class();

        // TODO think about better way to deal with multiple bams
        output_blob.setParams(new CompareParams(bam_filename.get(0).substring(0,Math.min(64,bam_filename.get(0).length())), wiggle, bed_filename));
        output_blob.setStats(new MapRatioRecordSum());

        // generate the output files
        PrintWriter JSON_writer = null;
        PrintWriter FP_writer = null;
        PrintWriter DEL_FP_writer = null;
        PrintWriter INV_FP_writer = null;
        PrintWriter INS_FP_writer = null;
        PrintWriter SEQ_FP_writer = null;
        PrintWriter TD_FP_writer = null;
        PrintWriter TUM_FP_writer = null;
        try {
            JSON_writer = new PrintWriter(out_prefix + "_report.json", "UTF-8");
            FP_writer = new PrintWriter(out_prefix + "_FP.SAM", "UTF-8");
            DEL_FP_writer = new PrintWriter(out_prefix + "_DEL_FP.SAM", "UTF-8");
            INV_FP_writer = new PrintWriter(out_prefix + "_INV_FP.SAM", "UTF-8");
            INS_FP_writer = new PrintWriter(out_prefix + "_INS_FP.SAM", "UTF-8");
            SEQ_FP_writer = new PrintWriter(out_prefix + "_SEQ_FP.SAM", "UTF-8");
            TD_FP_writer = new PrintWriter(out_prefix + "_TD_FP.SAM", "UTF-8");
            TUM_FP_writer = new PrintWriter(out_prefix + "_TUM_FP.SAM", "UTF-8");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        // read sam/bam file
        int num_read = 0;

        for(String filename : bam_filename) {
            log.info("Reading file: " + filename);
            SAMFileReader reader = new SAMFileReader(new File(filename));
            try {
                reader.setValidationStringency(ValidationStringency.LENIENT);

                SAMRecordIterator it = reader.iterator();
                while (it.hasNext()) {
                    SAMRecord rec = it.next();
                    num_read++;

                    if (num_read % 100000 == 0) {
                        log.info("Read " + num_read + " records ...");
                    }

                    String name = rec.getReadName();

                    // parse the name
                    // TODO need to check for errors here
                    SimulatedRead true_read = new SimulatedRead(name);

                    int pair_idx = getPairIdx(rec.getFirstOfPairFlag());

                    List<GenomeLocation> true_locs;
                    if (pair_idx == 0) {
                        true_locs = true_read.locs1;
                    } else {
                        true_locs = true_read.locs2;
                    }


                    if (!(intersector == null)) {
                        boolean contained_in_bed = false;
                        for (GenomeLocation loc : true_locs) {
                            if (intersector.contains(loc.chromosome, loc.location
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


                    HashSet<String> features = new HashSet<String>(4);
                    features.add("All");

                    // determine if the read really should be unmapped
                    if (true_locs.isEmpty()) {
                        true_unmapped = true;
                    } else {
                        true_unmapped = true;

                        // get types of genome features
                        // if only insertions and deletions then it is also unmapped
                        for (GenomeLocation loc : true_locs) {
                            BlockTypeOut feat = convert_feature_name(loc.feature);
                            features.add(feat.toString());

                            // TODO check this with marghoob
                            if (!(feat == BlockTypeOut.INS) && !(feat == BlockTypeOut.DEL)) {
                                true_unmapped = false;
                            }
                        }
                    }

                    if (true_unmapped) {
                        features.add("True_Unmapped");
                    }

                    output_blob.getStats().incT(features);

                    boolean unmapped = rec.getReadUnmappedFlag();
                    int mapping_quality = rec.getMappingQuality();

                    if (unmapped) {
                        // set mapping quality of unmapped reads to 255 since it is unknown
                        // otherwise it will not be counted in the simulation
                        mapping_quality = 255;

                        if (true_unmapped) {
                            // correctly aligned
                            output_blob.getStats().incTN(features, mapping_quality);
                        } else {
                            // incorrectly unmapped
                            output_blob.getStats().incFN(features, mapping_quality);
                        }
                    } else {
                        // check if the it mapped to the correct location
                        if (true_unmapped) {
                            output_blob.getStats().incFP(features, mapping_quality);

                            if (mapping_quality > 10) {
                                TUM_FP_writer.println(rec.getSAMString());

                            }
                        } else {
                            String chr = rec.getReferenceName();

                            // unclipped start adjusts for clipped bases
                            int pos = rec.getUnclippedStart();

                            boolean good_aln = false;

                            for (GenomeLocation loc : true_locs) {
                                BlockTypeOut feat = convert_feature_name(loc.feature);

                                if ((feat == BlockTypeOut.INS) || (feat == BlockTypeOut.DEL)) {
                                    continue;
                                }

                                if (chr.equals(loc.chromosome)) {

                                    if (Math.abs(loc.location - pos) <= wiggle) {
                                        good_aln = true;
                                    }
                                }
                            }

                            if (good_aln) {
                                output_blob.getStats().incTP(features, mapping_quality);
                            } else {
                                output_blob.getStats().incFP(features, mapping_quality);
                                if (mapping_quality > 10) {
                                    FP_writer.println(rec.getSAMString());
                                    if (features.contains(BlockTypeOut.DEL.toString())) {
                                        DEL_FP_writer.println(rec.getSAMString());
                                    }

                                    if (features.contains(BlockTypeOut.INV.toString())) {
                                        INV_FP_writer.println(rec.getSAMString());
                                    }

                                    if (features.contains(BlockTypeOut.INS.toString())) {
                                        INS_FP_writer.println(rec.getSAMString());
                                    }

                                    if (features.contains(BlockTypeOut.SEQ.toString())) {
                                        SEQ_FP_writer.println(rec.getSAMString());
                                    }

                                    if (features.contains(BlockTypeOut.DUP_TANDEM.toString())) {
                                        TD_FP_writer.println(rec.getSAMString());
                                    }
                                }
                            }
                        }
                    }
                }

            } finally {
                reader.close();
            }

        }

        log.info("Number of reads read: " + num_read);

        // output the statistics
        System.err.println("Statistics for all reads");
        System.err.println(output_blob.getStats());

        // output a JSON object
        ObjectMapper mapper = new ObjectMapper();
        try {
            JSON_writer.print(mapper.writeValueAsString(output_blob));
        } catch (Exception e) {
            e.printStackTrace();
        }

        JSON_writer.close();
        FP_writer.close();
        DEL_FP_writer.close();
        INV_FP_writer.close();
        INS_FP_writer.close();
        SEQ_FP_writer.close();
        TD_FP_writer.close();
        TUM_FP_writer.close();
    }

    /**
     * This encodes the variant type as a string.
     * WARNING: the string must not have spaces, or it will break the javascript later
     */
    public enum BlockTypeOut {
        SEQ("Sequence"), INS("Insertion"), DEL("Deletion"), INV("Inversion"), DUP_TANDEM("Tandem_Duplication"), UNKNOWN("Unknown");

        private final String name;

        private BlockTypeOut(String s) {
            name = s;
        }

        public String toString() {
            return name;
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
