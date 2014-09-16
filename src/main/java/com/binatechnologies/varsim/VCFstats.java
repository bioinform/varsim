package com.binatechnologies.varsim;

/**
 *
 */
import org.apache.log4j.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import java.util.HashMap;
import java.util.Map;

/**
 * @author johnmu
 */

class Stats_record {
    int[] bin_counts; // the last bin is for anything larger
    int total_count;
    int[] bin_breaks = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 29, 39, 49, 99,
            199, 399, 799, 1599, 3199, 6399, 12799, 25599, 51199, 102399, 500000, 1000000};
    int num_bins;

    Stats_record() {
        num_bins = bin_breaks.length + 1;
        bin_counts = new int[num_bins];
        total_count = 0;
    }

    /**
     * Increments the bin that val belongs in
     *
     * @param val value to be added
     */
    public void add(int val) {
        total_count++;
        for (int i = 0; i < bin_breaks.length; i++) {
            if (val <= bin_breaks[i]) {
                bin_counts[i]++;
                return;
            }
        }
        bin_counts[num_bins - 1]++;
    }

    public int getTotal_count() {
        return total_count;
    }

    public String toString() {
        return toString(bin_breaks[bin_breaks.length - 1] + 1);
    }

    public String toString(int max_len) {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        sb.append(1);
        sb.append(",");
        sb.append(bin_breaks[0]);
        sb.append("]");
        sb.append(':');
        sb.append(bin_counts[0]);
        sb.append('\n');
        for (int i = 1; i < bin_breaks.length; i++) {

            if (bin_breaks[i] > max_len) {
                break;
            }

            sb.append("[");
            sb.append(bin_breaks[i - 1] + 1);
            sb.append(",");
            sb.append(bin_breaks[i]);
            sb.append("]");
            sb.append(':');
            sb.append(bin_counts[i]);
            sb.append('\n');
        }
        if (bin_breaks[bin_breaks.length - 1] < max_len) {
            sb.append("[");
            sb.append(bin_breaks[bin_breaks.length - 1] + 1);
            sb.append(",");
            sb.append("inf");
            sb.append("]");
            sb.append(':');
            sb.append(bin_counts[num_bins - 1]);
            sb.append('\n');
        }
        return sb.toString();
    }
}

class Type_record {
    HashMap<Variant.Type, Stats_record> data;

    Type_record() {
        data = new HashMap<Variant.Type, Stats_record>();
    }

    /**
     * Increments the bin (of type) that val belongs in
     *
     * @param type type of bin
     * @param val  value to be incremented
     */
    public void add(Variant.Type type, int val) {
        if (data.containsKey(type)) {
            Stats_record rec = data.get(type);
            rec.add(val);
        } else {
            Stats_record rec = new Stats_record();
            rec.add(val);
            data.put(type, rec);
        }
    }

    public int getTotal_nonref() {
        int total = 0;
        for (Map.Entry<Variant.Type, Stats_record> entry : data.entrySet()) {
            if (entry.getKey() == Variant.Type.Reference)
                continue;
            total += entry.getValue().getTotal_count();
        }
        return total;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<Variant.Type, Stats_record> entry : data.entrySet()) {
            if (entry.getKey() == Variant.Type.Reference)
                continue;

            if (entry.getKey() == Variant.Type.SNP) {
                sb.append(entry.getKey().name());
                sb.append('\n');
                sb.append(entry.getValue().toString(1));
            } else {
                sb.append(entry.getKey().name());
                sb.append('\n');
                sb.append(entry.getValue());
            }
        }
        return sb.toString();
    }
}

/**
 * Stores the data for each haplotype separately
 */
class Parent_record {
    // 0 = paternal, 1 = maternal
    public final static int PATERNAL = 0;
    public final static int MATERNAL = 1;
    Type_record[] data;
    int total_count;

    Parent_record() {
        data = new Type_record[2];
        data[PATERNAL] = new Type_record();
        data[MATERNAL] = new Type_record();
        total_count = 0;
    }

    public void add(Variant var, BedFile bed_file) {
        int paternal_allele = var.paternal();
        int maternal_allele = var.maternal();

        boolean added = false;
        if (bed_file == null
                || bed_file.contains(var.getChr_name(), var.get_interval(paternal_allele))) {
            data[PATERNAL].add(var.getType(paternal_allele), var.max_len(paternal_allele));
            added = true;
        }

        if (bed_file == null
                || bed_file.contains(var.getChr_name(), var.get_interval(maternal_allele))) {
            data[MATERNAL].add(var.getType(maternal_allele), var.max_len(maternal_allele));
            added = true;
        }
        if (added) {
            total_count++;
        }
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Total count: ");
        sb.append(total_count);
        sb.append("\n");
        sb.append("Paternal\n");
        sb.append("Total: ");
        sb.append(data[PATERNAL].getTotal_nonref());
        sb.append("\n");
        sb.append(data[PATERNAL]);
        sb.append("\n");
        sb.append("Maternal\n");
        sb.append("Total: ");
        sb.append(data[MATERNAL].getTotal_nonref());
        sb.append("\n");
        sb.append(data[MATERNAL]);
        sb.append("\n");
        return sb.toString();
    }
}

public class VCFstats {
    private final static Logger log = Logger.getLogger(VCF2diploid.class.getName());

    BedFile bed_file;

    @Option(name = "-bed", usage = "BED file to restrict the analysis [Optional]",metaVar = "BED_file")
    String bed_filename = null;

    @Option(name = "-vcf", usage = "VCF file to analyse [Required]",metaVar = "VCF_file",required = true)
    String vcf_filename;

    /**
     * @param args command line parameters
     */
    public static void main(String[] args) {
        VCFstats runner = new VCFstats();
        runner.run(args);
    }

    private void run(String[] args) {
        String usage = "Compute the distribution of variants in a VCF file\n"
                + "Assume that there is only one sample in the VCF file\n";

        CmdLineParser cmd_parser = new CmdLineParser(this);

        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        cmd_parser.setUsageWidth(80);

        try {
            cmd_parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("java -jar vcfstats.jar [options...]");
            // print the list of available options
            cmd_parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }


        bed_file = null;

        if (bed_filename != null) {
            log.info("Reading BED file...");
            bed_file = new BedFile(bed_filename);
        }

        Parent_record data = new Parent_record();

        VCFparser parser = new VCFparser(vcf_filename, null, false);

        while (parser.hasMoreInput()) {
            Variant var = parser.parseLine();
            if (var == null) {
                continue;
            }
            if (var.maternal() == 0 && var.paternal() == 0) {
                continue;
            }
            data.add(var, bed_file);
        }

        System.out.println(data);
    }

}
