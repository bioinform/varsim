package com.binatechnologies.seqalto.varsim;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.log4j.Logger;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

/**
 * Compare two VCF files, output the TPR and FDR for various bins and variant types
 *
 * @author johnmu
 */

// TODO this class does not yet compare the contents of the variants
// TODO it also ignores genotypes for now
public class VCFcompare {
    private final static Logger log = Logger.getLogger(VCFcompare.class.getName());

    public static void main(String[] args) {
        VCFcompare runner = new VCFcompare();
        runner.run(args);
    }

    // end = true will add the indel to the end, other wise it will add to start
    private void add_indels(ArrayList<Variant> var_list, int[] diff, byte[] ref, byte[][] alt,
                            Variant var, int curr_pos, boolean end) {
        // add insertions or deletions for complex variants
        if (diff[0] == diff[1] && diff[0] != 0) {
            // homozygous
            if (diff[0] > 0) {
                // insertion
                if (Arrays.equals(alt[0], alt[1])) {
                    byte[] phase = {1, 1};

                    if(end){
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + ref.length, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[0], 0, diff[0]))},
                                phase, true, var.getVar_id(), ".", ""));
                    }else {
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[0], 0, diff[0]))},
                                phase, true, var.getVar_id(), ".", ""));
                    }
                } else {
                    byte[] phase = {0, 0};
                    if(end) {
                        phase[0] = 1;
                        phase[1] = 0;
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + ref.length, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[0], 0, diff[0]))},
                                phase, true, var.getVar_id(), ".", ""));
                        phase[0] = 0;
                        phase[1] = 1;
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + ref.length, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[1], 0, diff[1]))},
                                phase, true, var.getVar_id(), ".", ""));
                    }else{
                        phase[0] = 1;
                        phase[1] = 0;
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[0], 0, diff[0]))},
                                phase, true, var.getVar_id(), ".", ""));
                        phase[0] = 0;
                        phase[1] = 1;
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[1], 0, diff[1]))},
                                phase, true, var.getVar_id(), ".", ""));
                    }
                }
            } else if (diff[0] < 0) {
                // deletion
                byte[] phase = {1, 1};
                if(end) {
                    var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + alt[0].length, -diff[0],
                            Arrays.copyOfRange(ref, alt[0].length, alt[0].length-diff[0]), new FlexSeq[]{new FlexSeq()},
                            phase, true, var.getVar_id(), ".", ""));
                }else {
                    var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, -diff[0],
                            Arrays.copyOfRange(ref, 0, -diff[0]), new FlexSeq[]{new FlexSeq()},
                            phase, true, var.getVar_id(), ".", ""));
                }
            }
        } else {
            for (int a = 0; a < alt.length; a++) {
                if (diff[a] > 0) {
                    // insertion
                    byte[] phase = {0, 0};
                    phase[a] = 1;
                    if(end) {
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + ref.length, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[a], 0, diff[a]))},
                                phase, true, var.getVar_id(), ".", ""));
                    }else{
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 0, new byte[0],
                                new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt[a], 0, diff[a]))},
                                phase, true, var.getVar_id(), ".", ""));
                    }
                } else if (diff[a] < 0) {
                    // deletion
                    byte[] phase = {0, 0};
                    phase[a] = 1;
                    if(end) {
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + alt[a].length, -diff[a],
                                Arrays.copyOfRange(ref, alt[a].length, alt[a].length-diff[a]), new FlexSeq[]{new FlexSeq()},
                                phase, true, var.getVar_id(), ".", ""));
                    }else {
                        var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, -diff[a],
                                Arrays.copyOfRange(ref, 0, -diff[a]), new FlexSeq[]{new FlexSeq()},
                                phase, true, var.getVar_id(), ".", ""));
                    }
                }
            }
        }
    }

    private ArrayList<Variant> convert_var_to_var_list(Variant var) {
        ArrayList<Variant> var_list = convert_var_to_var_list(new Variant(var), false);
        ArrayList<Variant> var_list_end = convert_var_to_var_list(new Variant(var), true);
        if (var_list_end.size() < var_list.size()) {
            var_list = var_list_end;
        }
        return var_list;
    }

    //if end = true, we add indels to the end
    private ArrayList<Variant> convert_var_to_var_list(Variant var, boolean end) {
        ArrayList<Variant> var_list = new ArrayList<Variant>();

        //System.err.println("pat|mat: " + var.paternal() +"|"+ var.maternal());
        // if the variant is an MNP or SNP, break it dooooownnn

        boolean no_split = false;
        if (var.getType() == Variant.OverallType.SNP) {
            no_split = true;
        }
        if (var.paternal() == 0 && var.maternal() == 0) {
            no_split = true;
        }
        if (var.paternal() > 0 && var.getAlt(var.paternal()).getType() != FlexSeq.Type.SEQ) {
            no_split = true;
        }
        if (var.maternal() > 0 && var.getAlt(var.maternal()).getType() != FlexSeq.Type.SEQ) {
            no_split = true;
        }
        if (var.paternal() > 0 && var.getAlt(var.paternal()).length() == 0 && var.getRef().length == 0) {
            no_split = true;
        }
        if (var.maternal() > 0 && var.getAlt(var.maternal()).length() == 0 && var.getRef().length == 0) {
            no_split = true;
        }

        if (no_split) {
            var_list.add(var);
            return var_list;
        }


        if (var.getType(var.paternal()) != Variant.Type.Reference
                && var.getType(var.maternal()) != Variant.Type.Reference) {

            int[] allele = {var.get_allele(0), var.get_allele(1)};
            byte[][] alt = {var.getAlt(allele[0]).getSeq(), var.getAlt(allele[1]).getSeq()};
            byte[] ref = var.getRef();
            int curr_pos = var.position();

            // modify positions based on if ref matches alt
            int[] match_len = {0, 0};
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < Math.min(ref.length, alt[i].length); j++) {
                    if (alt[i][j] == ref[j]) {
                        match_len[i]++;
                    } else {
                        break;
                    }
                }
            }
            int min_match_len = Math.min(match_len[0], match_len[1]);
            //System.err.println("min_match_len: " + min_match_len);

            if (min_match_len > 0) {
                ref = Arrays.copyOfRange(ref, min_match_len, ref.length);
                for (int i = 0; i < 2; i++) {
                    alt[i] = Arrays.copyOfRange(alt[i], min_match_len, alt[i].length);
                }
                curr_pos += min_match_len;
            }

            int[] diff = {alt[0].length - ref.length, alt[1].length - ref.length};

            add_indels(var_list, diff, ref, alt, var, curr_pos,end);

            for (int i = 0; i < ref.length; i++, curr_pos++) {

                int[] idx = new int[2];
                if(end){
                    for(int j = 0;j<2;j++) {
                        if(i < ref.length + diff[j]){
                            idx[j] = i;
                        }else{
                            idx[j] = -1; // we are into deleted bases
                        }
                    }
                }else{
                    for(int j = 0;j<2;j++) {
                        idx[j] = i + diff[j];
                    }
                }

                if (idx[0] < 0 && idx[1] < 0) {
                    // both deleted
                } else if (idx[0] >= 0 && idx[1] < 0 && alt[0][idx[0]] != ref[i]) {
                    // one deleted, hence the other is homozygous
                    byte[] phase = {1, 1};
                    var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 1, new byte[]{ref[i]},
                            new FlexSeq[]{new FlexSeq(alt[0][idx[0]])}, phase, true, var.getVar_id(), ".", ""));
                } else if (idx[0] < 0 && idx[1] >= 0 && alt[1][idx[1]] != ref[i]) {
                    // one deleted, hence the other is homozygous
                    byte[] phase = {1, 1};
                    var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 1, new byte[]{ref[i]},
                            new FlexSeq[]{new FlexSeq(alt[1][idx[1]])}, phase, true, var.getVar_id(), ".", ""));
                } else if (idx[0] >= 0 && idx[1] < 0 && alt[0][idx[0]] == ref[i]) {
                    // ref call with del
                } else if (idx[0] < 0 && idx[1] >= 0 && alt[1][idx[1]] == ref[i]) {
                    // ref call with del
                } else if (alt[0][idx[0]] == ref[i] && alt[1][idx[1]] == ref[i]) {
                    // ref call
                } else if (alt[0][idx[0]] == alt[1][idx[1]]) {
                    // homozygous
                    byte[] phase = {1, 1};
                    var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 1, new byte[]{ref[i]},
                            new FlexSeq[]{new FlexSeq(alt[0][idx[0]])}, phase, true, var.getVar_id(), ".", ""));

                } else if (alt[0][idx[0]] != ref[i] && alt[1][idx[1]] != ref[i]) {
                    // het but both alt
                    byte[] phase = {1, 2};
                    var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 1, new byte[]{ref[i]},
                            new FlexSeq[]{new FlexSeq(alt[0][idx[0]]), new FlexSeq(alt[1][idx[1]])},
                            phase, true, var.getVar_id(), ".", ""));
                } else {
                    // het with one ref
                    for (int a = 0; a < 2; a++) {
                        if (alt[a][idx[a]] != ref[i]) {
                            byte[] phase = {0, 0};
                            phase[a] = 1;
                            var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 1, new byte[]{ref[i]},
                                    new FlexSeq[]{new FlexSeq(alt[a][idx[a]])}, phase, true, var.getVar_id(), ".", ""));
                        }
                    }
                }

            }

            var.set_allele(0, (byte) 0); // set to reference
            var.set_allele(1, (byte) 0); // set to reference

        } else {
            for (int a = 0; a < 2; a++) {
                int allele = var.get_allele(a);
                if (var.getType(allele) == Variant.Type.Complex
                        || var.getType(allele) == Variant.Type.MNP
                        || var.getType(allele) == Variant.Type.SNP) {
                    byte[] alt = var.getAlt(allele).getSeq();
                    byte[] ref = var.getRef();
                    int curr_pos = var.position();
                    int diff = alt.length - ref.length;

                    // add insertions or deletions for complex variants

                    if (diff > 0) {
                        // insertion
                        byte[] phase = {0, 0};
                        phase[a] = 1;
                        if(end){
                            var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + ref.length, 0, new byte[0],
                                    new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt, 0, diff))},
                                    phase, true, var.getVar_id(), ".", ""));
                        }else {
                            var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 0, new byte[0],
                                    new FlexSeq[]{new FlexSeq(Arrays.copyOfRange(alt, 0, diff))},
                                    phase, true, var.getVar_id(), ".", ""));
                        }
                    } else if (diff < 0) {
                        // deletion
                        byte[] phase = {0, 0};
                        phase[a] = 1;
                        if(end){
                            var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos + alt.length, -diff,
                                    Arrays.copyOfRange(ref, alt.length, alt.length-diff),
                                    new FlexSeq[]{new FlexSeq()},
                                    phase, true, var.getVar_id(), ".", ""));
                        }else {
                            var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, -diff,
                                    Arrays.copyOfRange(ref, 0, -diff),
                                    new FlexSeq[]{new FlexSeq()},
                                    phase, true, var.getVar_id(), ".", ""));
                        }
                    }

                    for (int i = 0; i < ref.length; i++) {
                        int idx;
                        if(end) {
                            if(i < ref.length + diff) {
                                idx = i;
                            }else{
                                idx = -1; // we are in a deleted region
                            }
                        }else{
                            idx = i + diff;
                        }

                        if (idx >= 0 && alt[idx] != ref[i]) {
                            byte[] phase = {0, 0};
                            phase[a] = 1;

                            var_list.add(new Variant(var.getChr_name(), var.chromosome(), curr_pos, 1, new byte[]{ref[i]},
                                    new FlexSeq[]{new FlexSeq(alt[idx])}, phase, true, var.getVar_id(), ".", ""));
                        }

                        curr_pos++;
                    }

                    var.set_allele(a, (byte) 0); // set to reference
                }
            }
        }

        if (!var.isRef()) {
            var_list.add(var);
        }

        return var_list;
    }

    private void run(String[] args) {
        String usage = "VCFcompare true_vcf new_vcf out_prefix overlap_ratio wiggle [bed_file]\n" +
                "Reciprocal overlap\n";
        // TODO use args4j or something similar
        boolean compare_genotypes = false;

        if (args.length < 5 || args.length > 6) {
            System.err.println(usage);
            System.exit(1);
        }

        // these are the statistics we "ideally" want to collect
        // number of variants correct (either genotype) (for each type)
        // number homozygous correct (for each type)
        // number heterozygous correct (for each type)
        // number homozygous genotype correct (for each type)
        // number heterozyous genotype correct (for each type)

        // load true VCF into interval tree
        System.err.println("Load True VCF");

        String true_vcf_filename = args[0];
        String new_vcf_filename = args[1];
        String out_prefix = args[2];
        Double overlap_ratio = Double.parseDouble(args[3]);
        int wiggle = Integer.parseInt(args[4]);

        String bed_filename = "";
        BedFile intersector = null;

        if (args.length == 6) {
            bed_filename = args[5];
            intersector = new BedFile(bed_filename);
        }

        /**
         * This is just for outputting to JSON
         */
        class output_class {
            CompareParams params;
            EnumStatsRatioCounter<Variant.OverallType> num_true_correct;

            output_class(CompareParams params, EnumStatsRatioCounter<Variant.OverallType> num_true_correct) {
                this.params = params;
                this.num_true_correct = num_true_correct;
            }

            output_class() {
            }

            public CompareParams getParams() {
                return params;
            }

            public void setParams(CompareParams params) {
                this.params = params;
            }

            public EnumStatsRatioCounter<Variant.OverallType> getNum_true_correct() {
                return num_true_correct;
            }

            public void setNum_true_correct(EnumStatsRatioCounter<Variant.OverallType> num_true_correct) {
                this.num_true_correct = num_true_correct;
            }
        }


        output_class output_blob = new output_class();

        output_blob.setParams(new CompareParams());
        output_blob.getParams().setBed_filename(bed_filename);
        output_blob.getParams().setNew_vcf_filename(new_vcf_filename);
        output_blob.getParams().setOverlap_percent(overlap_ratio);
        output_blob.getParams().setTrue_vcf_filename(true_vcf_filename);
        output_blob.getParams().setWiggle(wiggle);

        VCFparser true_parser = new VCFparser(true_vcf_filename, null, false);

        // allow duplicates, this is needed because insertions don't actually take up a location
        chrST<Variant> true_store = new chrST<Variant>(true);
        int num_read = 0;
        int num_added = 0;

        // this is for the original variants
        // it stores the total length of the original variant in bases
        // Still check for validation of canonical full variants
        ArrayList<Integer> full_validated_total = new ArrayList<Integer>();
        ArrayList<Variant> true_var_list = new ArrayList<Variant>();

        // For each true variant, if the number of bases validated is over a certain threshold
        // call it correct
        output_blob.setNum_true_correct(new EnumStatsRatioCounter<Variant.OverallType>());

        // For called variants, break down into canonical ones and count based on that
        // if any called variant overlaps a complex variant or MNP, count it as "complex"
        // otherwise, simple count them in their canonical forms


        // store true variants as canonical ones, but remember original form
        while (true_parser.hasMoreInput()) {
            Variant var = true_parser.parseLine();
            if (var == null) {
                System.err.println("skip line");
                continue;
            }

            Genotypes geno = var.getGeno();

            if (!geno.isNonRef()) {
                continue;
            }

            String chr_name = var.getChr_name();
            Variant.OverallType orig_type = var.getType();

            // determine max variant region
            // when comparing genotypes, we need to individually compare
            // to make sure they really overlap

            ArrayList<Variant> var_list = convert_var_to_var_list(new Variant(var));

            int total_len = 0;
            double max_len = 0;

            // add to interval tree
            for (Variant curr_var : var_list) {
                int curr_len = curr_var.max_len();

                if (curr_len > max_len) {
                    max_len = curr_len;
                }

                total_len += curr_len;

                Interval1D curr_var_reg = curr_var.get_geno_var_interval();
                curr_var.idx = num_added;
                curr_var.full_idx = num_read;
                curr_var.original_type = orig_type;

                //if(curr_var.original_type == Variant.OverallType.Complex){
                //    System.err.println("Set to complex... ");
                //}

                true_store.put(chr_name, curr_var_reg, curr_var);
                num_added++;
            }

            if (total_len >= 50 && max_len / total_len >= overlap_ratio && var_list.size() > 1) {
                // in this case we break down the variant into canoical forms since
                // the original variant was probably a large deletion with a small insertion

                for (Variant curr_var : var_list) {
                    int curr_len = curr_var.max_len();
                    full_validated_total.add(curr_len);
                    true_var_list.add(curr_var);
                    num_read++;
                }

            } else {
                full_validated_total.add(total_len);
                true_var_list.add(var);
                num_read++;
            }
        }

        log.info("Num read:  " + num_read);
        log.info("Num added: " + num_added);
        log.info("Num nodes: " + true_store.size());

        // this is for the split variants
        // set to true if the canonical original variant was validated true
        BitSet validated_true = new BitSet(num_added);

        // this is for the original variants
        // count of the number of bases validated for the original variant
        int[] full_validated_count = new int[num_read];

        // generate the output files
        PrintWriter TP_writer = null;
        PrintWriter FP_writer = null;
        PrintWriter FN_writer = null;
        PrintWriter JSON_writer = null;
        try {
            TP_writer = new PrintWriter(out_prefix + "_TP.vcf", "UTF-8");
            FP_writer = new PrintWriter(out_prefix + "_FP.vcf", "UTF-8");
            FN_writer = new PrintWriter(out_prefix + "_FN.vcf", "UTF-8");
            JSON_writer = new PrintWriter(out_prefix + "_report.json", "UTF-8");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


        // for this case we add to false positives if the variant is not validated.
        // However, do don't add to true positives, those that computed later

        log.info("Load New VCF");
        int num_new_vars = 0;
        // iterate over new VCF and collect stats
        VCFparser new_parser = new VCFparser(new_vcf_filename, null, false);

        while (new_parser.hasMoreInput()) {
            Variant var = new_parser.parseLine();

            if (var == null) {
                // System.err.println("Bad variant or not a variant line");
                continue;
            }

            Genotypes geno = var.getGeno();

            String chr_name = var.getChr_name();
            Interval1D var_reg = var.get_geno_interval();

            if (!(intersector == null || intersector.contains(chr_name, var_reg))) {
                continue;
            }

            // the overall type of the called variant
            Variant.OverallType curr_var_type = var.getType();

            // if called as complex variant convert to indel+snps
            ArrayList<Variant> var_list = convert_var_to_var_list(new Variant(var));

            double total_len = 0;
            double validated_len = 0;
            double max_len = 0;

            for (Variant curr_var : var_list) {
                total_len += curr_var.max_len();
                if (max_len < curr_var.max_len()) {
                    max_len = curr_var.max_len();
                }
            }

            // split up variants that are basically one big variant and one small one
            boolean compute_as_split = false;
            if (total_len >= 50 && max_len / total_len >= overlap_ratio && var_list.size() > 1){
                compute_as_split = true;
            }

            for (Variant curr_var : var_list) {

                // get genotype
                geno = curr_var.getGeno();
                result_comparator comp = new result_comparator(true_store, overlap_ratio, wiggle);

                if (curr_var.isHom()) {

                    comp.compare_variant(curr_var, geno.geno[0], validated_true);

                    dual_idx idx;
                    if (compare_genotypes) {
                        idx = comp.isHomMatch();
                    } else {
                        idx = comp.isMatch();
                    }

                    if (idx.idx >= 0) {
                        validated_true.set(idx.idx);
                        full_validated_count[idx.full_idx] += curr_var.max_len(); // this 'should' be overlap len
                        validated_len += curr_var.max_len();
                    }else if(compute_as_split){
                        output_blob.getNum_true_correct().addFP(curr_var.getType(), var.max_len());
                        FP_writer.println(var);
                    }

                } else {
                    // het
                    boolean matched = false;
                    for (int i = 0; i < 2; i++) {
                        byte allele = geno.geno[i];
                        if (allele > 0) {
                            comp.compare_variant(curr_var, allele, validated_true);
                        }
                    }

                    dual_idx idx;
                    if (compare_genotypes) {
                        idx = comp.isHetMatch();
                    } else {
                        idx = comp.isMatch();
                    }

                    if (idx.idx >= 0) {
                        validated_true.set(idx.idx);
                        full_validated_count[idx.full_idx] += curr_var.max_len(); // this 'should' be overlap len
                        validated_len += curr_var.max_len();
                    }else if(compute_as_split){
                        output_blob.getNum_true_correct().addFP(curr_var.getType(), var.max_len());
                        FP_writer.println(var);
                    }
                }
            }

            if (!compute_as_split && validated_len < (total_len*overlap_ratio)) {
                // this is a false positive!
                output_blob.getNum_true_correct().addFP(curr_var_type, var.max_len());
                FP_writer.println(var);
            }

            num_new_vars++;
        }

        log.info("Num new variants read: " + num_new_vars);

        // read through again and compute for the true variants


        num_read = 0;
        for (Variant var : true_var_list) {

            String chr_name = var.getChr_name();
            Interval1D curr_var_reg = var.get_geno_interval();

            if (intersector == null || intersector.contains(chr_name, curr_var_reg)) {
                int total_len = full_validated_total.get(num_read);
                int validated_len = full_validated_count[num_read];

                if (validated_len >= (overlap_ratio * total_len)) {
                    // validated
                    output_blob.getNum_true_correct().addTP(var.getType(), var.max_len());
                    TP_writer.println(var);
                } else {
                    FN_writer.println(var);
                }

                output_blob.getNum_true_correct().addT(var.getType(), var.max_len());
            }
            num_read++;
        }

        // output the stats
        System.err.println(output_blob.getNum_true_correct());

        ObjectMapper mapper = new ObjectMapper();
        mapper.configure(JsonGenerator.Feature.AUTO_CLOSE_TARGET, false);

        try {
            mapper.writeValue(JSON_writer, output_blob);
        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            TP_writer.close();
            FP_writer.close();
            FN_writer.close();
            JSON_writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        log.info("Done!"); // used to record the time
    }

    class CompareParams {
        String true_vcf_filename;
        String new_vcf_filename;
        Double overlap_ratio;
        int wiggle;
        String bed_filename;

        public CompareParams() {
        }

        public CompareParams(String true_vcf_filename, String new_vcf_filename, Double overlap_ratio, int wiggle, String bed_filename) {
            this.true_vcf_filename = true_vcf_filename;
            this.new_vcf_filename = new_vcf_filename;
            this.overlap_ratio = overlap_ratio;
            this.wiggle = wiggle;
            this.bed_filename = bed_filename;
        }


        public String getTrue_vcf_filename() {
            return true_vcf_filename;
        }

        public void setTrue_vcf_filename(String true_vcf_filename) {
            this.true_vcf_filename = true_vcf_filename;
        }

        public String getNew_vcf_filename() {
            return new_vcf_filename;
        }

        public void setNew_vcf_filename(String new_vcf_filename) {
            this.new_vcf_filename = new_vcf_filename;
        }

        public Double getOverlap_percent() {
            return overlap_ratio;
        }

        public void setOverlap_percent(Double overlap_ratio) {
            this.overlap_ratio = overlap_ratio;
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

    class dual_idx {
        public int idx;
        public int full_idx;

        dual_idx(int idx, int full_idx) {
            this.idx = idx;
            this.full_idx = full_idx;
        }

        dual_idx() {
            idx = -1;
            full_idx = -1;
        }

        public boolean equals(Object obj) {
            if (obj == null)
                return false;
            if (obj == this)
                return true;
            if (!(obj instanceof dual_idx))
                return false;

            // use EqualsBuilder?
            dual_idx temp = (dual_idx) obj;
            if (idx != temp.idx) {
                return false;
            }

            if (full_idx != temp.full_idx) {
                return false;
            }

            return true;
        }

        @Override
        public String toString() {
            return "dual_idx{" +
                    "idx=" + idx +
                    ", full_idx=" + full_idx +
                    '}';
        }
    }

    class result_comparator {

        chrST<Variant> _true_store; // true variants
        double _overlap_ratio;
        boolean _overlap_complex;
        int _wiggle;

        // Results to store
        // this stores the indexes of the true variants matched
        ArrayList<dual_idx> matches_hom = new ArrayList<dual_idx>();
        ArrayList<ArrayList<dual_idx>> matches_het = new ArrayList<ArrayList<dual_idx>>(2); // matches either parent

        //type_diff type_hom;
        //type_diff[] type_het = {new type_diff(),new type_diff()};

        public result_comparator(chrST<Variant> true_store, double overlap_ratio, int wiggle) {
            _true_store = true_store;
            _overlap_ratio = overlap_ratio;
            _wiggle = wiggle;

            matches_het.add(new ArrayList<dual_idx>());
            matches_het.add(new ArrayList<dual_idx>());
            _overlap_complex = false;
        }

        public dual_idx isHomMatch() {
            if (matches_hom.size() > 0) {
                return matches_hom.get(0);
            }
            return new dual_idx();
        }

        public dual_idx isHetMatch() {
            ArrayList<dual_idx> temp = new ArrayList<dual_idx>(matches_het.get(0));
            temp.retainAll(matches_het.get(1));
            if (temp.size() > 0) {
                return temp.get(0);
            } else if (matches_het.get(0).size() > 0 || matches_het.get(1).size() > 0) {

                if (matches_het.get(0).size() > matches_het.get(1).size()) {
                    return matches_het.get(0).get(0);
                } else {
                    return matches_het.get(1).get(0);
                }

            }
            return new dual_idx();
        }

        public dual_idx isMatch() {
            dual_idx idx = isHomMatch();
            if (idx.idx >= 0) {
                return idx;
            }
            idx = isHetMatch();
            if (idx.idx >= 0) {
                return idx;
            }
            return idx;
        }

        /**
         * Only compares one allele at a time
         *  - don't match variants in the bitset
         *  - if match set the bitset
         *
         * @param var       variant we want to compare
         * @param geno      allele of the variant to compare
         * @param validated BitSet that records the true variants that have already been validated
         */
        public void compare_variant(Variant var, int geno, BitSet validated) {
            double overlap_ratio = _overlap_ratio;
            // consider type to change overlap percent
            Variant.Type type = var.getType(geno);

            String chr_name = var.getChr_name();

            Interval1D orig_inter = var.get_var_interval(geno);

            //System.err.println("Comparing: " + var);

            // sometimes MNPs are called as SNPs?
            if (type == Variant.Type.SNP) {
                // handle SNPs differently
                // require SNP content to match
                Iterable<Variant> out = _true_store.getAll(chr_name, orig_inter, 0);

                byte val = var.getAlt(geno).getSeq()[0];

                int num_matches = 0;
                if (out != null) {
                    for (Variant true_var : out) {

                        boolean has_snp = false;
                        int idx = true_var.idx;
                        int full_idx = true_var.full_idx;

                        if (true_var.original_type == Variant.OverallType.Complex) {
                            //System.err.println("Overlap complex SNP!");
                            _overlap_complex = true;
                        }

                        if (validated.get(idx)) {
                            // skip ones already validated
                            continue;
                        }

                        // check genotype
                        if (true_var.isHom()) {
                            // position is correct, check genotype
                            if (true_var.getType(true_var.paternal()) == Variant.Type.SNP
                                    && var.position() == true_var.position()) {
                                if (val == true_var.getAlt(true_var.paternal()).getSeq()[0]) {
                                    matches_hom.add(new dual_idx(idx, full_idx));
                                }
                                has_snp = true;
                            }
                        } else {
                            for (int parent = 0; parent < 2; parent++) {
                                int allele = true_var.get_allele(parent);
                                if (allele > 0) {
                                    //byte[] alt = true_var.getAlt(allele).getSeq();
                                    if (true_var.getType(allele) == Variant.Type.SNP
                                            && var.position() == true_var.position()) {
                                        if (val == true_var.getAlt(allele).getSeq()[0]) {
                                            matches_het.get(parent).add(new dual_idx(idx, full_idx));
                                            //type_het[parent].update(Variant.Type.SNP, 0);
                                        }
                                        has_snp = true;
                                    }
                                }
                            }
                        }

                        if (has_snp) {
                            num_matches++;
                        }
                    }

                    if (num_matches > 1) {
                        log.info("Something strange, multiple SNP matches in true set: " + num_matches);
                    }
                }


            } else {
                // the rest
                Interval1D wiggle_inter = new Interval1D(orig_inter.low - _wiggle, orig_inter.high + _wiggle);
                Iterable<Variant> out = _true_store.getAll(chr_name, wiggle_inter, 0);


                if (out == null) {
                    // nothing found
                    return;
                }

                for (Variant true_var : out) {

                    int idx = true_var.idx;
                    int full_idx = true_var.full_idx;

                    if (true_var.original_type == Variant.OverallType.Complex) {
                        //System.err.println("Overlap complex Indel!");
                        _overlap_complex = true;
                    }

                    if (validated.get(idx)) {
                        // skip ones already validated
                        //System.err.println("Skip..." + idx);
                        continue;
                    }

                    for (int parent = 0; parent < 2; parent++) {
                        if (true_var.isHom() && parent == 1) {
                            break;
                        }

                        int allele = true_var.get_allele(parent);

                        if (allele == 0) {
                            continue;
                        }

                        if (type != true_var.getType(allele)) {
                            // need type to be the same
                            continue;
                        }


                        // check if the variant interval matches
                        if (orig_inter.intersects(true_var.get_var_interval(allele), overlap_ratio, _wiggle)) {
                            // it matches an allele!
                            // now check alternate allele length
                            int alt_len = var.getAlt(geno).length(); // TODO ignore copy number for now
                            int true_alt_len = true_var.getAlt(allele).length();
                            double ratio = (alt_len > 0) ? (true_alt_len / (double) alt_len) : 1.0;
                            double min_ratio = Math.min(ratio, 1 / ratio);
                            if (min_ratio >= overlap_ratio || Math.abs(alt_len - true_alt_len) < _wiggle) {
                                // yay, it is a match!
                                if (true_var.isHom()) {
                                    matches_hom.add(new dual_idx(idx, full_idx));
                                } else {
                                    matches_het.get(parent).add(new dual_idx(idx, full_idx));
                                }
                            }
                        }
                    }
                }

            }


        }
    }

}
