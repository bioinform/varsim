package com.bina.varsim.types;


import com.bina.varsim.types.stats.Type_record;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.types.variant.VariantOverallType;
import com.bina.varsim.types.variant.VariantType;

/**
 * Stores the data for each haplotype separately
 */
public class ParentRecord {
    // 0 = paternal, 1 = maternal
    public final static int PATERNAL = 0;
    public final static int MATERNAL = 1;
    Type_record<VariantType>[] data;
    Type_record<VariantOverallType> overall_data;
    int total_count;

    public ParentRecord() {
        data = (Type_record<VariantType>[]) new Type_record[2];
        data[PATERNAL] = new Type_record<>();
        data[MATERNAL] = new Type_record<>();
        overall_data = new Type_record<>();
        total_count = 0;
    }

    public void add(Variant var, BedFile bed_file) {
        int paternal_allele = var.getGoodPaternal();
        int maternal_allele = var.getGoodMaternal();

        boolean added = false;
        if (bed_file == null
                || bed_file.contains(var.getChr(), var.getAlternativeAlleleInterval(paternal_allele))) {
            data[PATERNAL].add(var.getType(paternal_allele), var.maxLen(paternal_allele));
            added = true;
        }

        if (bed_file == null
                || bed_file.contains(var.getChr(), var.getAlternativeAlleleInterval(maternal_allele))) {
            data[MATERNAL].add(var.getType(maternal_allele), var.maxLen(maternal_allele));
            added = true;
        }

        if (bed_file == null
                || bed_file.contains(var.getChr(), var.getGenotypeUnionAlternativeInterval())) {
            overall_data.add(var.getType(), var.maxLen());
            added = true;
        }

        if (added) {
            total_count++;
        }
    }

    public String toString() {
        return "Total count: " + total_count + "\n"
                + "Paternal\n"
                + "Total: " + data[PATERNAL].getTotal_nonref() + "\n"
                + data[PATERNAL] + "\n"
                + "Maternal\n"
                + "Total: " + data[MATERNAL].getTotal_nonref() + "\n"
                + data[MATERNAL] + "\n"
                + "Overall\n"
                + "Total: " + overall_data.getTotal_nonref() + "\n"
                + overall_data + "\n";
    }
}
