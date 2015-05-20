package com.bina.varsim.types;


/**
 * Stores the data for each haplotype separately
 */
public class Parent_record {
    // 0 = paternal, 1 = maternal
    public final static int PATERNAL = 0;
    public final static int MATERNAL = 1;
    Type_record<VariantType>[] data;
    Type_record<VariantOverallType> overall_data;
    int total_count;

    public Parent_record() {
        data = (Type_record<VariantType>[]) new Type_record[2];
        data[PATERNAL] = new Type_record<>();
        data[MATERNAL] = new Type_record<>();
        overall_data = new Type_record<>();
        total_count = 0;
    }

    public void add(Variant var, BedFile bed_file) {
        int paternal_allele = var.getgood_paternal();
        int maternal_allele = var.getgood_maternal();

        boolean added = false;
        if (bed_file == null
                || bed_file.contains(var.getChr(), var.get_interval(paternal_allele))) {
            data[PATERNAL].add(var.getType(paternal_allele), var.max_len(paternal_allele));
            added = true;
        }

        if (bed_file == null
                || bed_file.contains(var.getChr(), var.get_interval(maternal_allele))) {
            data[MATERNAL].add(var.getType(maternal_allele), var.max_len(maternal_allele));
            added = true;
        }

        if (bed_file == null
                || bed_file.contains(var.getChr(), var.get_geno_interval())) {
            overall_data.add(var.getType(), var.max_len());
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
