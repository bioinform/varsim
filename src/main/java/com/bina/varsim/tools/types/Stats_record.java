package com.bina.varsim.tools.types;

import com.bina.varsim.constants.Constant;

/**
 * @author johnmu
 */

public class Stats_record {
    int[] bin_counts; // the last bin is for anything larger
    int total_count;
    int sv_total_count;
    int[] bin_breaks = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 29, 39, 49, 99,
            199, 399, 799, 1599, 3199, 6399, 12799, 25599, 51199, 102399, 500000, 1000000};
    int num_bins;

    Stats_record() {
        num_bins = bin_breaks.length + 1;
        bin_counts = new int[num_bins];
        total_count = 0;
        sv_total_count = 0;
    }

    /**
     * Increments the bin that val belongs in
     *
     * @param val value to be added
     */
    public void add(int val) {
        total_count++;
        if (val >= Constant.SVLEN) sv_total_count++;

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

    public int getsvTotal_count() {
        return sv_total_count;
    }

    public String toString() {
        return toString(bin_breaks[bin_breaks.length - 1] + 1);
    }

    public String toString(int max_len) {
        StringBuilder sb = new StringBuilder();
        sb.append("Total: ").append(getTotal_count()).append("\n");
        sb.append("Total (>=" + Constant.SVLEN + "): ").append(getsvTotal_count()).append("\n");
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
