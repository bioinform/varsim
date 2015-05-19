package com.bina.varsim.tools.stats;

import com.bina.varsim.constants.Constant;

/**
 * This is for recording values in bins of various sizes, the bins are hard coded for now
 */
public class StatsRatioRecord {
    private RatioRecord[] bin_counts; // the last bin is for anything larger, this is the number correct
    private RatioRecord sum_count;
    private RatioRecord svSumCount;

    private int[] bin_breaks = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 29, 39, 49, 99,
            199, 399, 799, 1599, 3199, 6399, 12799, 25599, 51199, 102399, 500000, 1000000};

    StatsRatioRecord() {
        bin_counts = new RatioRecord[bin_breaks.length + 1];

        if (bin_breaks.length > 0) {
            bin_counts[0] = new RatioRecord(1, bin_breaks[0]);
        } else {
            bin_counts[0] = new RatioRecord(1, -1);
        }

        for (int i = 1; i < bin_counts.length; i++) {
            if (i < bin_breaks.length) {
                bin_counts[i] = new RatioRecord(bin_breaks[i - 1] + 1, bin_breaks[i]);
            } else {
                bin_counts[i] = new RatioRecord(bin_breaks[i - 1] + 1, -1);
            }
        }

        sum_count = new RatioRecord();
        svSumCount = new RatioRecord();
    }

    public void addTP(int val) {
        sum_count.incTP();
        if (val >= Constant.SVLEN) {
            svSumCount.incTP();
        }
        for (int i = 0; i < bin_breaks.length; i++) {
            if (val <= bin_breaks[i]) {
                bin_counts[i].incTP();
                return;
            }
        }
        bin_counts[bin_breaks.length].incTP();
    }

    public void addFP(int val) {
        sum_count.incFP();
        if (val >= Constant.SVLEN) {
            svSumCount.incFP();
        }
        for (int i = 0; i < bin_breaks.length; i++) {
            if (val <= bin_breaks[i]) {
                bin_counts[i].incFP();
                return;
            }
        }
        bin_counts[bin_breaks.length].incFP();
    }

    public void addT(int val) {
        sum_count.incT();
        if (val >= Constant.SVLEN) {
            svSumCount.incT();
        }
        for (int i = 0; i < bin_breaks.length; i++) {
            if (val <= bin_breaks[i]) {
                bin_counts[i].incT();
                return;
            }
        }
        bin_counts[bin_breaks.length].incT();
    }

    public RatioRecord[] getBin_counts() {
        return bin_counts;
    }

    public void setBin_counts(RatioRecord[] bin_counts) {
        this.bin_counts = bin_counts;
    }

    public RatioRecord getSum_count() {
        return sum_count;
    }

    public void setSum_count(RatioRecord sum_count) {
        this.sum_count = sum_count;
    }

    public int[] getBin_breaks() {
        return bin_breaks;
    }

    public String toString() {
        return toString(bin_breaks[bin_breaks.length - 1] + 1);
    }

    public String toString(int max_len) {
        StringBuilder sb = new StringBuilder();

        // write overall TPR/FDR;
        sb.append("TPR,FDR,TP,FP,T:\n");
        sb.append("ALL");
        sb.append(':');
        sb.append(sum_count);
        sb.append('\n');
        sb.append("[>=" + Constant.SVLEN + "]");
        sb.append(':');
        sb.append(svSumCount);
        sb.append('\n');

        for (RatioRecord bin_count : bin_counts) {

            if (bin_count.getLower() > max_len) {
                break;
            }

            if (bin_count.isEmpty()) {
                continue;
            }

            sb.append(bin_count.rangeStr());
            sb.append(':');
            sb.append(bin_count.toString());
            sb.append('\n');
        }
        return sb.toString();
    }
}
