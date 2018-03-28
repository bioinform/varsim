package com.bina.varsim.types.stats;

import com.bina.varsim.constants.Constant;
import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * This is for recording values in bins of various sizes, the bins are hard coded for now
 */
public class StatsRatioRecord {
    @JsonProperty(value = "bin_counts")
    private RatioRecord[] binCounts; // the last bin is for anything larger, this is the number correct
    @JsonProperty(value = "sum_count")
    private RatioRecord sumCount;
    @JsonProperty(value = "svSumCount")
    private RatioRecord svSumCount;
    @JsonProperty(value = "sum_per_base_count")
    private RatioRecord sumPerBaseCount;

    @JsonProperty(value = "bin_breaks")
    private int[] binBreaks = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 29, 39, 49, 99,
            199, 399, 799, 1599, 3199, 6399, 12799, 25599, 51199, 102399, 500000, 1000000};

    final private int svlen;
    StatsRatioRecord(int svlen) {
        binCounts = new RatioRecord[binBreaks.length + 1];

        if (binBreaks.length > 0) {
            binCounts[0] = new RatioRecord(1, binBreaks[0]);
        } else {
            binCounts[0] = new RatioRecord(1, -1);
        }

        for (int i = 1; i < binCounts.length; i++) {
            if (i < binBreaks.length) {
                binCounts[i] = new RatioRecord(binBreaks[i - 1] + 1, binBreaks[i]);
            } else {
                binCounts[i] = new RatioRecord(binBreaks[i - 1] + 1, -1);
            }
        }

        sumCount = new RatioRecord();
        svSumCount = new RatioRecord();
        sumPerBaseCount = new RatioRecord();
        this.svlen = svlen;
    }

    public void addTP(int val) {
        sumPerBaseCount.incTP(val);
        sumCount.incTP();
        if (val >= svlen) {
            svSumCount.incTP();
        }
        for (int i = 0; i < binBreaks.length; i++) {
            if (val <= binBreaks[i]) {
                binCounts[i].incTP();
                return;
            }
        }
        binCounts[binBreaks.length].incTP();
    }

    public void addFP(int val) {
        sumPerBaseCount.incFP(val);
        sumCount.incFP();
        if (val >= svlen) {
            svSumCount.incFP();
        }
        for (int i = 0; i < binBreaks.length; i++) {
            if (val <= binBreaks[i]) {
                binCounts[i].incFP();
                return;
            }
        }
        binCounts[binBreaks.length].incFP();
    }

    public void addT(int val, int referenceBases) {
        sumPerBaseCount.incT(referenceBases);
        sumCount.incT();
        if (val >= svlen) {
            svSumCount.incT();
        }
        for (int i = 0; i < binBreaks.length; i++) {
            if (val <= binBreaks[i]) {
                binCounts[i].incT();
                return;
            }
        }
        binCounts[binBreaks.length].incT();
    }

    /**
     * This only computes it for sum_per_base_count
     * @param numNonNReferenceBases
     */
    public void computeTN(int numNonNReferenceBases){
        int conditionNegative = numNonNReferenceBases - sumPerBaseCount.getT();
        sumPerBaseCount.setTN(conditionNegative - sumPerBaseCount.getFP());
    }

    public RatioRecord[] getBinCounts() {
        return binCounts;
    }

    public void setBinCounts(RatioRecord[] binCounts) {
        this.binCounts = binCounts;
    }

    public RatioRecord getSumCount() {
        return sumCount;
    }

    public void setSumCount(RatioRecord sumCount) {
        this.sumCount = sumCount;
    }

    public int[] getBinBreaks() {
        return binBreaks;
    }

    public RatioRecord getSumPerBaseCount() {
        return sumPerBaseCount;
    }

    public void setSumPerBaseCount(RatioRecord sumPerBaseCount) {
        this.sumPerBaseCount = sumPerBaseCount;
    }

    public String toString() {
        return toString(binBreaks[binBreaks.length - 1] + 1);
    }

    public String toString(int max_len) {
        StringBuilder sb = new StringBuilder();

        // write overall TPR/FDR;
        sb.append("TPR,FDR,TP,FP,T,F1:\n");
        sb.append("ALL");
        sb.append(':');
        sb.append(sumCount);
        sb.append('\n');
        sb.append("[>=" + svlen + "]");
        sb.append(':');
        sb.append(svSumCount);
        sb.append('\n');

        for (RatioRecord bin_count : binCounts) {

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
