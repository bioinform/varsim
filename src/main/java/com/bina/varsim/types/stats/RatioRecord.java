package com.bina.varsim.types.stats;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * Stores everything required to compute precision and recall
 */
public class RatioRecord {
    @JsonProperty(value = "_TP")
    private int TP = 0; // True Positive
    @JsonProperty(value = "_FP")
    private int FP = 0; // False Positive
    @JsonProperty(value = "_TN")
    private int TN = 0; // True Positive
    @JsonProperty(value = "_FN")
    private int FN = 0; // False Positive
    @JsonProperty(value = "T")
    private int T = 0; // True = True Positive + False Negative

    // these are inclusive
    @JsonProperty(value = "lower")
    private int lower = -1; // -1 means negative infinity or unknown
    @JsonProperty(value = "upper")
    private int upper = -1; // -1 means positive infinity or unknown

    public RatioRecord() {
    }

    public RatioRecord(int lower, int upper) {
        this.lower = lower;
        this.upper = upper;
    }

    public void incTP() {
        TP++;
    }

    public void incTP(final int TP) {
        this.TP += TP;
    }

    public void incTN() {
        TN++;
    }

    public void incFP() {
        FP++;
    }

    public void incFP(final int FP) {
        this.FP += FP;
    }

    public void incFN() {
        FN++;
    }

    public void incT() {
        T++;
    }

    public void incT(final int T) {
        this.T += T;
    }

    // recall
    @JsonProperty(value = "tpr")
    public double getTPR() {
        return TP / ((double) T);
    }
    @JsonProperty(value = "fdr")
    public double getFDR() {
        return FP / ((double) TP + FP);
    }

    // precision
    @JsonProperty(value = "ppv")
    public double getPPV() {
        return TP / ((double) TP + FP);
    }

    // specificiy
    @JsonProperty(value = "spc")
    public double getSPC() {
        return TN / ((double) TN + FP);
    }

    @JsonProperty(value = "empty")
    public boolean isEmpty() {
        return TP == 0 && FP == 0 && T == 0;
    }

    public int getTN() {
        return TN;
    }

    public void setTN(int TN) {
        this.TN = TN;
    }

    public int getFN() {
        return FN;
    }

    public void setFN(int FN) {
        this.FN = FN;
    }

    public int getTP() {
        return TP;
    }

    public void setTP(int TP) {
        this.TP = TP;
    }

    public int getFP() {
        return FP;
    }

    public void setFP(int FP) {
        this.FP = FP;
    }

    public int getT() {
        return T;
    }

    public void setT(int t) {
        this.T = t;
    }

    public int getLower() {
        return lower;
    }

    public void setLower(int lower) {
        this.lower = lower;
    }

    public int getUpper() {
        return upper;
    }

    public void setUpper(int upper) {
        this.upper = upper;
    }

    @JsonProperty(value = "f1")
    public double getF1() {
        double precision = 1 - getFDR();
        double recall = getTPR();

        return 2 * (precision * recall) / (precision + recall);
    }

    public String rangeStr() {
        StringBuilder sb = new StringBuilder();

        sb.append('[');
        if (lower < 0) {
            sb.append("-inf");
        } else {
            sb.append(lower);
        }
        sb.append(',');
        if (upper < 0) {
            sb.append("inf");
        } else {
            sb.append(upper);
        }
        sb.append(']');

        return sb.toString();
    }

    public String toString() {
        return String.format("%.4f,%.4f,%d,%d,%d,%.4f", 100 * getTPR(), 100 * getFDR(), TP, FP, T, getF1());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof RatioRecord)) return false;

        RatioRecord that = (RatioRecord) o;

        if (FN != that.FN) return false;
        if (FP != that.FP) return false;
        if (T != that.T) return false;
        if (TN != that.TN) return false;
        if (TP != that.TP) return false;
        if (lower != that.lower) return false;
        if (upper != that.upper) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = TP;
        result = 31 * result + FP;
        result = 31 * result + TN;
        result = 31 * result + FN;
        result = 31 * result + T;
        result = 31 * result + lower;
        result = 31 * result + upper;
        return result;
    }
}
