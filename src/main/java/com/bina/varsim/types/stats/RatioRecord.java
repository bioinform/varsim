package com.bina.varsim.types.stats;

/**
 * Stores everything required to compute precision and recall
 */
public class RatioRecord {
    private int _TP = 0; // True Positive
    private int _FP = 0; // False Positive
    private int _TN = 0; // True Positive
    private int _FN = 0; // False Positive
    private int _T = 0; // True = True Positive + False Negative

    // these are inclusive
    private int lower = -1; // -1 means negative infinity or unknown
    private int upper = -1; // -1 means positive infinity or unknown

    public RatioRecord() {
    }

    public RatioRecord(int lower, int upper) {
        this.lower = lower;
        this.upper = upper;
    }

    public void setT(int T) {
        _T = T;
    }

    public void incTP() {
        _TP++;
    }

    public void incTN() {
        _TN++;
    }

    public void incFP() {
        _FP++;
    }

    public void incFN() {
        _FN++;
    }

    public void incT() {
        _T++;
    }

    // recall
    public double getTPR() {
        return _TP / ((double) _T);
    }

    public double getFDR() {
        return _FP / ((double) _TP + _FP);
    }

    // precision
    public double getPPV() {
        return _TP / ((double) _TP + _FP);
    }

    public boolean isEmpty() {
        return _TP == 0 && _FP == 0 && _T == 0;
    }

    public int get_TN() {
        return _TN;
    }

    public void set_TN(int _TN) {
        this._TN = _TN;
    }

    public int get_FN() {
        return _FN;
    }

    public void set_FN(int _FN) {
        this._FN = _FN;
    }

    public int get_TP() {
        return _TP;
    }

    public void set_TP(int _TP) {
        this._TP = _TP;
    }

    public int get_FP() {
        return _FP;
    }

    public void set_FP(int _FP) {
        this._FP = _FP;
    }

    public int get_T() {
        return _T;
    }

    public void set_T(int _T) {
        this._T = _T;
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
        return String.format("%.4f,%.4f,%d,%d,%d,%.4f", 100 * getTPR(), 100 * getFDR(), _TP, _FP, _T, getF1());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        RatioRecord that = (RatioRecord) o;

        if (_FN != that._FN) return false;
        if (_FP != that._FP) return false;
        if (_T != that._T) return false;
        if (_TN != that._TN) return false;
        if (_TP != that._TP) return false;
        if (lower != that.lower) return false;
        if (upper != that.upper) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = _TP;
        result = 31 * result + _FP;
        result = 31 * result + _TN;
        result = 31 * result + _FN;
        result = 31 * result + _T;
        result = 31 * result + lower;
        result = 31 * result + upper;
        return result;
    }
}
