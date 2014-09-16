package com.binatechnologies.varsim;

import java.util.*;

/**
 * Created by johnmu on 3/28/14.
 */

/**
 * Stores everything required to compute precision and recall
 */
class ratio_record {
    private int _TP = 0; // True Positive
    private int _FP = 0; // False Positive
    private int _TN = 0; // True Positive
    private int _FN = 0; // False Positive
    private int _T = 0; // True = True Positive + False Negative

    // these are inclusive
    private int lower = -1; // -1 means negative infinity or unknown
    private int upper = -1; // -1 means positive infinity or unknown

    public ratio_record() {
    }

    public ratio_record(int lower, int upper) {
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
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%.4f", 100 * getTPR()));
        sb.append(',');
        sb.append(String.format("%.4f", 100 * getFDR()));
        sb.append(',');
        sb.append(_TP);
        sb.append(',');
        sb.append(_FP);
        sb.append(',');
        sb.append(_TN);
        sb.append(',');
        sb.append(_FN);
        sb.append(',');
        sb.append(_T);
        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ratio_record that = (ratio_record) o;

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


/**
 * Stores list of ration records that are filtered by some increasing value, eg. MAPQ
 * The filter value is stored in the "upper" field
 */
class RatioRecordSum {
    private ArrayList<ratio_record> data; // accending order
    private int T; // number of true

    public RatioRecordSum() {
        data = new ArrayList<ratio_record>();
        T = 0;
    }

    private void resize(int val) {
        int max_val = data.size() - 1;
        if (val > max_val) {
            int num_add = val - max_val;
            for (int i = 0; i < num_add; i++) {
                data.add(new ratio_record(0, max_val + i + 1));
            }
        }
    }

    public void incT() {
        T++;
    }

    public void incTP(int val) {
        resize(val);
        for (int i = 0; i <= val; i++) {
            data.get(i).incTP();
        }
    }

    public void incTN(int val) {
        resize(val);
        for (int i = 0; i <= val; i++) {
            data.get(i).incTN();
        }
    }

    public void incFP(int val) {
        resize(val);
        for (int i = 0; i <= val; i++) {
            data.get(i).incFP();
        }
    }

    public void incFN(int val) {
        resize(val);
        for (int i = 0; i <= val; i++) {
            data.get(i).incFN();
        }
    }

    public ArrayList<ratio_record> getData() {
        return data;
    }

    public void setData(ArrayList<ratio_record> data) {
        this.data = data;
    }

    public int getT() {
        return T;
    }

    public void setT(int t) {
        T = t;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < data.size(); i++) {
            ratio_record rec = data.get(i);
            rec.setT(T);
            sb.append(i);
            sb.append(",");
            sb.append(rec);
            sb.append('\n');
        }
        return sb.toString();
    }
}

/**
 * This maps variant type to a ratiorecordsum
 */
class MapRatioRecordSum {

    // use tree map so that the JSON is ordered...
    private TreeMap<String, RatioRecordSum> data;
    private RatioRecordSum all_data; // this records regardless of type

    public MapRatioRecordSum() {
        data = new TreeMap<String, RatioRecordSum>();
        all_data = new RatioRecordSum();
    }

    public void incTP(HashSet<String> a, int val) {
        for (String key : a) {
            RatioRecordSum count = data.get(key);
            if (count != null) {
                count.incTP(val);
            } else {
                RatioRecordSum contents = new RatioRecordSum();
                contents.incTP(val);
                data.put(key, contents);
            }
        }

        all_data.incTP(val);
    }

    public void incTP(String a, int val) {
        RatioRecordSum count = data.get(a);
        if (count != null) {
            count.incTP(val);
        } else {
            RatioRecordSum contents = new RatioRecordSum();
            contents.incTP(val);
            data.put(a, contents);
        }

        all_data.incTP(val);
    }

    public void incTN(HashSet<String> a, int val) {
        for (String key : a) {
            RatioRecordSum count = data.get(key);
            if (count != null) {
                count.incTN(val);
            } else {
                RatioRecordSum contents = new RatioRecordSum();
                contents.incTN(val);
                data.put(key, contents);
            }
        }

        all_data.incTN(val);
    }

    public void incTN(String a, int val) {
        RatioRecordSum count = data.get(a);
        if (count != null) {
            count.incTN(val);
        } else {
            RatioRecordSum contents = new RatioRecordSum();
            contents.incTN(val);
            data.put(a, contents);
        }

        all_data.incTN(val);
    }

    public void incFP(HashSet<String> a, int val) {
        for (String key : a) {
            RatioRecordSum count = data.get(key);
            if (count != null) {
                count.incFP(val);
            } else {
                RatioRecordSum contents = new RatioRecordSum();
                contents.incFP(val);
                data.put(key, contents);
            }
        }

        all_data.incFP(val);
    }

    public void incFP(String a, int val) {
        RatioRecordSum count = data.get(a);
        if (count != null) {
            count.incFP(val);
        } else {
            RatioRecordSum contents = new RatioRecordSum();
            contents.incFP(val);
            data.put(a, contents);
        }

        all_data.incFP(val);
    }


    public void incFN(HashSet<String> a, int val) {
        for (String key : a) {
            RatioRecordSum count = data.get(key);
            if (count != null) {
                count.incFN(val);
            } else {
                RatioRecordSum contents = new RatioRecordSum();
                contents.incFN(val);
                data.put(key, contents);
            }
        }

        all_data.incFN(val);
    }

    public void incFN(String a, int val) {
        RatioRecordSum count = data.get(a);
        if (count != null) {
            count.incFN(val);
        } else {
            RatioRecordSum contents = new RatioRecordSum();
            contents.incFN(val);
            data.put(a, contents);
        }

        all_data.incFN(val);
    }

    public void incT(HashSet<String> a) {
        for (String key : a) {
            RatioRecordSum count = data.get(key);
            if (count != null) {
                count.incT();
            } else {
                RatioRecordSum contents = new RatioRecordSum();
                contents.incT();
                data.put(key, contents);
            }
        }

        all_data.incT();
    }

    public void incT(String a) {
        RatioRecordSum count = data.get(a);
        if (count != null) {
            count.incT();
        } else {
            RatioRecordSum contents = new RatioRecordSum();
            contents.incT();
            data.put(a, contents);
        }

        all_data.incT();
    }

    public TreeMap<String, RatioRecordSum> getData() {
        return data;
    }

    public void setData(TreeMap<String, RatioRecordSum> data) {
        this.data = data;
    }

    public RatioRecordSum getAll_data() {
        return all_data;
    }

    public void setAll_data(RatioRecordSum all_data) {
        this.all_data = all_data;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        // Output for everything first
        sb.append(all_data.toString());

        sb.append("---------\n");

        // Output for each variant type
        for (Map.Entry<String, RatioRecordSum> entry : data.entrySet()) {
            sb.append(entry.getKey());
            sb.append('\n');
            sb.append(entry.getValue());
            sb.append('\n');
        }

        sb.append("---------\n");

        return sb.toString();
    }

}


/**
 * This is for recording values in bins of various sizes, the bins are hard coded for now
 */
class Stats_ratio_record {
    private ratio_record[] bin_counts; // the last bin is for anything larger, this is the number correct
    private ratio_record sum_count;

    private int[] bin_breaks = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 29, 39, 49, 99,
            199, 399, 799, 1599, 3199, 6399, 12799, 25599, 51199, 102399, 500000, 1000000};

    Stats_ratio_record() {
        bin_counts = new ratio_record[bin_breaks.length + 1];

        if (bin_breaks.length > 0) {
            bin_counts[0] = new ratio_record(1, bin_breaks[0]);
        } else {
            bin_counts[0] = new ratio_record(1, -1);
        }

        for (int i = 1; i < bin_counts.length; i++) {
            if (i < bin_breaks.length) {
                bin_counts[i] = new ratio_record(bin_breaks[i - 1] + 1, bin_breaks[i]);
            } else {
                bin_counts[i] = new ratio_record(bin_breaks[i - 1] + 1, -1);
            }
        }

        sum_count = new ratio_record();
    }

    public void addTP(int val) {
        sum_count.incTP();
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
        for (int i = 0; i < bin_breaks.length; i++) {
            if (val <= bin_breaks[i]) {
                bin_counts[i].incT();
                return;
            }
        }
        bin_counts[bin_breaks.length].incT();
    }

    public ratio_record[] getBin_counts() {
        return bin_counts;
    }

    public void setBin_counts(ratio_record[] bin_counts) {
        this.bin_counts = bin_counts;
    }

    public ratio_record getSum_count() {
        return sum_count;
    }

    public void setSum_count(ratio_record sum_count) {
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
        sb.append(sum_count);
        sb.append('\n');

        for (int i = 0; i < bin_counts.length; i++) {

            if (bin_counts[i].getLower() > max_len) {
                break;
            }

            if (bin_counts[i].isEmpty()) {
                continue;
            }

            sb.append(bin_counts[i].rangeStr());
            sb.append(':');
            sb.append(bin_counts[i].toString());
            sb.append('\n');
        }
        return sb.toString();
    }
}


/**
 * Maps variant type to a stats_ratio_record
 *
 * @param <Value> This is usually a string
 */
public class EnumStatsRatioCounter<Value extends Enum> {

    private HashMap<Value, Stats_ratio_record> data;
    private Stats_ratio_record all_data; // this records regardless of type

    public EnumStatsRatioCounter() {
        data = new HashMap<Value, Stats_ratio_record>();
        all_data = new Stats_ratio_record();
    }

    public void addTP(Value a, int len) {
        Stats_ratio_record count = data.get(a);
        if (count != null) {
            count.addTP(len);
        } else {
            Stats_ratio_record contents = new Stats_ratio_record();
            contents.addTP(len);
            data.put(a, contents);
        }

        all_data.addTP(len);
    }

    public void addFP(Value a, int len) {
        Stats_ratio_record count = data.get(a);
        if (count != null) {
            count.addFP(len);
        } else {
            Stats_ratio_record contents = new Stats_ratio_record();
            contents.addFP(len);
            data.put(a, contents);
        }

        all_data.addFP(len);
    }

    public void addT(Value a, int len) {
        Stats_ratio_record count = data.get(a);
        if (count != null) {
            count.addT(len);
        } else {
            Stats_ratio_record contents = new Stats_ratio_record();
            contents.addT(len);
            data.put(a, contents);
        }

        all_data.addT(len);
    }

    public HashMap<Value, Stats_ratio_record> getData() {
        return data;
    }

    public void setData(HashMap<Value, Stats_ratio_record> data) {
        this.data = data;
    }

    public Stats_ratio_record getAll_data() {
        return all_data;
    }

    public void setAll_data(Stats_ratio_record all_data) {
        this.all_data = all_data;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        // Output for everything first
        sb.append(all_data.toString());

        sb.append("---------\n");

        // Output for each variant type
        for (Map.Entry<Value, Stats_ratio_record> entry : data.entrySet()) {
            sb.append(entry.getKey());
            sb.append('\n');
            sb.append(entry.getValue());
            sb.append('\n');
        }

        sb.append("---------\n");

        return sb.toString();
    }

}
