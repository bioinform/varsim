package com.bina.varsim.types.stats;

import java.util.ArrayList;

/**
 * Stores list of ration records that are filtered by some increasing value, eg. MAPQ
 * The filter value is stored in the "upper" field
 */
public class RatioRecordSum {
    private ArrayList<RatioRecord> data; // accending order
    private int T; // number of true

    public RatioRecordSum() {
        data = new ArrayList<>();
        T = 0;
    }

    private void resize(int val) {
        int max_val = data.size() - 1;
        if (val > max_val) {
            int num_add = val - max_val;
            for (int i = 0; i < num_add; i++) {
                data.add(new RatioRecord(0, max_val + i + 1));
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

    public ArrayList<RatioRecord> getData() {
        return data;
    }

    public void setData(ArrayList<RatioRecord> data) {
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
            RatioRecord rec = data.get(i);
            rec.setT(T);
            sb.append(i);
            sb.append(",");
            sb.append(rec);
            sb.append('\n');
        }
        return sb.toString();
    }
}
