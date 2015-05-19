package com.bina.varsim.tools.stats;

import java.util.*;

/**
 * Created by johnmu on 3/28/14.
 */


/**
 * Maps variant type to a stats_ratio_record
 *
 * @param <Value> This is usually a string
 */
public class EnumStatsRatioCounter<Value extends Enum> {

    private HashMap<Value, StatsRatioRecord> data;
    private StatsRatioRecord all_data; // this records regardless of type

    public EnumStatsRatioCounter() {
        data = new HashMap<>();
        all_data = new StatsRatioRecord();
    }

    public void addTP(Value a, int len) {
        StatsRatioRecord count = data.get(a);
        if (count != null) {
            count.addTP(len);
        } else {
            StatsRatioRecord contents = new StatsRatioRecord();
            contents.addTP(len);
            data.put(a, contents);
        }

        all_data.addTP(len);
    }

    public void addFP(Value a, int len) {
        StatsRatioRecord count = data.get(a);
        if (count != null) {
            count.addFP(len);
        } else {
            StatsRatioRecord contents = new StatsRatioRecord();
            contents.addFP(len);
            data.put(a, contents);
        }

        all_data.addFP(len);
    }

    public void addT(Value a, int len) {
        StatsRatioRecord count = data.get(a);
        if (count != null) {
            count.addT(len);
        } else {
            StatsRatioRecord contents = new StatsRatioRecord();
            contents.addT(len);
            data.put(a, contents);
        }

        all_data.addT(len);
    }

    public HashMap<Value, StatsRatioRecord> getData() {
        return data;
    }

    public void setData(HashMap<Value, StatsRatioRecord> data) {
        this.data = data;
    }

    public StatsRatioRecord getAll_data() {
        return all_data;
    }

    public void setAll_data(StatsRatioRecord all_data) {
        this.all_data = all_data;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        // Output for everything first
        sb.append(all_data.toString());

        sb.append("---------\n");

        // Output for each variant type
        for (Map.Entry<Value, StatsRatioRecord> entry : data.entrySet()) {
            sb.append(entry.getKey());
            sb.append('\n');
            sb.append(entry.getValue());
            sb.append('\n');
        }

        sb.append("---------\n");

        return sb.toString();
    }

}
