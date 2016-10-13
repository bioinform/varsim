package com.bina.varsim.types.stats;

import com.bina.varsim.types.variant.INonReference;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by johnmu on 3/28/14.
 */


/**
 * Maps variant type to a stats_ratio_record
 *
 * @param <Value> This is usually a string
 */
public class EnumStatsRatioCounter<Value extends Enum & INonReference> {

    private TreeMap<Value, StatsRatioRecord> data;
    private StatsRatioRecord all_data; // this records regardless of type

    public EnumStatsRatioCounter() {
        //why TreeMap used rather than HashMap?
        data = new TreeMap<>();
        all_data = new StatsRatioRecord();
    }

    //TODO: rename incTP to something easier to understand
    public void incTP(Value a, int len) {
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

    public void incFP(Value a, int len) {
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

    public void incT(Value a, int len) {
        StatsRatioRecord count = data.get(a);
        if (count != null) {
            count.addT(len, a.isNonReference() ? 0 : len);
        } else {
            StatsRatioRecord contents = new StatsRatioRecord();
            contents.addT(len, a.isNonReference() ? 0 : len);
            data.put(a, contents);
        }

        all_data.addT(len, a.isNonReference() ? 0 : len);
    }

    public TreeMap<Value, StatsRatioRecord> getData() {
        return data;
    }

    public void setData(TreeMap<Value, StatsRatioRecord> data) {
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
