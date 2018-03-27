package com.bina.varsim.types.stats;

import com.bina.varsim.types.variant.INonReference;
import com.fasterxml.jackson.annotation.JsonProperty;

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

    @JsonProperty(value = "data")
    private TreeMap<Value, StatsRatioRecord> data;
    @JsonProperty(value = "all_data")
    private StatsRatioRecord allData; // this records regardless of type
    final private int svlen;

    public EnumStatsRatioCounter(int svlen) {
        //why TreeMap used rather than HashMap?
        data = new TreeMap<>();
        //instantiate bins(1:1,2:2,...), some objects for stats
        this.svlen = svlen;
        allData = new StatsRatioRecord(this.svlen);
    }

    //TODO: rename incTP to something easier to understand
    public void incTP(Value a, int len) {
        StatsRatioRecord count = data.get(a);
        if (count != null) {
            count.addTP(len);
        } else {
            StatsRatioRecord contents = new StatsRatioRecord(svlen);
            contents.addTP(len);
            data.put(a, contents);
        }

        allData.addTP(len);
    }

    public void incFP(Value a, int len) {
        StatsRatioRecord count = data.get(a);
        if (count != null) {
            count.addFP(len);
        } else {
            StatsRatioRecord contents = new StatsRatioRecord(svlen);
            contents.addFP(len);
            data.put(a, contents);
        }

        allData.addFP(len);
    }

    public void incT(Value a, int len) {
        StatsRatioRecord count = data.get(a);
        if (count != null) {
            count.addT(len, a.isNonReference() ? 0 : len);
        } else {
            StatsRatioRecord contents = new StatsRatioRecord(svlen);
            contents.addT(len, a.isNonReference() ? 0 : len);
            data.put(a, contents);
        }

        allData.addT(len, a.isNonReference() ? 0 : len);
    }

    public TreeMap<Value, StatsRatioRecord> getData() {
        return data;
    }

    public void setData(TreeMap<Value, StatsRatioRecord> data) {
        this.data = data;
    }

    public StatsRatioRecord getAllData() {
        return allData;
    }

    public void setAllData(StatsRatioRecord allData) {
        this.allData = allData;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        // Output for everything first
        sb.append(allData.toString());

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
