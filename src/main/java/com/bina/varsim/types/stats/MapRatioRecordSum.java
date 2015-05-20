package com.bina.varsim.types.stats;

import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

/**
 * This maps variant type to a ratiorecordsum
 */
public class MapRatioRecordSum {

    // use tree map so that the JSON is ordered...
    private TreeMap<String, RatioRecordSum> data;
    private RatioRecordSum all_data; // this records regardless of type

    public MapRatioRecordSum() {
        data = new TreeMap<>();
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
