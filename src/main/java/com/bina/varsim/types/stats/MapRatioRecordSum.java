package com.bina.varsim.types.stats;

import java.util.*;

/**
 * This maps variant type to a ratiorecordsum
 */
public class MapRatioRecordSum {

    // use tree map so that the JSON is ordered...
    private EnumMap<EventTypesForStats, RatioRecordSum> data;
    private RatioRecordSum all_data; // this records regardless of type

    public MapRatioRecordSum() {
        data = new EnumMap<>(EventTypesForStats.class);
        all_data = new RatioRecordSum();
    }

    public void incStat(final Set<EventTypesForStats> a, final int val, final StatsNamespace stat) {
        switch(stat) {
            case TP:
                incTP(a, val);
                break;
            case FP:
                incFP(a, val);
                break;
            case TN:
                incTN(a, val);
                break;
            case FN:
                incFN(a, val);
                break;
            case T:
                incT(a);
                break;
            default:
        }
    }

    public void incTP(Set<EventTypesForStats> a, int val) {
        for (EventTypesForStats key : a) {
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

    public void incTN(Set<EventTypesForStats> a, int val) {
        for (EventTypesForStats key : a) {
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

    public void incFP(Set<EventTypesForStats> a, int val) {
        for (EventTypesForStats key : a) {
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

    public void incFN(Set<EventTypesForStats> a, int val) {
        for (EventTypesForStats key : a) {
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

    public void incT(Set<EventTypesForStats> a) {
        for (EventTypesForStats key : a) {
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

    public EnumMap<EventTypesForStats, RatioRecordSum> getData() {
        return data;
    }

    public void setData(EnumMap<EventTypesForStats, RatioRecordSum> data) {
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
        for (Map.Entry<EventTypesForStats, RatioRecordSum> entry : data.entrySet()) {
            sb.append(entry.getKey());
            sb.append('\n');
            sb.append(entry.getValue());
            sb.append('\n');
        }

        sb.append("---------\n");

        return sb.toString();
    }
    public enum EventTypesForStats {
        All, True_Unmapped,
        Sequence,
        Insertion,
        Deletion,
        Inversion,
        Tandem_Duplication,
        Unknown;
    }

}
