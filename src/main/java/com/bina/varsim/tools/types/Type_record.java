package com.bina.varsim.tools.types;

import com.bina.varsim.types.Variant;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by johnmu on 5/19/15.
 */
public class Type_record<T extends Enum> {
    HashMap<T, Stats_record> data;

    Type_record() {
        data = new HashMap<>();
    }

    /**
     * Increments the bin (of type) that val belongs in
     *
     * @param type type of bin
     * @param val  value to be incremented
     */
    public void add(T type, int val) {
        if (data.containsKey(type)) {
            Stats_record rec = data.get(type);
            rec.add(val);
        } else {
            Stats_record rec = new Stats_record();
            rec.add(val);
            data.put(type, rec);
        }
    }

    public int getTotal_nonref() {
        int total = 0;
        for (Map.Entry<T, Stats_record> entry : data.entrySet()) {
            total += entry.getValue().getTotal_count();
        }
        return total;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<T, Stats_record> entry : data.entrySet()) {
            sb.append(entry.getKey().name()).append("\n");
            sb.append('\n');
            if (entry.getKey() == Variant.Type.SNP) {
                sb.append(entry.getValue().toString(1));
            } else {
                sb.append(entry.getValue());
            }
        }
        return sb.toString();
    }
}
