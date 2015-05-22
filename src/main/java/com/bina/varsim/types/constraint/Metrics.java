package com.bina.varsim.types.constraint;

/**
 * Created by johnmu on 5/20/15.
 */
public enum Metrics {
    F1, TPR, FDR;

    public static Metrics fromString(String s) {
        for (Metrics type : Metrics.values()) {
            if (s.equalsIgnoreCase(type.toString())) {
                return type;
            }
        }
        throw new RuntimeException(String.format("Invalid metric type: %s\n Valid ones are %s", s, Metrics.allToString()));
    }

    public static String allToString() {
        StringBuilder sb = new StringBuilder();
        for (Metrics type : Metrics.values()) {
            sb.append(type.toString()).append(',');
        }
        return sb.toString();
    }
}
