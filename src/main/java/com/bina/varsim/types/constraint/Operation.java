package com.bina.varsim.types.constraint;

/**
 * Created by johnmu on 5/20/15.
 */
public enum Operation {
    GT, // Greater than
    LT; // Less than

    public static Operation fromString(String s) {
        for (Operation type : Operation.values()) {
            if (s.equalsIgnoreCase(type.toString())) {
                return type;
            }
        }
        throw new RuntimeException(String.format("Invalid operation type: %s\n Valid ones are %s", s, Operation.allToString()));
    }

    public static String allToString() {
        StringBuilder sb = new StringBuilder();
        for (Operation type : Operation.values()) {
            sb.append(type.toString()).append(',');
        }
        return sb.toString();
    }
}
