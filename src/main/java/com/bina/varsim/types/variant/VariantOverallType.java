package com.bina.varsim.types.variant;

/**
 * Created by johnmu on 5/20/15.
 */ // Type for whole variant
public enum VariantOverallType {
    Reference, SNP, Deletion, Insertion, Inversion, Tandem_Duplication, Complex;

    public static VariantOverallType fromString(String s){
        for (VariantOverallType type : VariantOverallType.values()) {
            if(s.equalsIgnoreCase(type.toString())){
                return type;
            }
        }

        // Shorter ones
        switch(s.toUpperCase()){
            case "DEL":
                return Deletion;
            case "INS":
                return Insertion;
            case "SNP":
                return SNP;
            case "REF":
                return Reference;
            case "INV":
                return Inversion;
            case "TD_DUP":
                return Tandem_Duplication;
            case "CLPX":
                return Complex;
            default:
                throw new RuntimeException(String.format("Invalid variant type: %s\n Valid ones are %s", s, VariantOverallType.allToString()));
        }
    }

    public static String allToString() {
        StringBuilder sb = new StringBuilder();
        for (VariantOverallType type : VariantOverallType.values()) {
            sb.append(type.toString()).append(',');
        }
        return sb.toString();
    }
}
