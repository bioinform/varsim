package com.bina.varsim.types.variant;

/**
 * Created by johnmu on 5/20/15.
 */ // Type for whole variant
public enum VariantOverallType implements INonReference{
    Reference(false), SNP(false), Deletion(false), Insertion(true), Inversion(false), Tandem_Duplication(false), Complex(false),
    Translocation_Duplication(false), Translocation_Deletion(false), Interspersed_Duplication(false);

    boolean nonReference;

    VariantOverallType(boolean nonReference){
        this.nonReference = nonReference;
    }

    @Override
    public boolean isNonReference() {
        return nonReference;
    }

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
            case "ISP_DUP":
                return Interspersed_Duplication;
            case "TD_DUP":
                return Tandem_Duplication;
            case "TRA_DUP":
                return Translocation_Duplication;
            case "TRA_DEL":
                return Translocation_Deletion;
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
