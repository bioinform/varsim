package com.bina.varsim.types;

import java.rmi.UnexpectedException;
import java.util.*;
import java.util.concurrent.Exchanger;

/**
 * Created by guoy28 on 10/5/16.
 * CLass for VCF INFO field.
 * TODO: this class should be augmented to handle FORMAT field
 * TODO: constructors of this class should honor specifications outlined in VCF header
 */
public class VCFInfo {
    private Map<String, VCFInfoElement> info2Value;

    /**
     * parse INFO field string
     * store each key value pair in a map
     * @param infoString
     */
    public VCFInfo(String infoString) throws UnexpectedException {
        this.info2Value = new HashMap<String, VCFInfoElement>();
        String[] infos = infoString.split(";");
        for (int i = 0; i < infos.length; i++) {
            String[] keyAndValue = infos[i].split("=");
            if (keyAndValue.length > 1) {
                this.info2Value.put(keyAndValue[0], new VCFInfoElement(keyAndValue[0], keyAndValue[1]));
            } else {
                //must be boolean or flag
                this.info2Value.put(keyAndValue[0], new VCFInfoElement());
            }
        }
    }

    public Object getValue(String id) {
        return this.info2Value.containsKey(id) ? this.info2Value.get(id).getValue() : null;
    }
    private class VCFInfoElement {
        private String[] stringFields;
        //TODO: Integer should be changed to Long if varsim is used for large genome.
        private int[] numberFields;
        private Boolean flagValue;
        private String type;

        /**
         * parse comma separated value and store it
         * in proper types
         * @param id
         * @param value
         */
        public VCFInfoElement(String id, String value) throws UnexpectedException {
            this.type = getType(id);
            String[] valueArray = value.split(",");
            switch(this.type) {
                case "Integer":
                    numberFields = new int[valueArray.length];
                    for (int i = 0; i < valueArray.length; i++) {
                        numberFields[i] = Integer.parseInt(valueArray[i]);
                    }
                    break;
                case "String":
                    stringFields = valueArray;
                    break;
                default:
                    throw new UnexpectedException("ERROR: only Integer and String supported for INFO field (" + id + ").");
            }
        }

        /**
         * store id as a boolean field
         */
        public VCFInfoElement() {
            this.type = "Boolean";
            this.flagValue = true;
        }

        /**
         * return appropriate values based on types
         * @return return should be casted
         */
        public Object getValue() {
            switch (this.type) {
                case "Integer":
                    return this.numberFields;
                case "String":
                    return this.stringFields;
                case "Boolean":
                    return this.flagValue;
                default:
                    return null;
            }
        }
    }

    /**
     * return hard-coded type for some INFO IDs (including some reserved IDs in VCF
     * spec)
     * TODO: replace hard-coded infoID-type mapping with VCF header defined mapping
     *
     "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Length of variant\">\n" +
     "##INFO=<ID=POS2,Number=A,Type=Integer,Description=\"1-based Start position of source sequence\">\n" +
     "##INFO=<ID=END2,Number=A,Type=Integer,Description=\"1-based End position of source sequence\">\n" +
     "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n" +
     "##INFO=<ID=CHR2,Number=A,Type=String,Description=\"Chromosome of source sequence\">\n" +
     "##INFO=<ID=TRASUBTYPE,Number=A,Type=String,Description=\"Subtype of translocation event:" +
     " source sequence deleted (SELFISHNESS); source sequence accepted (CHIVALRY).\">\n"

     * @param infoID
     * @return
     */
    public static String getType(String infoID) {
        if (infoID.equals("SVLEN") || infoID.equals("POS2") || infoID.equals("END2") || infoID.equals("END")
                || infoID.equals("DP")) {
            return "Integer";
        } else if (infoID.equals("SVTYPE") || infoID.equals("CHR2") || infoID.equals("TRASUBTYPE")) {
            return "String";
        } else {
            throw new IllegalArgumentException("ERROR: unrecognized INFO ID (" + infoID + ").");
        }
    }
}
