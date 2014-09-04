package com.binatechnologies.seqalto.varsim;

//--- Java imports ---

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;


/**
 * Reads a DGV database flat file... not really sure if this format is stable
 */
public class DGVparser extends variantFileParser {

    private SimpleReference _reference;

    /**
     * This does not read the file, it just initialises the reading
     *
     * @param fileName  DGV flat file filename
     * @param reference Reference genome
     */
    public DGVparser(String fileName, SimpleReference reference) {
        try {
            _br = new BufferedReader(new InputStreamReader(decompressStream(fileName)));
            readLine(); // skip the first line
            readLine();
        } catch (IOException ex) {
            System.err.println("Can't open file " + fileName);
            System.err.println(ex.toString());
        }

        _reference = reference;
    }


    // returns null if the line is not a variant
    public Variant parseLine() {
        String line = _line;
        readLine();

        if (line == null || line.length() == 0)
            return null;

        String[] ll = line.split("\t");

        /*
         format
         0 1 2 3 4 5 6
         variantaccession chr start end varianttype variantsubtype reference
         7 8 9
         pubmedid method platform
         10 11 12 13 14
         mergedvariants supportingvariants mergedorsample frequency samplesize
         15 16
         observedgains observedlosses
         17 18 19
         cohortdescription genes samples
        */

        // determine variant type
        Variant.Type type;
        String varianttype = ll[4];
        String variantsubtype = ll[5];

        if (varianttype.equals("CNV")) {
            if (variantsubtype.equals("Gain")) {
                type = Variant.Type.Tandem_Duplication;
            } else if (variantsubtype.equals("Loss")) {
                type = Variant.Type.Deletion;
            } else if (variantsubtype.equals("CNV")) {
                type = Variant.Type.Tandem_Duplication;
            } else if (variantsubtype.equals("Duplication")) {
                type = Variant.Type.Tandem_Duplication;
            } else if (variantsubtype.equals("Insertion")) {
                type = Variant.Type.Insertion;
            } else if (variantsubtype.equals("Deletion")) {
                type = Variant.Type.Deletion;
            } else {
                // System.err.println("Unknown CNV Variant");
                return null;
            }
        } else if (varianttype.equals("OTHER")) {
            if (variantsubtype.equals("Tandem Duplication")) {
                type = Variant.Type.Tandem_Duplication;
            } else if (variantsubtype.equals("Inversion")) {
                type = Variant.Type.Inversion;
            } else {
                // System.err.println("Unknown OTHER Variant");
                return null;
            }
        } else {
            // System.err.println("Unknown Variant");
            return null;
        }

        // TODO right now we treat all gains as tandem duplications
        int observedgains = 0;
        if (ll[15].length() > 0) {
            observedgains = Integer.parseInt(ll[15]);
        }

        // TODO right now we treat any number of losses as a complete loss (ie.
        // deletion)
        // int observedlosses = Integer.parseInt(ll[16]);

        String var_id = ll[0];
        String chr_name = ll[1];
        int chr_idx = getChromIndex(chr_name);

        if (chr_idx <= 0) {
            return null;
        }

        int start_loc = Integer.parseInt(ll[2]);
        int end_loc = Integer.parseInt(ll[3]);

        String REF;
        FlexSeq[] alts = new FlexSeq[1];

        switch (type) {
            case Deletion:
                // reference sequence is the deletion
                // reference sequence always includes an extra character...
                byte[] temp = _reference.byteRange(chr_idx, start_loc, end_loc);
                if (temp != null) {
                    REF = new String(temp);
                } else {
                    System.err.println("Error: Invalid range");
                    System.err.println(" " + line);
                    return null;
                }
                alts[0] = new FlexSeq();
                break;
            case Insertion:
                REF = "";
                //System.err.println("Insertion: " +start_loc+ "," +end_loc+ "," +
                //(end_loc-start_loc+1));
                // TODO this is suspect... should sample from distribution of deletions.. maybe ok for now
                alts[0] = new FlexSeq(FlexSeq.Type.INS, end_loc - start_loc + 1);
                break;
            case Tandem_Duplication:
                REF = "";
                if (observedgains < 2) {
                    observedgains = 2;
                }
                //System.err.println("DUP: " +start_loc+ "," +end_loc+ "," +
                //(end_loc-start_loc+1));

                alts[0] = new FlexSeq(FlexSeq.Type.DUP, end_loc - start_loc + 1,
                        observedgains);
                break;
            case Inversion:
                REF = "";
                //System.err.println("INV: " +start_loc+ "," +end_loc+ "," +
                //(end_loc-start_loc+1));
                alts[0] = new FlexSeq(FlexSeq.Type.INV, end_loc - start_loc + 1);
                break;
            default:
                return null;
        }

        // Upper casing
        REF = REF.toUpperCase();

        byte[] refs = new byte[REF.length()];
        for (int i = 0; i < REF.length(); i++) {
            refs[i] = (byte) REF.charAt(i);
        }

        // Check

        if (REF.length() == 1 && alts[0].length() == 1) {
            // SNP
            return null; // TODO ?
        } else if (REF.length() == 0 && alts[0].length() == 0) {
            System.err.println("Skipping invalid record:");
            System.err.println(" " + line);
            return null;
        }

        byte[] phase = {1, 1};
        return new Variant(chr_name, chr_idx, start_loc, refs.length, refs,
                alts, phase, false, var_id, "PASS", String.valueOf((char) _reference
                .byteAt(chr_idx, start_loc - 1))
        );
    }

}
