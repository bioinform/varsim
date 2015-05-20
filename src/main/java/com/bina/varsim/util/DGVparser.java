package com.bina.varsim.util;

//--- Java imports ---


import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.VariantType;
import com.bina.varsim.types.Variant;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Random;


/**
 * Reads a DGV database flat file... not really sure if this format is stable
 */
public class DGVparser extends GzFileParser<Variant> {
    private final static Logger log = Logger.getLogger(DGVparser.class.getName());

    Random _rand = null;

    private SimpleReference _reference;

    /**
     * This does not read the file, it just initialises the reading
     *
     * @param fileName  DGV flat file filename
     * @param reference Reference genome
     */
    public DGVparser(String fileName, SimpleReference reference, Random rand) {
        _rand = rand;

        try {
            _br = new BufferedReader(new InputStreamReader(decompressStream(fileName)));
            readLine(); // skip the first line
            readLine();
        } catch (IOException ex) {
            log.error("Can't open file " + fileName);
            log.error(ex.toString());
        }

        _reference = reference;
    }


    // returns null if the line is not a variant
    public Variant parseLine() {
        String line = _line;
        readLine();

        if (line == null || line.length() == 0) {
            return null;
        }

        if (line.charAt(0) == '#') {
            return null;
        }

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
        VariantType type;
        String varianttype = ll[4];
        String variantsubtype = ll[5];

        switch (varianttype) {
            case "CNV":
                switch (variantsubtype) {
                    case "Gain":
                        type = VariantType.Tandem_Duplication;
                        break;
                    case "Loss":
                        type = VariantType.Deletion;
                        break;
                    case "CNV":
                        type = VariantType.Tandem_Duplication;
                        break;
                    case "Duplication":
                        type = VariantType.Tandem_Duplication;
                        break;
                    case "Insertion":
                        type = VariantType.Insertion;
                        break;
                    case "Deletion":
                        type = VariantType.Deletion;
                        break;
                    default:
                        return null;
                }
                break;
            case "OTHER":
                switch (variantsubtype) {
                    case "Tandem Duplication":
                        type = VariantType.Tandem_Duplication;
                        break;
                    case "Inversion":
                        type = VariantType.Inversion;
                        break;
                    default:
                        return null;
                }
                break;
            default:
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
        ChrString chr = new ChrString(ll[1]);

        int start_loc = Integer.parseInt(ll[2]);
        int end_loc = Integer.parseInt(ll[3]);

        String REF;
        FlexSeq[] alts = new FlexSeq[1];

        switch (type) {
            case Deletion:
                // reference sequence is the deletion
                // reference sequence always includes an extra character...
                byte[] temp = _reference.byteRange(chr, start_loc, end_loc);
                if (temp != null) {
                    REF = new String(temp);
                } else {
                    log.error("Error: Invalid range");
                    log.error(" " + line);
                    return null;
                }
                alts[0] = new FlexSeq();
                break;
            case Insertion:
                REF = "";
                // TODO this is suspect... should sample from distribution of deletions.. maybe ok for now
                alts[0] = new FlexSeq(FlexSeq.Type.INS, end_loc - start_loc + 1);
                break;
            case Tandem_Duplication:
                REF = "";
                if (observedgains < 2) {
                    observedgains = 2;
                }
                alts[0] = new FlexSeq(FlexSeq.Type.DUP, end_loc - start_loc + 1,
                        observedgains);
                break;
            case Inversion:
                REF = "";
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
            log.error("Skipping invalid record:");
            log.error(" " + line);
            return null;
        }

        byte[] phase = {1, 1};
        return new Variant(chr, start_loc, refs.length, refs,
                alts, phase, false, var_id, "PASS", String.valueOf((char) _reference
                .byteAt(chr, start_loc - 1)), _rand);
    }

}
