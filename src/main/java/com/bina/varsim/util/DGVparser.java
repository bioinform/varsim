package com.bina.varsim.util;

//--- Java imports ---


import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.variant.VariantType;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.types.variant.alt.Alt;
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

    public final static String DGV_HEADER_START = "#";

    public final static String DGV_COLUMN_SEPARATOR = "\t";

    Random rand = null;

    private SimpleReference reference;

    /**
     * This does not read the file, it just initialises the reading
     *
     * @param fileName  DGV flat file filename
     * @param reference Reference genome
     */
    public DGVparser(String fileName, SimpleReference reference, Random rand) {
        this.rand = rand;

        try {
            bufferedReader = new BufferedReader(new InputStreamReader(decompressStream(fileName)));
            readLine(); // skip the first line
            readLine();
        } catch (IOException ex) {
            log.error("Can't open file " + fileName);
            log.error(ex.toString());
        }

        this.reference = reference;
    }

    /**
     *
     * @param dgvVariantType DGV Variant type
     * @param dgvVariantSubtype DGV Variant subtype
     * @return VariantType. Return value is null if cannot map to the right VariantType
     */
    public VariantType getVariantType(final String dgvVariantSubtype) {

      /* as of 12/16/2016 the following subtypes are found in dgvdb
        complex
                deletion
        duplication
                gain
        gain+loss
        insertion
                inversion
        loss
        mobile element insertion
        novel sequence insertion
        sequence alteration
        tandem duplication

        subtypes not handled by dgvparser will be ignored
        */
      switch (dgvVariantSubtype) {
        case "gain":
          return VariantType.Tandem_Duplication;
        case "loss":
          return VariantType.Deletion;
        case "cnv":
          return VariantType.Tandem_Duplication;
        case "duplication":
          return VariantType.Tandem_Duplication;
        case "insertion":
        case "mobile element insertion":
        case "novel sequence insertion":
          return VariantType.Insertion;
        case "deletion":
          return VariantType.Deletion;
        case "tandem duplication":
          return VariantType.Tandem_Duplication;
        case "inversion":
          return VariantType.Inversion;
        default:
          throw new IllegalArgumentException("Unrecognized subtype: " + dgvVariantSubtype);
      }
    }

    /**
     *
     * @param line
     * @return True if the line is a header-line, false otherwise
     */
    public boolean isHeaderLine(final String line) {
        return line.startsWith(DGV_HEADER_START);
    }
    
    // returns null if the line is not a variant
    public Variant parseLine() {
      String line = this.line;
      try {
        readLine();

        if (line == null || line.length() == 0) {
          return null;
        }

        if (isHeaderLine(line)) {
          return null;
        }

        String[] ll = line.split(DGV_COLUMN_SEPARATOR);

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
        final String var_id = ll[0];
        final ChrString chr = new ChrString(ll[1]);
        final int start_loc = Integer.parseInt(ll[2]);
        final int end_loc = Integer.parseInt(ll[3]);

        // determine variant type
        VariantType type = getVariantType(ll[5].toLowerCase());

        // TODO right now we treat all gains as tandem duplications
	if (ll[15].length() == 0)
	    throw new RuntimeException("no observedgain");
        int observedgains = Integer.parseInt(ll[15]);

        String REF;
        Alt[] alts = new Alt[1];

        switch (type) {
          case Deletion:
            // reference sequence is the deletion
            // reference sequence always includes an extra character...
            byte[] temp = reference.byteRange(chr, start_loc, end_loc);
            if (temp != null) {
              REF = new String(temp);
            } else {
              log.error("Error: Invalid range");
              log.error(" " + line);
              return null;
            }
            alts[0] = new Alt(new FlexSeq());
            break;
          case Insertion:
            REF = "N";
            // TODO this is suspect... should sample from distribution of deletions.. maybe ok for now
            alts[0] = new Alt(new FlexSeq(FlexSeq.Type.INS, end_loc - start_loc + 1));
            break;
          case Tandem_Duplication:
            REF = "N";
            alts[0] = new Alt(new FlexSeq(FlexSeq.Type.TANDEM_DUP, end_loc - start_loc + 1,observedgains+1));
            break;
          case Inversion:
            REF = "N";
            alts[0] = new Alt(new FlexSeq(FlexSeq.Type.INV, end_loc - start_loc + 1));
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
        return new Variant.Builder().chr(chr).pos(start_loc).referenceAlleleLength(refs.length).
                ref(refs).alts(alts).phase(phase).isPhased(false).varId(var_id).filter("PASS").
                refDeleted(String.valueOf((char) reference.byteAt(chr, start_loc - 1))).randomNumberGenerator(rand).build();
      } catch (Exception e) {
        log.warn(e.getMessage());
        log.warn("line: " + line);
        return null;
      }
    }

}
