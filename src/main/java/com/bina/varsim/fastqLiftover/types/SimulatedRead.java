package com.bina.varsim.fastqLiftover.types;

import com.google.common.base.Joiner;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

public class SimulatedRead {
    private final static Logger log = Logger.getLogger(SimulatedRead.class.getName());
    private static final Joiner joiner = Joiner.on(",").skipNulls();
    private static final int MAX_SEQ_SELECT = 2;
    public List<GenomeLocation> locs1; // List of locations of the alignment of the first end
    public List<GenomeLocation> locs2; // List of locations of the alignment of the second end
    public List<GenomeLocation> origLocs1; // Original alignment location as reported by the read simulator for the first end
    public List<GenomeLocation> origLocs2; // Original alignment location as reported by the read simulator for the second end
    public int random1 = 0;
    public int random2 = 0;
    public int seqErrors1 = 0;
    public int snps1 = 0;
    public int indels1 = 0;
    public int seqErrors2 = 0;
    public int snps2 = 0;
    public int indels2 = 0;
    public String readId = "";
    public int fragment = 1;
    public int laneId = 0;
    public String sequence;
    public String quality;
    public int alignedBases1 = 0;
    public int alignedBases2 = 0;

    public SimulatedRead() {
        locs1 = new ArrayList<>();
        locs2 = new ArrayList<>();
        origLocs1 = new ArrayList<>();
        origLocs2 = new ArrayList<>();
    }

    /* Parses the VarSim generated read-name and fills the fields, except for the quality and sequence */
    public SimulatedRead(final String readName) {
        final String fields[] = readName.split("[:/]", -1);
        random1 = decodeInt(fields[4]);
        random2 = decodeInt(fields[5]);
        seqErrors1 = decodeInt(fields[6]);
        snps1 = decodeInt(fields[7]);
        indels1 = decodeInt(fields[8]);
        seqErrors2 = decodeInt(fields[9]);
        snps2 = decodeInt(fields[10]);
        indels2 = decodeInt(fields[11]);
        readId = fields[12];
        laneId = Integer.parseInt(fields[13]);

        locs1 = parseLocations(fields[0]);
        locs2 = parseLocations(fields[1]);
        origLocs1 = parseLocations(fields[2]);
        origLocs2 = parseLocations(fields[3]);
    }

    public List<GenomeLocation> parseLocations(final String locationsString) {
        final String fields[] = locationsString.split(",", -1);
        List<GenomeLocation> locations = new ArrayList<>();
        for (String field : fields) {
            if (!("".equals(field))) {
                locations.add(new GenomeLocation(field));
            }
        }
        return locations;
    }

    public int getLength() {
        return sequence.length();
    }

    public String encodeInt(int n) {
        return (n == 0) ? "" : Integer.toString(n);
    }

    public int decodeInt(String n) {
        return "".equals(n) ? 0 : Integer.parseInt(n);
    }

    public String encodeInts(int... numbers) {
        StringBuilder sb = new StringBuilder();
        for (int number : numbers) {
            sb.append(":");
            if (number != 0) {
                sb.append(number);
            }
        }
        return sb.toString();
    }

    public String getName() {
        return joiner.join(locs1) + ":" + joiner.join(locs2) + ":" + joiner.join(origLocs1) + ":" + joiner.join(origLocs2)
                + encodeInts(random1, random2, seqErrors1, snps1, indels1, seqErrors2, snps2, indels2) + ":" + readId + ":" + laneId + "/" + fragment;

    }

    private List<GenomeLocation> selectLocs(final List<GenomeLocation> locs) {
        List<GenomeLocation> selected = new ArrayList<>();
        int seqSelected = 0;
        for (GenomeLocation g : locs) {
            if ("S".equals(g.feature) || "".equals(g.feature)) {
                if (seqSelected < MAX_SEQ_SELECT) {
                    selected.add(g);
                    seqSelected++;
                }
            } else {
                selected.add(g);
            }
        }
        return selected;
    }

    public String getShortName() {
        return joiner.join(selectLocs(locs1)) + ":" + joiner.join(selectLocs(locs2)) + ":" + joiner.join(origLocs1) + ":" + joiner.join(origLocs2)
                + encodeInts(random1, random2, seqErrors1, snps1, indels1, seqErrors2, snps2, indels2) + ":" + readId + ":" + laneId + "/" + fragment;

    }

    public String toString() {
        String name = getName();
        if (name.length() > 255) {
            log.warn("Read name " + name + " too long " + name.length() + " chars. Using only a few locations instead.");
            name = getShortName();
        }
        return "@" + name + "\n" + sequence + "\n+\n" + quality;
    }
}
