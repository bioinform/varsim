package com.bina.varsim.fastqLiftover.types;

import java.util.ArrayList;

public class MafRecord {
    private final String header_;
    private final ArrayList<MafEntry> entries_;
    public MafRecord(String h, ArrayList<String> as) {
        header_ = new String(h);
        entries_ = new ArrayList<MafEntry>(as.size());
        for (int ii = 0; ii < as.size(); ++ii) {
            entries_.add(new MafEntry(as.get(ii)));
        }
    }

    public String header() {
        return header_;
    }

    public int size() {
        return entries_.size();
    }

    public MafEntry get(int ii) {
        return entries_.get(ii);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(header() + '\n');
        for (MafEntry entry : entries_) {
            sb.append(entry.src + " "
                    + Integer.toString(entry.start0) + " "
                    + Integer.toString(entry.size) + " "
                    + (entry.strand ? "+" : "-") + " "
                    + Integer.toString(entry.srcSize) + " "
                    + Integer.toString(entry.srcSize) + "\n");
        }
        return sb.toString();
    }

    static public class MafEntry {
        //https://cgwb.nci.nih.gov/FAQ/FAQformat.html#format5
        public final String src;
        public final int start0; //0-based
        public final int size;
        public final boolean strand;
        public final int srcSize;
        public final String text;

        public MafEntry(String entry) {
            String[] fields = entry.split("\\s+");
            if (fields.length != 7) throw new RuntimeException("Unexpected MAF line: " + entry);
            if (!fields[0].equals("s")) throw new RuntimeException("Unexpected MAF line: " + entry);
            src = fields[1];
            start0 = Integer.parseInt(fields[2]);
            size = Integer.parseInt(fields[3]);
            strand = fields[4].equals("+"); // true if "+", false if "-"
            if (!strand && !fields[4].equals("-")) throw new RuntimeException("Unexpected MAF line: " + entry);
            srcSize = Integer.parseInt(fields[5]);
            text = fields[6];
        }
    }
}
