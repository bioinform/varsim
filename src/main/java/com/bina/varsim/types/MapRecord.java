package com.bina.varsim.types;

import java.util.Hashtable;
import java.util.StringJoiner;

import static com.bina.varsim.tools.simulation.VCF2diploid.DELETED_BASE;

/**
 * Created by guoy28 on 10/4/16.
"#Len\tHOST_chr\tHOST_pos\tREF_chr\tREF_pos\tDIRECTION\tFEATURE\tVAR_ID"
 this is more like a struct :)
 I think MapRecord stands for records in MFF file (the map file)

 TODO: MapRecord contains essentially same information as Variant, these two classes can be consolidated.
 */
public class MapRecord {
    //len for length of alternative allele
    public int len = 0;
    public String hostChr = "";
    public int hostPos = 0;
    public String refChr = "";
    public int refPos = 0;
    public boolean isForward = true; // true is forward
    public String feature = ".";
    public String varId = ".";

    private MapRecord() {}
    /**
     * create a new map file record based on the variant or sequence object stored at insertPosition2Sequence[idx]
     * from the perspective of good programming practice, it's better to iterate over all variants
     * in insertPosition2Sequence rather than using idx.
     * after the record is created, advance host genome and reference indices so we know where next
     * event occurs. For some variants, e.g. TANDEM_DUP, several records will be created on the fly
     * and got appended to output string sb.
     * TODO: separate creation of new map file records with appending to output string (for simplification)
     *
     * @param sb output string
     * @param idx position of examination
     * @param chr_name
     * @param ref_chr_name
     * @param hostRefIdx
     * @param genome
     * @param insertPosition2Sequence
     * @return
     */
    public static MapRecord generateNewMapRecord(StringBuilder sb, int idx, String chr_name, String ref_chr_name, HostRefIdx hostRefIdx,
                      byte[] genome, Hashtable<Integer, FlexSeq> insertPosition2Sequence) {
        MapRecord currentMapRecord = new MapRecord();
        currentMapRecord.hostChr = chr_name;
        currentMapRecord.refChr = ref_chr_name;

        // compute what the first line is and store as object

        // if it is inserted, we copy the varId to the getReferenceAlleleLength
        boolean inserted = false;
        String varId = ".";
        if (insertPosition2Sequence.containsKey(idx + 1)) {
            inserted = true;
            // insertion at currentMapRecord location
            FlexSeq insertion = insertPosition2Sequence.get(idx + 1);
            varId = insertion.getVar_id();

            //System.err.println("Check type: " + ins.getType());

            // need to check what kind of insertion it is
            switch (insertion.getType()) {
                case SEQ:
                case INS:
                    currentMapRecord.hostPos = hostRefIdx.hostIdx;
                    currentMapRecord.refPos = hostRefIdx.refIdx - 1;
                    currentMapRecord.feature = "INS";
                    currentMapRecord.isForward = true;
                    currentMapRecord.len = insertion.var_length();
                    currentMapRecord.varId = varId;

                    break;
                case TRANSLOCATION:
                    currentMapRecord.hostPos = hostRefIdx.hostIdx;
                    currentMapRecord.refPos = Math.min(insertion.getPos2(), insertion.getEnd2()) - 1;
                    currentMapRecord.refChr = insertion.getChr2().toString();
                    currentMapRecord.feature = "TRANSLOCATION";
                    currentMapRecord.isForward = insertion.getPos2() <= insertion.getEnd2() ? true : false;
                    currentMapRecord.len = insertion.var_length();
                    currentMapRecord.varId = varId;
                    break;
                case INV:
                    // TODO: treat inversion like MNP
                    currentMapRecord.hostPos = hostRefIdx.hostIdx;
                    /*position of inversion on reference genome is questionable
                    inversion does not change length of reference, so is it
                    still necessary to set its position on reference as 1bp
                    before the event?
                    */
                    currentMapRecord.refPos = hostRefIdx.refIdx - 1;
                    currentMapRecord.feature = "INV";
                    //why direction is false (negative strand)?
                    currentMapRecord.isForward = false;
                    currentMapRecord.len = insertion.var_length();
                    currentMapRecord.varId = varId;

                    break;
                case DUP:
                    // need to replicate several blocks
                    int cn = insertion.getCopy_num();

                    // first build one
                    currentMapRecord.hostPos = hostRefIdx.hostIdx;
                    currentMapRecord.refPos = hostRefIdx.refIdx - 1;
                    currentMapRecord.feature = "DUP_TANDEM";
                    currentMapRecord.isForward = true;
                    currentMapRecord.len = insertion.length();
                    currentMapRecord.varId = varId;

                    // iterate
                    for (int i = 1; i < cn; i++) {
                        hostRefIdx.adjust_idx(currentMapRecord);
                        sb.append(currentMapRecord);
                        sb.append('\n');
                        currentMapRecord = new MapRecord();
                        currentMapRecord.hostChr = chr_name;
                        currentMapRecord.refChr = ref_chr_name;
                        currentMapRecord.hostPos = hostRefIdx.hostIdx;
                        currentMapRecord.refPos = hostRefIdx.refIdx - 1;
                        currentMapRecord.feature = "DUP_TANDEM";
                        currentMapRecord.isForward = true;
                        currentMapRecord.len = insertion.length();
                        currentMapRecord.varId = varId;
                    }

                    break;
            }

            // output it
            hostRefIdx.adjust_idx(currentMapRecord);
            sb.append(currentMapRecord);
            sb.append('\n');

            currentMapRecord = new MapRecord();
            currentMapRecord.hostChr = chr_name;
            currentMapRecord.refChr = ref_chr_name;
        }

        if (genome[idx] == DELETED_BASE) {
            // deleted base
            currentMapRecord.hostPos = hostRefIdx.hostIdx - 1;
            currentMapRecord.refPos = hostRefIdx.refIdx;
            currentMapRecord.feature = "DEL";
            currentMapRecord.isForward = true;
            /*
            although length is 1 here, currentMapRecord is just the beginning
            of a block, length may increase later (after current
            method returns).
             */
            currentMapRecord.len = 1;
            if (inserted) {
                currentMapRecord.varId = varId;
            }
        } else {
            // regular sequence
            currentMapRecord.hostPos = hostRefIdx.hostIdx;
            currentMapRecord.refPos = hostRefIdx.refIdx;
            currentMapRecord.feature = "SEQ";
            currentMapRecord.isForward = true;
            /*
            although length is 1 here, currentMapRecord is just the beginning
            of a block, length may increase later (after current
            method returns).
             */
            currentMapRecord.len = 1;
        }
        return currentMapRecord;
    }

    public String toString() {
        StringJoiner joiner = new StringJoiner("\t");
        joiner.add(Integer.toString(len));
        joiner.add(hostChr);
        joiner.add(Integer.toString(hostPos));
        joiner.add(refChr);
        joiner.add(Integer.toString(refPos));
        if (isForward) {
            joiner.add("+");
        } else {
            joiner.add("-");
        }
        joiner.add(feature);
        joiner.add(varId);

        return joiner.toString();
    }
}

