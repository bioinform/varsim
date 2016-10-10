package com.bina.varsim.types;

/**
 * Created by guoy28 on 10/4/16.
 */
/**
 * class for storing host genome and
 * reference genome index
 */
public class HostRefIdx {
    public int hostIdx;
    public int refIdx;
    public HostRefIdx() {
        this.hostIdx = 0;
        this.refIdx = 0;
    }
    /**
     * adjust indeces(positions) of host genome and reference genome
     * based on variant type
     * @param currentMapRecord current record in map file
     */
    public void adjust_idx(MapRecord currentMapRecord) {
        switch (currentMapRecord.feature) {
            case "SEQ":
                this.hostIdx += currentMapRecord.len;
                this.refIdx += currentMapRecord.len;
                break;
            case "DEL":
                this.refIdx += currentMapRecord.len;
                break;
            case "INS":
                this.hostIdx += currentMapRecord.len;
                break;
            case "DUP_TANDEM":
                this.hostIdx += currentMapRecord.len;
                break;
            case "INV":
                this.hostIdx += currentMapRecord.len;
                break;
            case "TRANSLOCATION":
                this.hostIdx += currentMapRecord.len;
        }
    }
}
