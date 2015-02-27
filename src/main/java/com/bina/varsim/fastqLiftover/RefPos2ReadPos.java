package com.bina.varsim.fastqLiftover;

import java.lang.Math;

public class RefPos2ReadPos {
    /**
     * get the base-b index of the read corresponding to the base-b index of the reference
     *
     * @param pos base-b index of reference
     * @param b base of pos, typically 0 or 1
     * @return the base-b index of read sequence
     */
    public int get(int pos, int b) {
        final int index = pos - b - shift0_;
        if(index<0) {
            throw new RuntimeException(  Integer.toString(pos) + "< ["
                                       + Integer.toString(shift0_+b) + ","
                                       + Integer.toString(shift0_+map0_.length+b) + ")");
        }
        if(index>=map0_.length) {
            throw new RuntimeException(  "[" + Integer.toString(shift0_+b)
                                       + "," + Integer.toString(shift0_+map0_.length+b)
                                       + ") < " + Integer.toString(pos) );
        }
        return map0_[index] + b;
    }

    public RefPos2ReadPos(final MafRecord.MafEntry ref, final MafRecord.MafEntry seq) {
        final int alignment_length = ref.text.length();
        if( alignment_length != seq.text.length() ) throw new RuntimeException("corrupted alignment strings");
        map0_ = new int[ref.size];
        shift0_ = ref.start0;

        int r_pos = 0;
        int s_count = 0;
        for (int ii = 0, s_pos = s_count; ii < alignment_length ; ++ii) {
            if( ref.text.charAt(ii) != '-' ) map0_[r_pos++] = s_pos;
            if( seq.text.charAt(ii) != '-' ) s_pos = Math.min(++s_count,seq.size-1);
        }
        if (r_pos != ref.size) throw new RuntimeException("not enough ref positions in MAF alignment");
        if (s_count > seq.size) throw new RuntimeException("too many seq positions in MAF alignment");
    }
    final private int shift0_; //overall shift to key index
    final private int[] map0_; //0-base to 0-base mapping
}
