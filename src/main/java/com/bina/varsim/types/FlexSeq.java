package com.bina.varsim.types;

import java.io.UnsupportedEncodingException;
import java.util.Arrays;

/**
 * This can represent a byte[] or a more flexible sequence such as insertion with unknown contents, dup, inv, etc...
 *
 * @author johnmu
 */

public class FlexSeq {
    Type type;
    int length;
    int copyNumber;
    byte[] sequence;
    String variantId = "";
    ChrString chr2;
    int pos2;
    int end2;
    int referenceAlleleLength;

    /**
     * Empty FlexSeq, default to a SEQ type
     */
    public FlexSeq() {
        // empty sequence
        type = Type.SEQ;
        length = 0;
        copyNumber = 1;
        sequence = new byte[0];
    }

    /**
     * Single byte FlexSeq, eg. for a SNP
     *
     * @param seq single byte
     */
    public FlexSeq(final byte seq) {
        // regular byte array
        sequence = new byte[1];
        sequence[0] = seq;
        type = Type.SEQ;
        length = 1;
        copyNumber = 1;
    }

    /**
     * Byte array FlexSeq, this for an indel like thing
     *
     * @param seq byte array of sequence to be stored
     */
    public FlexSeq(final byte[] seq) {
        // regular byte array
        sequence = seq.clone();
        type = Type.SEQ;
        length = seq.length;
        copyNumber = 1;
    }

    /**
     * FlexSeq with unknown sequence but known type, eg inversion
     *
     * @param t   Type of the sequence
     * @param len Length of the sequence
     */
    public FlexSeq(final Type t, int len) {
        if (len < 0) {
            len = 0;
        }
        // regular byte array
        sequence = null;
        type = t;
        length = len;
        copyNumber = 1;
    }

    /**
     * FlexSeq with unknown sequence but known type and copy number, eg duplication
     *
     * @param t        Type of the sequence
     * @param len      Length of the sequence
     * @param copy_num Copy number of sequence
     */
    public FlexSeq(final Type t, final int len, final int copy_num) {
        // regular byte array
        this(t, len);
        copyNumber = copy_num;
    }


    public FlexSeq(final FlexSeq b) {
        type = b.type;
        length = b.length;
        copyNumber = b.copyNumber;
        if (b.sequence == null) {
            sequence = null;
        } else {
            sequence = b.sequence.clone();
        }
    }

    /**
     * Return byte at some index in the sequence
     *
     * @param i index in sequence
     * @return
     */
    public byte charAt(final int i) {
        if (sequence != null) {
            return sequence[i];
        } else {
            return 0;
        }
    }

    /**
     * @param beginIndex start of substring
     * @return array from beginIndex to end
     */
    public byte[] substring(final int beginIndex) {
        return Arrays.copyOfRange(sequence, beginIndex, sequence.length);
    }

    /**
     * @param beginIndex start of substring (inclusive)
     * @param endIndex   end of substring (exclusive)
     * @return substring from [beginIndex,endIndex)
     */
    public byte[] substring(final int beginIndex, final int endIndex) {
        return Arrays.copyOfRange(sequence, beginIndex, endIndex);
    }

    /**
     * @return length of the sequence
     */
    public int length() {
        return length;
    }

    /**
     * @return length of the sequence accounting for copy numbers
     */
    public int var_length() {
        if (type == Type.DUP) {
            return length * copyNumber;
        }
        return length;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FlexSeq flexSeq = (FlexSeq) o;

        if (copyNumber != flexSeq.copyNumber) return false;
        if (length != flexSeq.length) return false;
        if (!Arrays.equals(sequence, flexSeq.sequence)) return false;
        if (type != flexSeq.type) return false;
        if (variantId != null ? !variantId.equals(flexSeq.variantId) : flexSeq.variantId != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = type != null ? type.hashCode() : 0;
        result = 31 * result + length;
        result = 31 * result + copyNumber;
        result = 31 * result + (sequence != null ? Arrays.hashCode(sequence) : 0);
        result = 31 * result + (variantId != null ? variantId.hashCode() : 0);
        return result;
    }

    public int getCopy_num() {
        return copyNumber;
    }

    public void setCopy_num(int cn) {
        copyNumber = cn;
    }

    public Type getType() {
        return type;
    }

    public void setType(final Type t) {
        type = t;
    }

    public String getVar_id() {
        return variantId;
    }

    public void setVar_id(final String var_id) {
        variantId = var_id;
    }

    public void setLength(int l) {
        length = l;
    }

    public byte[] getSeq() {
        return sequence;
    }

    // has the actual sequence
    public boolean isSeq() {
        return (type == Type.SEQ);
    }

    public String toString() {
        if (sequence == null) {
            switch (type) {
                case DUP:
                    return "<DUP:TANDEM>";
                case INS:
                    return "<INS>";
                case INV:
                    return "<INV>";
                case DEL:
                    return "<DEL>";
                case TRA:
                    return "<TRA>";
                default:
                    break;
            }
        } else {
            try {
                return new String(sequence, "US-ASCII");
            } catch (UnsupportedEncodingException e) {
                e.printStackTrace();
            }
        }
        return "";
    }

    public boolean equals(final FlexSeq seq) {
        if (type != seq.type) {
            return false;
        }
        if (length != seq.length) {
            return false;
        }
        if (copyNumber != seq.copyNumber) {
            return false;
        }

        if (!Arrays.equals(sequence, seq.sequence)) {
            return false;
        }

        return true;
    }

    /**
     * SEQ means the byte sequence is given
     * INV = inversion
     * Deletion = delelton (unknown reference )
     * DUP = duplication
     * Insertion = insertion (unknown sequence)
     * TRA = translocation (unknown sequence, sequence is defined at variant level)
     */
    public enum Type {
        SEQ, INV, DUP, INS, DEL, TRA
    }
    public void setChr2(ChrString chr2) {
        this.chr2 = chr2;
    }
    public void setPos2(int pos2) {
        this.pos2 = pos2;
    }

    public ChrString getChr2() {
        return chr2;
    }

    public int getPos2() {
        return pos2;
    }

    public int getEnd2() {
        return end2;
    }

    public void setEnd2(final int end2) {
        this.end2 = end2;
    }

    public void setReferenceAlleleLength(int referenceAlleleLength) {
        this.referenceAlleleLength = referenceAlleleLength;
    }

    public int getReferenceAlleleLength() {

        return referenceAlleleLength;
    }
}
