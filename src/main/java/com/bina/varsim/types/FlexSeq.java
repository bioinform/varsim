package com.bina.varsim.types;

import com.bina.varsim.types.variant.Variant;

import java.io.UnsupportedEncodingException;
import java.util.Arrays;

/**
 * This can represent a byte[] or a more flexible sequence such as insertion with unknown contents, dup, inv, etc...
 * immutable class
 *
 * @author johnmu
 */

public final class FlexSeq {
    //TODO make all fields final if possible
    Type type;
    int length;
    int copyNumber = 1;
    byte[] sequence;
    String variantId = "";
    ChrString chr2;
    int pos2;
    int end2;
    int referenceAlleleLength;
    volatile int hashCode;

    /*
    use Builder pattern to make code more readable and maintainable.
     */
    public static class Builder {
        Type type;
        int length;
        int copyNumber = 1;
        byte[] sequence;
        String variantId = "";
        ChrString chr2;
        int pos2;
        int end2;
        int referenceAlleleLength;

        public Builder() {}
        public Builder type(final Type t) {
            this.type = t;
            return this;
        }
        public Builder length(final int l) {
            this.length = l;
            return this;
        }
        public Builder referenceAlleleLength(final int referenceAlleleLength) {
            this.referenceAlleleLength = referenceAlleleLength;
            return this;
        }
        public Builder copyNumber(final int c) {
            this.copyNumber = c;
            return this;
        }
        public Builder sequence(byte[] s) {
            this.sequence = s.clone();
            return this;
        }
        public Builder variantId(final String i) {
            this.variantId = i;
            return this;
        }
        public Builder chr2(ChrString c) {
            this.chr2 = c;
            return this;
        }
        public Builder pos2(final int p) {
          this.pos2 = p;
            return this;
        }
        public Builder end2(final int e) {
          this.end2 = e;
            return this;
        }
        public FlexSeq build() {
            return new FlexSeq(this);
        }
    }
    private FlexSeq(Builder b) {
        this.type = b.type;
        this.length = b.length;
        this.copyNumber = b.copyNumber;
        this.sequence = b.sequence;
        this.variantId = b.variantId;
        this.chr2 = b.chr2;
        this.pos2 = b.pos2;
        this.end2 = b.end2;
        this.referenceAlleleLength = b.referenceAlleleLength;
    }
    /**
     * Empty FlexSeq, default to a SEQ type
     */
    public FlexSeq() {
        // empty sequence
        type = Type.SEQ;
        length = 0;
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
    public byte byteAt(final int i) {
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
    public int varLength() {
        if (type == Type.DUP) {
            return length * copyNumber;
        }
        return length;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof FlexSeq)) return false;

        FlexSeq flexSeq = (FlexSeq) o;

        if (copyNumber != flexSeq.copyNumber) return false;
        if (length != flexSeq.length) return false;
        if (type != flexSeq.type) return false;
        if (variantId != null ? !variantId.equals(flexSeq.variantId) : flexSeq.variantId != null) return false;
        if (!Arrays.equals(sequence, flexSeq.sequence)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = hashCode;
        if (result != 0) return result;
       /*initial hash code cannot be zero, otherwise, we cannot
       tell the difference from hashCode if the first fields
       are zero.
        */
        result = type != null ? type.hashCode() : 17;
        result = 31 * result + length;
        result = 31 * result + copyNumber;
        result = 31 * result + (sequence != null ? Arrays.hashCode(sequence) : 0);
        result = 31 * result + (variantId != null ? variantId.hashCode() : 0);
        hashCode = result;
        return result;
    }

    public int getCopyNumber() {
        return copyNumber;
    }

    public Type getType() {
        return type;
    }

    public String getVariantId() {
        return variantId;
    }

    public byte[] getSequence() {
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

    public ChrString getChr2() {
        return chr2;
    }

    public int getPos2() {
        return pos2;
    }

    public int getEnd2() {
        return end2;
    }

    public int getReferenceAlleleLength() {

        return referenceAlleleLength;
    }
}
