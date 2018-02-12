package com.bina.varsim.types;

import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This can represent a byte[] or a more flexible sequence such as insertion with unknown contents, dup, inv, etc...
 * immutable class
 *
 * @author johnmu
 */

public final class FlexSeq {
    //TODO make all fields final if possible
    public static Set<Byte> LEGAL_SEQUENCE_CHARS =
            Stream.of((byte)'A',(byte)'T',(byte)'C',(byte)'G',
                    (byte)'a',(byte)'t',(byte)'c',(byte)'g',
                    (byte)'N',(byte)'n').collect(Collectors.toCollection(HashSet<Byte>::new));
    Type type;
    int length;
    int copyNumber = 1;
    byte[] sequence;
    String variantId = "";
    ChrString chr2;
    int pos2;
    int end2;
    int referenceAlleleLength;
    boolean isinv;
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
        boolean isinv;

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
        public Builder isinv(final boolean b) {
            this.isinv = b;
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
        checkSequence(sequence);
        this.variantId = b.variantId;
        this.chr2 = b.chr2;
        this.pos2 = b.pos2;
        this.end2 = b.end2;
        this.referenceAlleleLength = b.referenceAlleleLength;
        this.isinv = b.isinv;
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
        checkSequence(sequence);
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
        checkSequence(sequence);
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
        /*temporarily lift this restriction*/
        /*if (len < 0) {
            len = 0;
        }*/
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
        isinv = b.isinv;
        if (b.sequence == null) {
            sequence = null;
        } else {
            sequence = b.sequence.clone();
        }
        checkSequence(sequence);
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

    public boolean isInversed() {
        return isinv;
    }

    /**
     * @return length of the sequence accounting for copy numbers
     */
    public int varLength() {
        if (type == Type.TANDEM_DUP || type == Type.ISP_DUP || type == Type.TRA_DUP) {
            return length * copyNumber;
        }
        return length;
    }

    /**
     * check if sequence contains illegal characters
     * only ATCGatcgNn are allowed
     * @param s
     * @return
     */
    private void checkSequence(byte[] s) {
       if (s != null) {
           for (byte b : s) {
               if (!LEGAL_SEQUENCE_CHARS.contains(b)) {
                   throw new IllegalArgumentException("Found illegal character: " + b);
               }
           }
       }
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
                case TANDEM_DUP:
                    return "<DUP:TANDEM>";
                case TRA_DUP:
                    return "<DUP:TRA>";
                case ISP_DUP:
                    return "<DUP:ISP>";
                case INS:
                    return "<INS>";
                case INV:
                    return "<INV>";
                case DEL:
                    return "<DEL>";
                case TRA_DEL:
                    return "<DEL:TRA>";
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
     * INV = simple inversion
     * DEL = delelton (unknown reference )
     * TRA_DEL = delelton in translocation
     * TANDEM_DUP = tandem duplication
     * INS = insertion (unknown sequence)
     * ISP_DUP = interspersed duplication
     * TRA_DUP = translocation duplication (essentially interspersed duplication)
     */
    public enum Type {
        SEQ, INV, TANDEM_DUP, TRA_DUP, ISP_DUP, INS, DEL, TRA_DEL
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
