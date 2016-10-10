package com.bina.varsim.types;

import java.io.UnsupportedEncodingException;
import java.util.Arrays;

/**
 * This can represent a byte[] or a more flexible sequence such as insertion with unknown contents, dup, inv, etc...
 *
 * @author johnmu
 */

public class FlexSeq {
    Type _type;
    int _len;
    int _copy_num;
    byte[] _seq;
    String _var_id = "";
    ChrString chr2;
    int pos2;
    int end2;
    int referenceAlleleLength;

    /**
     * Empty FlexSeq, default to a SEQ type
     */
    public FlexSeq() {
        // empty sequence
        _type = Type.SEQ;
        _len = 0;
        _copy_num = 1;
        _seq = new byte[0];
    }

    /**
     * Single byte FlexSeq, eg. for a SNP
     *
     * @param seq single byte
     */
    public FlexSeq(byte seq) {
        // regular byte array
        _seq = new byte[1];
        _seq[0] = seq;
        _type = Type.SEQ;
        _len = 1;
        _copy_num = 1;
    }

    /**
     * Byte array FlexSeq, this for an indel like thing
     *
     * @param seq byte array of sequence to be stored
     */
    public FlexSeq(byte[] seq) {
        // regular byte array
        _seq = seq.clone();
        _type = Type.SEQ;
        _len = seq.length;
        _copy_num = 1;
    }

    /**
     * FlexSeq with unknown sequence but known type, eg inversion
     *
     * @param t   Type of the sequence
     * @param len Length of the sequence
     */
    public FlexSeq(Type t, int len) {
        if (len < 0) {
            len = 0;
        }
        // regular byte array
        _seq = null;
        _type = t;
        _len = len;
        _copy_num = 1;
    }

    /**
     * FlexSeq with unknown sequence but known type and copy number, eg duplication
     *
     * @param t        Type of the sequence
     * @param len      Length of the sequence
     * @param copy_num Copy number of sequence
     */
    public FlexSeq(Type t, int len, int copy_num) {
        // regular byte array
        this(t, len);
        _copy_num = copy_num;
    }


    public FlexSeq(FlexSeq b) {
        _type = b._type;
        _len = b._len;
        _copy_num = b._copy_num;
        if (b._seq == null) {
            _seq = null;
        } else {
            _seq = b._seq.clone();
        }
    }

    /**
     * Return byte at some index in the sequence
     *
     * @param i index in sequence
     * @return
     */
    public byte charAt(int i) {
        if (_seq != null) {
            return _seq[i];
        } else {
            return 0;
        }
    }

    /**
     * @param beginIndex start of substring
     * @return array from beginIndex to end
     */
    public byte[] substring(int beginIndex) {
        return Arrays.copyOfRange(_seq, beginIndex, _seq.length);
    }

    /**
     * @param beginIndex start of substring (inclusive)
     * @param endIndex   end of substring (exclusive)
     * @return substring from [beginIndex,endIndex)
     */
    public byte[] substring(int beginIndex, int endIndex) {
        return Arrays.copyOfRange(_seq, beginIndex, endIndex);
    }

    /**
     * @return length of the sequence
     */
    public int length() {
        return _len;
    }

    /**
     * @return length of the sequence accounting for copy numbers
     */
    public int var_length() {
        if (_type == Type.DUP) {
            return _len * _copy_num;
        }
        return _len;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FlexSeq flexSeq = (FlexSeq) o;

        if (_copy_num != flexSeq._copy_num) return false;
        if (_len != flexSeq._len) return false;
        if (!Arrays.equals(_seq, flexSeq._seq)) return false;
        if (_type != flexSeq._type) return false;
        if (_var_id != null ? !_var_id.equals(flexSeq._var_id) : flexSeq._var_id != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = _type != null ? _type.hashCode() : 0;
        result = 31 * result + _len;
        result = 31 * result + _copy_num;
        result = 31 * result + (_seq != null ? Arrays.hashCode(_seq) : 0);
        result = 31 * result + (_var_id != null ? _var_id.hashCode() : 0);
        return result;
    }

    public int getCopy_num() {
        return _copy_num;
    }

    public void setCopy_num(int cn) {
        _copy_num = cn;
    }

    public Type getType() {
        return _type;
    }

    public void setType(Type t) {
        _type = t;
    }

    public String getVar_id() {
        return _var_id;
    }

    public void setVar_id(String var_id) {
        _var_id = var_id;
    }

    public void setLength(int l) {
        _len = l;
    }

    public byte[] getSeq() {
        return _seq;
    }

    // has the actual sequence
    public boolean isSeq() {
        return (_type == Type.SEQ);
    }

    public String toString() {
        if (_seq == null) {
            switch (_type) {
                case DUP:
                    return "<DUP:TANDEM>";
                case INS:
                    return "<INS>";
                case INV:
                    return "<INV>";
                case DEL:
                    return "<DEL>";
                case TRANSLOCATION:
                    return "<TRA>";
                default:
                    break;
            }
        } else {
            try {
                return new String(_seq, "US-ASCII");
            } catch (UnsupportedEncodingException e) {
                e.printStackTrace();
            }
        }
        return "";
    }

    public boolean equals(final FlexSeq seq) {
        if (_type != seq._type) {
            return false;
        }
        if (_len != seq._len) {
            return false;
        }
        if (_copy_num != seq._copy_num) {
            return false;
        }

        if (!Arrays.equals(_seq, seq._seq)) {
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
     * TRANSLOCATION = translocation (unknown sequence, sequence is defined at variant level)
     */
    public enum Type {
        SEQ, INV, DUP, INS, DEL, TRANSLOCATION
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

    public void setEnd2(int end2) {
        this.end2 = end2;
    }

    public void setReferenceAlleleLength(int referenceAlleleLength) {
        this.referenceAlleleLength = referenceAlleleLength;
    }

    public int getReferenceAlleleLength() {

        return referenceAlleleLength;
    }
}
