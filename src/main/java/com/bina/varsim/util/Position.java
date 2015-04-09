package com.bina.varsim.util;

/**
 * base class for base-0 or base-1 location
 */
public abstract class Position extends Number{
    private final long pos_;

    // only derived class can call constructor
    protected Position(long pos) { pos_ = pos; }

    /**
     * @return the base in which position is represented
     */
    public abstract int base();

    /**
     * @return the location in base-b
     */
    public long asBase(long b) { return pos_ - this.base() + b; }

    /**
     * @return the string representation in base-b
     */
    public String toStringAsBase(long b) { return Long.toString(this.asBase(b)); }

    @Override
    public int intValue() {
        if (pos_ > Integer.MAX_VALUE) throw new RuntimeException("position " + pos_ + " out of range");
        return (int)pos_;
    }

    @Override
    public long longValue() { return pos_; }

    @Override
    public float floatValue() { throw new UnsupportedOperationException(); }

    @Override
    public double doubleValue() { throw new UnsupportedOperationException(); }

    @Override
    public byte byteValue() {
        if (pos_ > Byte.MAX_VALUE) throw new RuntimeException("position " + pos_ + " out of range");
        return (byte) pos_;
    }

    @Override
    public short shortValue() {
        if (pos_ > Short.MAX_VALUE) throw new RuntimeException("position " + pos_ + " out of range");
        return (short) pos_;
    }
}
