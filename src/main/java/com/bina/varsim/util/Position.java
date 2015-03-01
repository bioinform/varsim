package com.bina.varsim.util;

/**
 * Created by bayo on 2/28/15.
 * base class for base-0 or base-1 location
 * in C++ this can be done with template<char val> class Position; I'm not so sure about java
 */
public abstract class Position extends Number{
    private long pos_;

    // only derived class can call constructor
    // in C++ this whole mess can be avoided with template<char val> class Position;
    // I'm not so sure about java
    protected Position(long pos) {
        pos_ = pos;
    }

    public abstract int base();

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
