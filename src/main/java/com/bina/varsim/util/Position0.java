package com.bina.varsim.util;

/**
 * A class to store a base-0 location.
 */
public final class Position0 extends Position {
    public static final int BASE = 0;

    public Position0(long pos) {
        super(pos);
        if (pos < 0) throw new RuntimeException("position " + pos + " is smaller than 0");
    }

    public Position0(Position other) {
        super(other.longValue() - other.base() + BASE);
    }

    @Override
    public int base() { return BASE; }
}
