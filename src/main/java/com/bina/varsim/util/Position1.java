package com.bina.varsim.util;

/**
 * A class to store a base-1 location.
 */
public final class Position1 extends Position {
    private static final Integer BASE = 1;

    public Position1(long pos) {
        super(pos);
        if (pos < 1) throw new RuntimeException("position " + pos + " is smaller than 1");
    }

    public Position1(Position other) {
        super(other.longValue() - other.base() + BASE);
    }

    @Override
    public int base() { return BASE; }
}
