package com.bina.varsim.util;

/**
 * A class to store a base-0 location.
 */
public final class Position0 extends Position {
    public static final int BASE = 0;

    public Position0(long pos) {
        super(pos);
        if (pos < this.base()) throw new RuntimeException("position " + pos + " is smaller than base "+this.base());
    }

    public Position0(Position other) {
        super(other.longValue() - other.base() + BASE);
    }

    /**
     * @return a new instance of Position0 representing the provided base-b position p
     */
    public static Position0 valueOf(long p, int b) {
        return new Position0(p - b + BASE);
    }

    /**
     * @return a new instance of Position0 representing the provided base-b position p
     */
    public static Position0 valueOf(String p, int b) {
        return new Position0(Long.valueOf(p) - b + BASE);
    }

    @Override
    public int base() { return BASE; }
}
