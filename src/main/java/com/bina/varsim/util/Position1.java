package com.bina.varsim.util;

/**
 * A class to store a base-1 location.
 */
public final class Position1 extends Position {
    private static final Integer BASE = 1;

    public Position1(long pos) {
        super(pos);
        if (pos < this.base()) throw new RuntimeException("position " + pos + " is smaller than base "+this.base());
    }

    public Position1(Position other) {
        super(other.longValue() - other.base() + BASE);
    }

    /**
     * @return a new instance of Position1 representing the provided base-b position p
     */
    public static Position1 valueOf(long p, int b) {
        return new Position1(p - b + BASE);
    }

    /**
     * @return a new instance of Position1 representing the provided base-b position p
     */
    public static Position1 valueOf(String p, int b) {
        return new Position1(Long.valueOf(p) - b + BASE);
    }

    @Override
    public int base() { return BASE; }
}
