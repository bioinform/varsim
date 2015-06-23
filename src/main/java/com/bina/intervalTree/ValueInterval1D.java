package com.bina.intervaltree;

/**
 * Created by johnmu on 12/17/14.
 *
 * @author johnmu
 */
public class ValueInterval1D<Value> extends SimpleInterval1D {
    private Value data;

    public ValueInterval1D(long left, long right) {
        super(left, right);
        data = null;
    }

    public ValueInterval1D(long left, long right, Value data) {
        super(left, right);
        this.data = data;
    }

    public ValueInterval1D(SimpleInterval1D reg, Value data) {
        super(reg);
        this.data = data;
    }

    public Value get() {
        return data;
    }

    public void set(Value data) {
        this.data = data;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        ValueInterval1D that = (ValueInterval1D) o;

        if (data != null ? !data.equals(that.data) : that.data != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (data != null ? data.hashCode() : 0);
        return result;
    }
}
