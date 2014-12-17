package com.binatechnologies.varsim;

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

    public ValueInterval1D(long left, long right, Value data){
        super(left, right);
        this.data = data;
    }

    public ValueInterval1D(SimpleInterval1D reg, Value data){
        super(reg);
        this.data = data;
    }

    public Value get() {
        return data;
    }

    public void set(Value data) {
        this.data = data;
    }


}
