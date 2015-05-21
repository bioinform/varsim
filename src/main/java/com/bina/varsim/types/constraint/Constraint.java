package com.bina.varsim.types.constraint;

import com.bina.varsim.intervalTree.SimpleInterval1D;
import com.bina.varsim.types.variant.VariantOverallType;
import org.apache.commons.lang3.StringUtils;

import java.util.Collection;

/**
 * Created by johnmu on 5/20/15.
 */
public class Constraint {
    protected final VariantOverallType varType;
    protected final SimpleInterval1D range;
    protected final Metrics metric;
    protected final Operation operation;
    protected final double cutoff;

    /**
     * Will look something like Insertion/20-30/F1/GT/0.99
     * Means Insertions of size 20-30 (inclusive) must have F1 score greater than 0.99
     * @param constraintArg
     */
    public Constraint(String constraintArg) {
        String[] tokens = StringUtils.split(constraintArg, '/');


    }
}
