package com.bina.varsim.types.constraint;

import com.bina.intervalTree.SimpleInterval1D;
import com.bina.varsim.types.variant.VariantOverallType;
import org.apache.commons.lang3.StringUtils;

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
     *
     * @param constraintArg
     */
    public Constraint(String constraintArg) {
        String[] tokens = StringUtils.split(constraintArg, '/');

        if(tokens.length != 5) throw new RuntimeException(String.format("Invalid constraint: %s", constraintArg));

        // variant type
        varType = VariantOverallType.fromString(tokens[0]);

        // Size range
        range = parseRange(tokens[1]);

        // Metric for comparison
        metric = Metrics.fromString(tokens[2]);

        // Operation for comparison
        operation = Operation.valueOf(tokens[3]);

        // Cutoff to throw an error
        cutoff = Double.parseDouble(tokens[4]);
    }

    public static SimpleInterval1D parseRange(String s) {
        String[] tokens = StringUtils.split(s, '-');

        if (tokens.length != 2) {
            throw new RuntimeException(String.format("Invalid Range: %s", s));
        }

        return new SimpleInterval1D(Long.parseLong(tokens[0]), Long.parseLong(tokens[1]));
    }

    @Override
    public String toString() {
        return "Constraint{" +
                "varType=" + varType +
                ", range=" + range +
                ", metric=" + metric +
                ", operation=" + operation +
                ", cutoff=" + cutoff +
                '}';
    }

    public VariantOverallType getVarType() {
        return varType;
    }

    public SimpleInterval1D getRange() {
        return range;
    }

    public Metrics getMetric() {
        return metric;
    }

    public Operation getOperation() {
        return operation;
    }

    public double getCutoff() {
        return cutoff;
    }

    // This method checks if a given value satisfies a condition based on the
    // operation and cutoff point provided.
    public boolean isSatisfied(double value) {
        switch (operation) {
            case GT:
                return value > cutoff;
            case LT:
                return value < cutoff;
            default:
                throw new RuntimeException("Null operation");
        }
    }
}
