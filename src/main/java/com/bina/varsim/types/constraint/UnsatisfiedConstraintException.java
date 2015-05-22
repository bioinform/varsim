package com.bina.varsim.types.constraint;

import java.util.Collection;

/**
 * Created by johnmu on 5/21/15.
 */
public class UnsatisfiedConstraintException extends Exception{
    public static class valuePair{
        final Constraint constraint;
        final double actualValue;

        public valuePair(Constraint constraint, double actualValue) {
            this.constraint = constraint;
            this.actualValue = actualValue;
        }

        public Constraint getConstraint() {
            return constraint;
        }

        public double getActualValue() {
            return actualValue;
        }
    }

    final protected Collection<valuePair> constraints;

    public UnsatisfiedConstraintException(Collection<valuePair> constraints){
        this.constraints = constraints;
    }

    public Collection<valuePair> getConstraints() {
        return constraints;
    }
}
