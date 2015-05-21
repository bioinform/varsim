package com.bina.varsim.types.constraint;

import java.util.Collection;

/**
 * Created by johnmu on 5/21/15.
 */
public class UnsatisfiedConstraintException extends Exception{
    final protected Collection<Constraint> constraints;
    public UnsatisfiedConstraintException(Collection<Constraint> constraints){
        this.constraints = constraints;
    }

    public Collection<Constraint> getConstraints() {
        return constraints;
    }
}
