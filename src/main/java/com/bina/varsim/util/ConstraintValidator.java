package com.bina.varsim.util;

import com.bina.varsim.types.constraint.Constraint;
import com.bina.varsim.types.constraint.ConstraintRecord;
import com.bina.varsim.types.stats.StatsNamespace;
import com.bina.varsim.types.variant.VariantOverallType;

import java.util.ArrayList;
import java.util.Collection;

/**
 * This collects statistics for each constraint and determines if the constraint is satisfied
 */
public class ConstraintValidator {
    final ArrayList<ConstraintRecord> constraintRecords = new ArrayList<>();

    public ConstraintValidator() {
    }

    public ConstraintValidator(Collection<Constraint> constraints){
        for (Constraint constraint : constraints) {
            constraintRecords.add(new ConstraintRecord(constraint));
        }
    }

    public boolean isValid(){
        for (ConstraintRecord record : constraintRecords) {
            if(!record.isValid()) return false;
        }
        return true;
    }

    public void inc(StatsNamespace stats, VariantOverallType type, long len){
        for (ConstraintRecord record : constraintRecords) {
            record.inc(stats, type, len);
        }
    }
}
