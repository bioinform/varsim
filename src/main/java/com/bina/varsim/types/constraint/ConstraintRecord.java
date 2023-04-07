package com.bina.varsim.types.constraint;

import com.bina.varsim.types.stats.RatioRecord;
import com.bina.varsim.types.stats.StatsNamespace;
import com.bina.varsim.types.variant.VariantOverallType;

/**
 * A RatioRecord for each constraint, also able to validate the constraint
 */
public class ConstraintRecord {
    final RatioRecord stats;
    final Constraint constraint;

    public ConstraintRecord(Constraint constraint) {
        this.constraint = constraint;
        stats = new RatioRecord();
    }

    public ConstraintRecord(String constraintArg) {
        this.constraint = new Constraint(constraintArg);
        stats = new RatioRecord();
    }

    public void inc(StatsNamespace stat, VariantOverallType type, long len){
        if(constraint.getVarType() == type && constraint.getRange().contains(len)){
            switch (stat) {
                case TP:
                    stats.incTP();
                    break;
                case FP:
                    stats.incFP();
                    break;
                case TN:
                    stats.incTN();
                    break;
                case FN:
                    stats.incFN();
                    break;
                case T:
                    stats.incT();
                    break;
                default:
                    throw new RuntimeException("Null stats");
            }
        }
    }


    /**
     * Check if the constraint is satisfied by the stats
     * @return
     */
    public boolean isValid(){
        return constraint.isSatisfied(getStatsValue());
    }

    public double getStatsValue(){
        switch (constraint.getMetric()){
            case F1:
                double precision = 1 - stats.getFDR();
                double recall = stats.getTPR();
                return 2 * (precision * recall)/(precision + recall);
            case TPR:
                return stats.getTPR();
            case FDR:
                return stats.getFDR();
            default:
                throw new RuntimeException("Null metric");
        }
    }

    public Constraint getConstraint() {
        return constraint;
    }
}
