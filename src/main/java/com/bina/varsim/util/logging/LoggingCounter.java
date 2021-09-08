package com.bina.varsim.util.logging;

import org.apache.log4j.Logger;

public class LoggingCounter {
    private int countDown;
    private final int maxWarnings;
    private final Logger log = Logger.getLogger(LoggingCounter.class.getName());
    public LoggingCounter(int maxWarnings) {
        this.countDown =  maxWarnings;
        this.maxWarnings =  maxWarnings;
    }

    public boolean isCountLeftAndDecrement() {
        if (this.countDown > 0) {
            this.countDown--;
            if (this.countDown == 0) {
                log.warn("Reached max number of warnings (" + this.maxWarnings +
                        "), no more warnings.");
            }
            return true;
        }
        return false;
    }
}
