package com.bina.varsim;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.IOException;

/**
 * Created by mohiyudm on 12/9/16.
 */
public abstract class VarSimTool {
    public final String VERSION = "VarSim" + " " + getClass().getPackage().getImplementationVersion();

    @Option(name = "-version", usage = "Print version", help=true)
    protected boolean printVersion = false;

    public abstract void run(String[] args) throws IOException;
}
