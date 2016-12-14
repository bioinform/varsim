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

    @Option(name = "-version", usage = "Print version", help = true)
    protected boolean printVersion = false;

    @Option(name = "-h", usage = "Print help message", help = true, aliases = {"-help"})
    protected boolean help = false;

    protected String command = "";

    protected String description = "";

    public VarSimTool(final String command, final String description) {
        this.command = command;
        this.description = description;
    }

    public abstract void run(String[] args) throws IOException;

    public String getCommand() {
        return command;
    }

    public String getDescription() {
        return description;
    }

    public void printUsage(final CmdLineParser parser) {
        System.err.println(VERSION);
        System.err.println("java -jar VarSim.jar " + getCommand() + " [options...]");
        System.err.println(getDescription());
        // print the list of available options
        parser.printUsage(System.err);
    }

    public boolean parseArguments(final String[] args) {
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            printUsage(parser);
            return false;
        }

        if (help) {
            printUsage(parser);
            return false;
        }

        if (printVersion) {
            System.out.println(VERSION);
            return false;
        }

        return true;
    }
}
