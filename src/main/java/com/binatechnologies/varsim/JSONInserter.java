package com.binatechnologies.varsim;

import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.util.ArrayList;

/**
 * Created by johnmu on 12/3/14.
 */
public class JSONInserter {
    private final static Logger log = Logger.getLogger(JSONInserter.class.getName());

    @Argument(usage = "One or more JSON files from VarSim",metaVar = "json_files ...",required = true)
    private ArrayList<String> json_filename = new ArrayList<>();

    @Option(name = "-html", usage = "HTML to insert the JSON into" ,metaVar = "file",required = true)
    String true_vcf_filename;

    public void run(String[] args){
        String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
        String usage = "Inserts n JSON files to one HTML to create n HTML files\n";

        CmdLineParser parser = new CmdLineParser(this);

        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java -jar VarSim.jar json_inserter [options...] vcf_files ...");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }
    }

    public static void main(String[] args){
        new JSONInserter().run(args);
    }
}
