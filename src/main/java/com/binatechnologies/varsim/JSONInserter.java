package com.binatechnologies.varsim;

import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by johnmu on 12/3/14.
 */
public class JSONInserter {
    private final static Logger log = Logger.getLogger(JSONInserter.class.getName());

    @Argument(usage = "One or more JSON files from VarSim",metaVar = "json_files ...",required = true)
    private ArrayList<String> json_filename = new ArrayList<>();

    @Option(name = "-html", usage = "VarSim HTML to insert the JSON into" ,metaVar = "file",required = true)
    File html_file;

    public void run(String[] args) throws IOException{
        String VERSION = "VarSim " + getClass().getPackage().getImplementationVersion();
        String usage = "Inserts n JSON files to one HTML to create n HTML files\n";

        CmdLineParser parser = new CmdLineParser(this);

        parser.setUsageWidth(80);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.println(VERSION);
            System.err.println(e.getMessage());
            System.err.println("java -jar VarSim.jar json_inserter [options...] json_files ...");
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println(usage);
            return;
        }

        String varsim_html = FileUtils.readFileToString(html_file);

        for(String filename : json_filename){
            // generate output_filename
            String outName = html_file.getName() + filename.split(".")[0];
            File file = new File(outName);

            FileUtils.writeStringToFile(file, blah);
        }
    }

    public static void main(String[] args) throws IOException{
        new JSONInserter().run(args);
    }
}
