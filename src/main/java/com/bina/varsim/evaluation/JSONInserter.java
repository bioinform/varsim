package com.bina.varsim.evaluation;

import org.apache.commons.lang3.StringEscapeUtils;
import org.apache.commons.lang3.text.StrSubstitutor;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * Created by johnmu on 12/3/14.
 */
public class JSONInserter {
    private final static Logger log = Logger.getLogger(JSONInserter.class.getName());

    @Argument(usage = "One or more JSON files from VarSim", metaVar = "json_files ...", required = true)
    private ArrayList<String> jsonFilename = new ArrayList<>();

    @Option(name = "-html", usage = "VarSim HTML to insert the JSON into", metaVar = "file", required = true)
    File html_file;

    public void run(String[] args) throws IOException {
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

        for (String filename : jsonFilename) {
            // generate output_filename
            String[] prefixes = filename.split("\\.");
            String prefix;
            if (prefixes.length == 0) {
                prefix = filename;
            } else {
                prefix = prefixes[0];
            }
            String outName = prefix + "_" + html_file.getName();
            File outFile = new File(outName);

            // construct string to insert
            String jsonStr = StringEscapeUtils.escapeEcmaScript(FileUtils.readFileToString(new File(filename)).trim().replace("\n", "").replace("\r", ""));
            FileUtils.writeStringToFile(outFile, insertJSON(varsim_html, jsonStr));
        }
    }

    public static String insertJSON(String varsim_html, String jsonStr) {
        TreeMap<String, String> lookup = new TreeMap<>();
        lookup.put("varsim_data", "var varsim_data = \"" + StringEscapeUtils.escapeEcmaScript(jsonStr) + "\";");
        return new StrSubstitutor(lookup, "<!--", "-->").replace(varsim_html);
    }

    public static void main(String[] args) throws IOException {
        new JSONInserter().run(args);
    }
}
