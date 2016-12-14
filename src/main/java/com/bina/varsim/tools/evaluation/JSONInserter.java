package com.bina.varsim.tools.evaluation;

import com.bina.varsim.VarSimTool;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringEscapeUtils;
import org.apache.commons.lang3.text.StrSubstitutor;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * Created by johnmu on 12/3/14.
 */
public class JSONInserter extends VarSimTool {
    private final static Logger log = Logger.getLogger(JSONInserter.class.getName());
    @Option(name = "-html", usage = "VarSim HTML to insert the JSON into", metaVar = "file", required = true)
    File html_file;
    @Argument(usage = "One or more JSON files from VarSim", metaVar = "json_files ...", required = true)
    private ArrayList<String> jsonFilename = new ArrayList<>();

    public JSONInserter(final String command, final String description) {
        super(command, description);
    }

    public static String insertJSON(String varsim_html, String jsonStr) {
        TreeMap<String, String> lookup = new TreeMap<>();
        lookup.put("varsim_data", "var varsim_data = \'" + StringEscapeUtils.escapeEcmaScript(jsonStr) + "\';");
        return new StrSubstitutor(lookup, "<!--", "-->").replace(varsim_html);
    }

    public static void main(String[] args) throws IOException {
        new JSONInserter("", "").run(args);
    }

    public void run(String[] args) throws IOException {
        if (!parseArguments(args)) {
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
            String jsonStr = FileUtils.readFileToString(new File(filename)).trim().replace("\n", "").replace("\r", "");
            FileUtils.writeStringToFile(outFile, insertJSON(varsim_html, jsonStr));
        }
    }
}
