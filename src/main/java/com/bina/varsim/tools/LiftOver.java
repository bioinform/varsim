package com.bina.varsim.tools;

import com.bina.varsim.VarSimTool;
import com.bina.varsim.VarSimToolNamespace;
import com.bina.varsim.fastqLiftover.types.GenomeInterval;
import com.bina.varsim.fastqLiftover.types.MapBlock;
import com.bina.varsim.fastqLiftover.types.MapBlocks;
import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.ReadMapBlock;
import htsjdk.tribble.annotation.Strand;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.util.Collection;

/**
 * liftover files to desired coordinates
 */
public class LiftOver extends VarSimTool {
  private final static Logger log = Logger.getLogger(LiftOver.class.getName());

  @Option(name = "-input", usage = "input file to be lifted", metaVar = "file", required = true)
  String input;

  @Option(name = "-map", usage = "map file", metaVar = "file", required = true)
  String mapForLiftover;

  @Option(name = "-prefix", usage = "output prefix", required = true)
  String prefix;

  public LiftOver(final String command, final String description) {
    super(command, description);
  }

  public static void main(String[] args) {
    new LiftOver("", VarSimToolNamespace.LiftOver.description).run(args);
  }

  /**
   * Main method
   */
  public void run(String[] args){
    if (!parseArguments(args)) {
      return;
    }
    log.info("command: " + String.join(" ", args));
    if (!input.endsWith(".bed")) {
      throw new IllegalArgumentException("Right now only takes .bed files.");
    }
    String outfile = prefix + ".bed";
    try (BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outfile));
         BufferedReader fileToLift = new BufferedReader(new FileReader(input))) {

        MapBlocks mapBlocks = new MapBlocks(new File(mapForLiftover));
        log.info("Assume forward strand if not set.");

        String line;
        while ((line = fileToLift.readLine()) != null) {
          String[] fields = line.split("\t");
          GenomeInterval interval = new GenomeInterval(new ChrString(fields[0]), Integer.parseInt(fields[1]), Integer.parseInt(fields[2]),
                  fields.length >= 6 ? Strand.decode(fields[5]) : Strand.FORWARD, MapBlock.BlockType.UNKNOWN);
          Collection<ReadMapBlock> liftedReadMapBlocks = mapBlocks.liftOverGenomeInterval(interval, 1);
          for (ReadMapBlock i : liftedReadMapBlocks) {
            GenomeInterval liftedInterval = i.getMapInterval();
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder.append(liftedInterval.getChromosome().toString());
            stringBuilder.append("\t");
            stringBuilder.append(liftedInterval.getStart());
            stringBuilder.append("\t");
            stringBuilder.append(liftedInterval.getEnd());
            // try to write whatever original input file has
            if (fields.length >= 4) {
              stringBuilder.append("\t");
              stringBuilder.append(fields[3]);
            }
            if (fields.length >= 5) {
              stringBuilder.append("\t");
              stringBuilder.append(fields[4]);
            }
            if (fields.length >= 6) {
              stringBuilder.append("\t");
              stringBuilder.append(liftedInterval.getStrand().toString());
            }
            if (fields.length >= 7) {
              for (int j = 6; j < fields.length; j++) {
                stringBuilder.append("\t");
                stringBuilder.append(fields[j]);
              }
            }
            stringBuilder.append("\n");
            outputWriter.write(stringBuilder.toString());
          }
        }
    } catch(IOException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Unable to write to " + outfile + " or read " + mapForLiftover + " or " + input);
    }
  }
}
