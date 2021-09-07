package com.bina.varsim.util;

//--- Java imports ---

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;
import com.bina.varsim.types.Sequence;
import com.bina.varsim.types.VCFInfo;
import com.bina.varsim.types.variant.Variant;
import com.bina.varsim.types.variant.alt.Alt;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.rmi.UnexpectedException;
import java.util.*;

import static com.bina.varsim.constants.Constant.MAX_WARNING_REPEAT;
import static com.bina.varsim.types.VCFInfo.getType;

public class VCFparser extends GzFileParser<Variant> {
    public static final String DEFAULT_FILTER = "."; //default value for many columns
    private final static Logger log = Logger.getLogger(VCFparser.class.getName());

    private Random random = null;

    private int sampleIndex = -1;

    public String getSampleId() {
        return sampleId;
    }

    private String sampleId = null;
    private boolean isPassFilterRequired = false;
    private boolean chromLineSeen = false;
    private int illegalPhasingWarningCount = 0;

    public VCFparser() {
        sampleIndex = 10; // the first sample
    }

    /**
     * Reads a VCF file line by line
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param id       ID of individual, essentially selects a column, null to use first ID column
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, String id, boolean pass, Random rand) {
        this(new File(fileName), id, pass, rand);
        log.info("Reading " + fileName);
    }


    public VCFparser(File file, String id, boolean pass, Random rand) {
        random = rand;
        try {
            bufferedReader = new BufferedReader(new InputStreamReader(decompressStream(file)));
            readLine();
        } catch (Exception ex) {
            log.error("Can't open file " + file.getName());
            log.error(ex.toString());
        }
        sampleId = id;
        isPassFilterRequired = pass;

        if (StringUtils.isEmpty(sampleId)) {
            sampleIndex = 10; // the first sample
        }
    }

    /**
     * Reads a VCF file line by line
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param id       ID of individual, essentially selects a column, null to use first ID column
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, String id, boolean pass) {
        this(fileName, id, pass, null);
    }

    public VCFparser(File file, String id, boolean pass) {
        this(file, id, pass, null);
    }

    //TODO: remove unused constructor
    /**
     * Reads a VCF file line by line, if there are multiple individuals, takes the first one
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, boolean pass, Random rand) {
        this(fileName, null, pass, rand);
    }

    //TODO: remove unused constructor
    /**
     * Reads a VCF file line by line, if there are multiple individuals, takes the first one
     *
     * @param fileName VCF file, doesn't have to be sorted or indexed
     * @param pass     If true, only output pass lines
     */
    public VCFparser(String fileName, boolean pass) {
        this(fileName, null, pass, null);
    }

    public VCFparser(File file, boolean pass) {
        this(file, null, pass, null);
    }

    /**
     * finds where "GT" or similar is in the VCF string so that the genotype can be read
     *
     * @param record The column in VCF that contains GT, CN, etc..
     * @param key    The key to be found "GT" or "CN" or other
     * @return the index of the key
     */
    private int getFormatKeyIndex(final String record, final String key) {
        StringTokenizer words = new StringTokenizer(record, ":");
        int ret = 0;
        while (words.hasMoreTokens()) {
            if (words.nextToken().equals(key))
                return ret;
            else
                ret++;
        }
        return -1;
    }

    /**
     * Takes genotype string and splits it into alleles, supports a maximum of two
     *
     * @param geno genotype string corresponding to the GT tag (field index 9)
     * @param vals Return the genotypes read here [paternal,maternal]
     * @param chr  chromosome we are dealing with, some are haploid, need to assign parent
     * @return true if the variant is phased
     */
    boolean isPhased(String geno, byte[] vals, ChrString chr) {
        boolean isPhased = false;

        geno = geno.trim();
        boolean strangePhase = false;
        String[] ll = StringUtilities.fastSplit(geno, "/|");

        if (ll.length == 1) {
            // phase is only a single number, for haploid chromosomes
            byte val = (byte) (ll[0] == "." ? -1 : StringUtilities.parseInt(ll[0]));
            if (chr.isX()) {
                vals[0] = -1; //paternal missing
                vals[1] = val; // maternal
                isPhased = true;
            } else if (chr.isY()) {
                vals[0] = val; // paternal
                vals[1] = -1; //maternal missing
                isPhased = true;
            } else if (chr.isMT()) {
                vals[0] = -1; //paternal missing
                vals[1] = val;
                isPhased = true;
            } else {
                vals[0] = vals[1] = val;
            }
        } else if (ll.length == 2) {
            // this is the case where phase looks like "1|0" or "10|4"
            int c1 = -1;
            int c2 = -1;
            char phasing = geno.charAt(ll[0].length());
            try {
                c1 = ll[0] == "." ? -1 : StringUtilities.parseInt(ll[0]);
                c2 = ll[1] == "." ? -1 : StringUtilities.parseInt(ll[1]);
            } catch (NumberFormatException e) {
                strangePhase = true;
            }

            if (phasing == '|') {
                isPhased = true;
            }
            if ((c1 >= 0) == (c2 >= 0)) {
                vals[0] = (byte) c1;
                vals[1] = (byte) c2;
            } else {
                strangePhase = true;
            }
        } else {
                strangePhase = true;
        }

        if (strangePhase) {
            if (illegalPhasingWarningCount < MAX_WARNING_REPEAT) {
                log.warn("Unrecognized phasing '" + geno + "'.");
                illegalPhasingWarningCount++;
                if (illegalPhasingWarningCount == MAX_WARNING_REPEAT) {
                    log.warn("Reached max number of warnings (" + MAX_WARNING_REPEAT +
                    ") for unrecognized phasing. No more warnings.");
                }
            }
            vals[0] = -1;
            vals[1] = -1;
            isPhased = false;
        }

        return isPhased;
    }

    /**
     * takes a line from a VCF file, parse it,
     * return a Variant object
     *
     * right now the meta-info lines (beginning with ##) are
     * not tied with data line parsing. This will be corrected
     * in the future (perhaps with help of HTSJDK).
     * @param line
     * @return
     */
    public Variant processLine(String line) throws UnexpectedException {

        // try to determine the column we should read for the genotype
        String[] toks = StringUtilities.fastSplit(line, "\t");
        if (line.startsWith("#")) {
            if (!StringUtils.isEmpty(sampleId) && line.startsWith("#CHROM")) {
                chromLineSeen = true;
                int index = 0;
                for (String tok : toks) {
                    index++;
                    if (tok.equals(sampleId))
                        sampleIndex = index;
                }
            } else if (StringUtils.isEmpty(sampleId)) {
                sampleIndex = 10; // the first sample
                if (line.startsWith("#CHROM") && toks.length >= sampleIndex) {
                    sampleId = toks[sampleIndex - 1];
                }
            }
            return null;
        }

        // If we cannot determine, then use the first one
        if (sampleIndex < 0 && !chromLineSeen) {
            sampleIndex = 10;
        } else if (sampleIndex < 0) {
            sampleIndex = 10;
            log.warn("Warning!!! ID (" + sampleId + ") does not exist... ");
        }


        int index = 0, genotypeIndex = -1, copyNumberIndex = -1;
        int pos = -1;
        ChrString chr = null;
        String REF = "", FILTER = "", ALT = "", variantId = "";
        String QUAL = ".";
        String phase = ".", copyNumber = "0/0", infoString = "", FORMAT;
        String[] sampleInfo;
        for (String tok : toks) {
            index++;
            if (index == 1) { // Parsing chromosome
                chr = new ChrString(tok);
            } else if (index == 2) // Parsing position
                pos = StringUtilities.parseInt(tok);
            else if (index == 3) // Parsing position
                variantId = tok;
            else if (index == 4) // Parsing reference allele
                REF = tok;
            else if (index == 5) // Parsing alternative allele
                ALT = tok;
            else if (index == 6) // Parsing alternative allele
                QUAL = tok;
            else if (index == 7) // FILTER field
                FILTER = tok;
            else if (index == 8) // INFO field
                infoString = tok;
            else if (index == 9) { // Output format
                FORMAT = tok;
                genotypeIndex = getFormatKeyIndex(FORMAT, "GT");
                copyNumberIndex = getFormatKeyIndex(FORMAT, "CN");
            } else if (index == sampleIndex) { // phased or unphased genotype
                sampleInfo = StringUtilities.fastSplit(tok, ":");
                if (genotypeIndex >= 0) {
                    phase = sampleInfo[genotypeIndex];
                }
                if (copyNumberIndex >= 0) {
                    copyNumber = sampleInfo[copyNumberIndex];
                }

                break;
            }
        }

        // unknown chromosome
        // TODO: throw an exception for unknown chromosome name
        if (chr == null) {
            log.warn("Bad chromosome name.");
            return null;
        }

        if (isPassFilterRequired && !(FILTER.contains("PASS") || FILTER.equals(DEFAULT_FILTER))) {
            log.warn("line is filtered out.");
            return null; // Filtered out
        }
        // parse the phased or unphased genotype
        byte[] genotypeArray = new byte[2]; // paternal-maternal
        boolean isGenotypePhased = isPhased(phase, genotypeArray, chr);


        if (genotypeIndex >= 0 && genotypeArray[0] == 0 && genotypeArray[1] == 0) {
            log.warn("All ALT alleles are reference sequences.");
            return null; // reference alleles... ignore them for now....
        }

        if (!Sequence.isNormalSequence(REF)) {
            log.warn("only ATCGN (case-insensitive) allowed for REF column.");
            return null; //
        }

        // determine copy-number
        // TODO need to be able to deal with unphased copy-numbers?
        byte[] copyNumberArray = new byte[2]; // paternal-maternal
        boolean isCopyNumberPhased;

        if (copyNumberIndex >= 0) {
            isCopyNumberPhased = isPhased(copyNumber, copyNumberArray, chr);
            if (isCopyNumberPhased != isGenotypePhased) {
                // TODO maybe don't throw error, this is not standard format
                // anyways
                log.warn("Inconsistent copy number.");
                return null;
            }
        }

        // Upper casing
        REF = REF.toUpperCase();
        ALT = ALT.toUpperCase();

        String deletedReference = "";
        VCFInfo info = new VCFInfo(infoString);

        /*if symbolic alleles are present,
        make sure # of alleles equal # of
        SV lengths. For non-symbolic alleles, SV lengths
        are not really used or checked.
         */
        /*!!!!!!!!!!
        CAUTION: we assume symbolic alleles are not mixed
        with non-symbolic alleles.
         */
        if (ALT.indexOf('<') != -1) {
            String[] alternativeAlleles = StringUtilities.fastSplit(ALT, ",");
            int[] svlen = info.getValue("SVLEN", int[].class);
            if (alternativeAlleles.length != svlen.length) {
                throw new IllegalArgumentException("ERROR: number of symbolic alleles is unequal to number of SV lengths.\n" + line);
            }
            for (int i = 0; i < alternativeAlleles.length; i++) {
                if (!alternativeAlleles[i].startsWith("<")) {
                    throw new IllegalArgumentException("ERROR: symbolic alleles are mixed with non-symbolic alleles.\n" + line);
                }
            }
        }
        Alt[] alts = null;
        try {
            alts = string2Alt(ALT);
        } catch (IllegalArgumentException e) {
            log.warn("ALT column is malformated: " + e.getMessage());
            return null;
        }

      if (alts[0].getSymbolicAllele() != null) {
          int[] end = info.getValue("END", int[].class);
          //SVLEN for alternative allele length
          int[] svlen =  info.getValue("SVLEN", int[].class);
          int[] end2 =  info.getValue("END2", int[].class);
          int[] pos2 =  info.getValue("POS2", int[].class);
          Boolean isinv =  info.getValue("ISINV", Boolean.class);
          Boolean isLengthImprecise =  info.getValue("IMPRECISE_LENGTH", Boolean.class);
          isLengthImprecise = isLengthImprecise == null? false : isLengthImprecise;
          String[] traid =  info.getValue("TRAID", String[].class);
          String[] chr2 =  info.getValue("CHR2", String[].class);
          deletedReference = REF;
          byte[] refs = new byte[0];
          pos++; //1-based start

          if (Alt.SVType.SVSubtype.TRA.equals(alts[0].getSymbolicAllele().getMinor())) {
            if (traid == null || traid.length == 0) {
                throw new IllegalArgumentException("ERROR: <*:TRA> must have TRAID in INFO field.\n" + line);
            }
          }

          //why no alts()? because it's not reference-based operation
          //alts() will save a deep copy of alts, rathern than reference, so we
          //cannot assign it until all changes are done.
          Variant.Builder template = new Variant.Builder().chr(chr).pos(pos).
                  ref(refs).phase(genotypeArray).isPhased(isGenotypePhased).
                  varId(variantId).filter(FILTER).refDeleted(deletedReference).
                  isLengthImprecise(isLengthImprecise).
                  randomNumberGenerator(random);

          if (alts[0].getSymbolicAllele().getMajor() == Alt.SVType.INV) {
              // inversion SV
              if (svlen.length > 0) {
                  for (int i = 0; i < svlen.length; i++) {
                      if (i > 0 && svlen[i] != svlen[i - 1]) {
                          log.warn("Right now VarSim does not support multiple SVLEN for <INV>.");
                          return null;
                      }
                      int alternativeAlleleLength = Math.max(Math.abs(svlen[i]), 1);
                      alts[i].setSeq(new FlexSeq(FlexSeq.Type.INV, alternativeAlleleLength));
                  }
                  // TODO this assumes only one alt
                return template.referenceAlleleLength(Math.abs(svlen[0])).alts(alts).build();
                  //TODO: this assumes only one alt, might not be true
              } else if (end != null && end.length > 0 && end[0] > 0) {
                  //assume only one END
                  int alternativeAlleleLength = Math.max(Math.abs(end[0] - pos + 1), 1);
                  alts[0].setSeq(new FlexSeq(FlexSeq.Type.INV, alternativeAlleleLength));

                  return template.referenceAlleleLength(alternativeAlleleLength).alts(alts).build();
              } else {
                  log.error("No length information for INV, skipping...");
                  return null;
              }
          } else if (alts[0].getSymbolicAllele().getMajor() == Alt.SVType.DUP &&
                    ((alts[0].getSymbolicAllele().getMinor() != Alt.SVType.SVSubtype.TRA &&
                            alts[0].getSymbolicAllele().getMinor() != Alt.SVType.SVSubtype.ISP &&
                            info.getValue("POS2", getType("POS2")) == null) ||
                     alts[0].getSymbolicAllele().getMinor() == Alt.SVType.SVSubtype.TANDEM)) {
              if (svlen.length > 0) {
                  for (int i = 0; i < svlen.length; i++) {
                      if (i > 0 && svlen[i] != svlen[i - 1]) {
                          log.warn("Right now VarSim does not handle multiple SVLEN for <DUP>.");
                          return null;
                      }
                      // TODO this is temporary, how to encode copy number?
                      int currentCopyNumber = 2;
                      for (int j = 0; j < 2; j++) {
                          if ((i + 1) == genotypeArray[j]) {
                            /*
                            if i = 0, genotype[0] = 1, genotype[1] = 1
                            copyNumberArray[0] = 3, copyNumberArray[1] = 2
                            then currentCopyNumber = 2.
                            what does currentCopyNumber mean in real world?
                             */
                              if (copyNumberArray[j] > 0) {
                                  currentCopyNumber = copyNumberArray[j];
                              }
                          }
                      }

                      int alternativeAlleleLength = Math.max(Math.abs(svlen[i]), 1);

                      alts[i].setSeq(new FlexSeq(FlexSeq.Type.TANDEM_DUP, alternativeAlleleLength, currentCopyNumber));
                  }

                  return template.referenceAlleleLength(Math.abs(svlen[0])).alts(alts).build();
                  //TODO: this assumes only one alt, which might not be true
              } else if (end != null && end.length > 0 && end[0] > 0) {
                  int alternativeAlleleLength = Math.max(Math.abs(end[0] - pos + 1), 1);
                  alts[0].setSeq(new FlexSeq(FlexSeq.Type.TANDEM_DUP, alternativeAlleleLength, Math.max(
                          copyNumberArray[0], copyNumberArray[1])));

                  return template.referenceAlleleLength(alternativeAlleleLength).alts(alts).build();
              } else {
                  log.error("No length information for DUP:TANDEM, skipping...");
                  return null;
              }
          } else if (alts[0].getSymbolicAllele().getMajor() == Alt.SVType.INS) {
              // insertion SV

              if (svlen.length > 0) {
                  for (int i = 0; i < svlen.length; i++) {
                    //if SVLEN=0, we take it as infinitely long
                      int alternativeAlleleLength = svlen[i] == 0 ? Integer.MAX_VALUE : Math.max(Math.abs(svlen[i]), 1);
                      alts[i].setSeq(new FlexSeq(FlexSeq.Type.INS, alternativeAlleleLength));
                  }
                  return template.referenceAlleleLength(0).alts(alts).build();
                //TODO, remove this as END should be equal to POS for insertion
              } else if (end != null && end.length > 0 && end[0] > 0) {
                  int alternativeAlleleLength = Math.max(Math.abs(end[0] - pos), 1);
                  alts[0].setSeq(new FlexSeq(FlexSeq.Type.INS, alternativeAlleleLength));
                  return template.referenceAlleleLength(0).alts(alts).build();
              } else {
                  log.error("No length information for INS, skipping...");
                  return null;
              }
          } else if (alts[0].getSymbolicAllele().getMajor() == Alt.SVType.DEL) {
              // deletion SV (maybe part of a translocation)
              // but... we don't have the reference... so we add some random sequence?
              template = template.traid(traid == null ? null : traid[0]);

              if (svlen.length > 0) {
                  for (int i = 0; i < svlen.length; i++) {
                      if (i > 0 && svlen[i] != svlen[i-1]) {
                          log.warn("Right now VarSim does not handle multiple SVLEN for <DEL>.");
                          return null;
                      }
                      // deletion has no alt
                      if (Alt.SVType.SVSubtype.TRA == alts[i].getSymbolicAllele().getMinor()) {
                          alts[i].setSeq(new FlexSeq(FlexSeq.Type.TRA_DEL, svlen[i]));
                      } else {
                          alts[i].setSeq(new FlexSeq(FlexSeq.Type.DEL, 0));
                      }
                  }
                  return template.alts(alts).referenceAlleleLength(Math.abs(svlen[0])).build();
              } else if (end != null && end.length > 0 && end[0] > 0) {
                //END is just one value, whereas there could be multiple alternative alleles with different svlens
                  //so END is in general not a good way to get lengths
                  int alternativeAlleleLength = end[0] - pos + 1;
                  if (Alt.SVType.SVSubtype.TRA == alts[0].getSymbolicAllele().getMinor()) {
                      alts[0].setSeq(new FlexSeq(FlexSeq.Type.TRA_DEL, -alternativeAlleleLength));
                  } else {
                      alts[0].setSeq(new FlexSeq(FlexSeq.Type.DEL, 0));
                  }
                  return template.alts(alts).referenceAlleleLength(alternativeAlleleLength).build();
              } else {
                  log.error("No length information for DEL, skipping...");
                  return null;
              }
              //TODO major SVTYPE actually does not allow TRA
          } else if ( alts[0].getSymbolicAllele().getMajor() == Alt.SVType.DUP &&
                  (alts[0].getSymbolicAllele().getMinor() == Alt.SVType.SVSubtype.TRA || alts[0].getSymbolicAllele().getMinor() == Alt.SVType.SVSubtype.ISP ||
                  info.getValue("POS2", getType("POS2")) != null)) {
              //translocation SV DUP or interspersed DUP

              if (svlen.length > 0) {
                  //0 is for reference allele
                  //alternative allele is numbered 1,2,... per VCFSpec
                  for (int altAlleleIndex = 1; altAlleleIndex <= svlen.length; altAlleleIndex++) {
                      int currentCopyNumber = 2;
                    /*
                    implicit assumption here: genotype[0] == genotype[1] => copyNumberArray[0] == copyNumberArray[1]
                     */
                      //check paternal
                      if (altAlleleIndex == genotypeArray[0]) {
                          currentCopyNumber = copyNumberArray[0];
                      }
                      //check maternal
                      if (altAlleleIndex == genotypeArray[1]) {
                          currentCopyNumber = copyNumberArray[1];
                      }
                      currentCopyNumber = Math.max(1, currentCopyNumber);
                      //allow svlen to be negative
                      int altAllelelength = Math.max(Math.abs(svlen[altAlleleIndex - 1]), 1);
                    /*
                    a translocation is decomposed into a duplication (a special translocation) and a deletion
                     */
                      alts[altAlleleIndex - 1].setSeq(new FlexSeq(
                              alts[altAlleleIndex - 1].getSymbolicAllele().getMinor() == Alt.SVType.SVSubtype.TRA ? FlexSeq.Type.TRA_DUP : FlexSeq.Type.ISP_DUP,
                              altAllelelength, currentCopyNumber));
                  }

                  //TODO: there could be multiple SVLEN values, but for now we only use one
                  //make sure all SVLEN values are equal if we only use the first one
                  //per VCFv4.1 spec, SVLEN should be length of alternative allele rather than
                  //reference allele
                  for (int i = 1; i < svlen.length; i++) {
                      if (svlen[i] != svlen[0]) {
                          throw new IllegalArgumentException("ERROR: SVLEN values not equal.\n" + line);
                      }
                  }
                  //pos is incremented by 1, so it becomes 1-based start
                  return template.referenceAlleleLength(0).alts(alts).
                          chr2(ChrString.string2ChrString(chr2)).pos2(pos2).end2(end2).isinv(isinv).traid(traid == null? null : traid[0]).build();
                  //TODO: this assumes only one alt, which might not be true
              } else {
                  log.error("No length information for DUP:TRA or DUP:ISP, skipping...");
                  return null;
              }
          } else {
              // imprecise variant
              log.warn("Imprecise variant.");
              return null;
          }
      } else if (alts[0].getSeq() != null){
          //ALT field contains actual sequence
            // Check
            for (int i = 0; i < alts.length; i++) {
                if (REF.length() == 1 && alts[i].length() == 1) {
                    // SNP
                } else if (REF.length() == 0 || alts[i].length() == 0) {
                    log.warn("Skipping invalid record.");
                    return null;
                }
            }
            /* Adjustment of first base
             basically if first base of ref and alt match, first base of
             ref and alt will both be removed, pos will increment by 1 to
             account for the removal.

            this is updated to account for multiple matching reference bases
            e.g. ref=AT alt=ATC (VCFv4.1 spec requires only 1-bp before event, but it might
            not be the case all the time.

            */

            while (REF.length() > 0) {
                boolean same = true;
                for (int i = 0; i < alts.length; i++) {
                    if (alts[i].length() == 0
                            || REF.charAt(0) != alts[i].byteAt(0)) {
                        same = false;
                        break;
                    }
                }
                if (same) {
                    pos++;
                    deletedReference = String.valueOf(REF.charAt(0));
                    REF = REF.substring(1);

                    //System.err.println(varId + " before :" + deletedReference);

                    for (int i = 0; i < alts.length; i++) {
                        alts[i].setSeq(new FlexSeq(alts[i].getSeq().substring(1)));
                    }
                } else {
                    break;
                }
            }

            // TODO this needs to be done
            // but if we want to preserve the original VCF record, then this
            // needs modification
            String clippedSequence = "";
            if (REF.length() > 0) {
                int referenceAlleleLength = REF.length();

                int minClipLength = Integer.MAX_VALUE;
                for (int i = 0; i < alts.length; i++) {
                    int alternativeAlleleLength = alts[i].length();

                    //what does clipLength represent?
                    int clipLength = 0;
                    for (int j = 0; j < Math.min(alternativeAlleleLength, referenceAlleleLength); j++) {
                        //this is based on the assumption that all characters are ASCII characters
                        if (REF.charAt(referenceAlleleLength - j - 1) != alts[i].byteAt(alternativeAlleleLength - j - 1)) {
                            clipLength = j;
                            break;
                        }
                        clipLength = j + 1;
                    }

                    minClipLength = Math.min(clipLength, minClipLength);
                }

                /*
                apparently this code block is part of normalization.
                is this working properly, though? e.g. it converts
                CGTG,CG => GT,""

                what this is doing is clip off tailing characters shared by both
                 REF and ALT.

                 make sure length after subtracting clip is nonnegative
                 */
                if (minClipLength > 0) {
                    clippedSequence = REF.substring(referenceAlleleLength - minClipLength, referenceAlleleLength);
                    REF = REF.substring(0, Math.max(0, referenceAlleleLength - minClipLength));
                    for (int i = 0; i < alts.length; i++) {
                        if (!clippedSequence.equals(new String(alts[i].getSeq().substring(alts[i].getSeq().length() - minClipLength, alts[i].getSeq().length())))) {
                            log.warn("Right clipping is initiated, but the clipped sequences are different for REF, ALT.");
                            return null;
                        }
                        alts[i].setSeq(new FlexSeq(alts[i].getSeq().substring(0,
                                Math.max(0, alts[i].length() - minClipLength))));
                    }
                }
            }

            byte[] refs = new byte[REF.length()];

            for (int i = 0; i < REF.length(); i++) {
                refs[i] = (byte) REF.charAt(i);
            }


            /*return new Variant(chr, pos, refs.length, refs, alts,
                    genotypeArray, isGenotypePhased, variantId, FILTER, deletedReference, random);
                    */
            return new Variant.Builder().chr(chr).pos(pos).referenceAlleleLength(refs.length).
                    ref(refs).alts(alts).phase(genotypeArray).isPhased(isGenotypePhased).
                    varId(variantId).filter(FILTER).qual(QUAL).info(infoString).refDeleted(deletedReference).
                    randomNumberGenerator(random).clippedSequence(clippedSequence).build();
        } else {
          // breakend
          log.warn("breakend is not handled directly now.");
          return null;
      }
    }

    /**
     * extract header from a VCF
     * must be run before parseLine, otherwise may return nothing
     */
    public String extractHeader() {
        StringBuilder stringBuilder = new StringBuilder();
        while (line != null && line.startsWith("#")) {
            stringBuilder.append(line);
            stringBuilder.append("\n");
            readLine();
        }

        //readLine automatically close filehandle when no more lines
        if (line != null) {
            try {
                this.bufferedReader.close(); //right now use this naive method to prevent more reading
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        return stringBuilder.toString();
    }

    public Variant parseLine() {
        /*
        TODO: line reading should be handled in a loop calling
         */
        String line = this.line;
        readLine();

        if (line == null || line.length() == 0) {
            log.info("blank line");
            return null;
        }

        Variant variant = null;

        try {
            variant = processLine(line);
        } catch (Exception e) {
            //TODO: right now just be lazy, die on any error
            log.fatal(e.getMessage() + "\n" + line);
            e.printStackTrace();
            System.exit(255);
        }

        if (variant == null && !line.startsWith("#")) {
            log.warn("Returned null variant for line " + line);
        }

        return variant;
    }

    /**
     * convert ALT field to Alt[] array
     * AT,ATG => ..
     * <DUP>,<DUP> => ..
     *
     * @param ALT
     * @return
     */
    public Alt[] string2Alt(String ALT){
        String[] altsString = StringUtilities.fastSplit(ALT, ",");
        Alt[] alts = new Alt[altsString.length];
        for (int i = 0; i < alts.length; i++) {
            alts[i] = Alt.altFactory(altsString[i]);
        }
        return alts;
    }

}
