package com.bina.varsim.util;

import com.bina.varsim.types.ChrString;
import com.google.common.collect.ImmutableList;

import java.util.List;
import java.util.StringJoiner;

/**
 * Created by guoy28 on 3/2/18.
 */
public class VCFWriter {
  /**
   * generate VCF file header
   *
   * a little explanation about translocation-related meta-info lines:
   * each locus involved in translocations is modeled independently as
   * a cut-paste event with slight variations for different types of
   * translocations. Basically, a region at locus A (the sink) will be
   * cut (deleted), and a region at locus B (the source) will be placed
   * at the sink. The placement may be: complete transfer; complete transfer
   * with inversion; no transfer (one-way or unbalanced translocation).
   *
   * @param referenceFileName reference file name
   * @param sampleNames list of sample names
   * @return
   */
  public static String generateVCFHeader(final String referenceFileName, final ImmutableList<String> sampleNames) {
    StringBuilder VCFHeader = new StringBuilder();
    VCFHeader.append("##fileformat=VCFv4.3\n" +
            "##reference=" + referenceFileName + "\n" +
                /*
                SVLEN is for alternative allele in truth VCF
                 */
            "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n" +
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n" +
                /*if POS2<=END2, then another sequence is inserted at positive strand
                if POS2>=END2, then reversed sequence is inserted at negative strand (insert with inversion)
                 */
            "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"1-based start position of source sequence\">\n" +
            "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"1-based end position of source sequence\">\n" +
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"1-based end position of variant\">\n" +
            "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome of source sequence\">\n" +
            "##INFO=<ID=ISINV,Number=1,Type=Flag,Description=\"whether a duplication is inverted\">\n" +
            "##INFO=<ID=TRAID,Number=1,Type=String,Description=\"translocation ID\">\n" +
            "##INFO=<ID=IMPRECISE_LENGTH,Number=1,Type=Flag,Description=\"SVLEN is imprecise\">\n" +
            "##INFO=<ID=VARIANT_OVERALL_TYPE,Number=1,Type=String,Description=\"Overall variant type\">\n" +
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                /*CN is defined as Integer in VCF4.1,4.3, making it impossible to specify multiple CN values
                here we changed it to String to allow such behavior.
                 */
            //TODO: this will be changed back to VCF4.3 format later.
            "##FORMAT=<ID=CN,Number=1,Type=String,Description=\"Copy number genotype.\">\n" +
            "##ALT=<ID=DEL,Description=\"Deletion\">\n" +
            "##ALT=<ID=DEL:TRA,Description=\"Deletion in translocation\">\n" +
            "##ALT=<ID=DUP,Description=\"Duplication\">\n" +
            "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n" +
            "##ALT=<ID=DUP:ISP,Description=\"Interspersed duplication\">\n" +
            "##ALT=<ID=DUP:TRA,Description=\"Duplication in translocation\">\n" +
            "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n" +
            "##ALT=<ID=INV,Description=\"Inversion\">\n");
    if (referenceFileName != null) {
      SimpleReference ref = new SimpleReference(referenceFileName);
      for (ChrString contig : ref.keySet()) {
        VCFHeader.append("##contig=<ID=" + contig + ",length=" + ref.getRefLen(contig) + ">\n");
      }
    }
    StringJoiner joiner = new StringJoiner("\t");
    for (String id : sampleNames) {
      joiner.add(id);
    }
    VCFHeader.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" +
            (sampleNames.isEmpty() ? "" : "\t") + joiner.toString() + "\n");
    return VCFHeader.toString();
  }
}
