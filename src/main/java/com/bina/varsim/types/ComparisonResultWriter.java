package com.bina.varsim.types;

import com.google.common.collect.ImmutableList;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

/**
 * Created by guoy28 on 12/15/16.
 */
public enum ComparisonResultWriter {
  TP_WRITER("_TP.vcf", FileType.VCF),
  UNKNOWN_TP_WRITER("_unknown_TP.vcf", FileType.VCF),
  FP_WRITER("_FP.vcf", FileType.VCF),
  UNKNOWN_FP_WRITER("_unknown_FP.vcf", FileType.VCF),
  FN_WRITER("_FN.vcf", FileType.VCF),
  JSON_WRITER("_report.json", FileType.JSON);

  private static final String encoding = "UTF-8";
  private FileType fileType;
  private final String suffix;
  private ComparisonResultWriter(String suffix, FileType fileType) {
    this.suffix = suffix;
    this.fileType = fileType;
  }

  public PrintWriter getWriter(String prefix) throws FileNotFoundException, UnsupportedEncodingException{
    PrintWriter writer = new PrintWriter(prefix + suffix, encoding);
    return writer;
  }
}
