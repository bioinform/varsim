package com.bina.varsim.types;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

/**
 * Created by guoy28 on 12/15/16.
 */
public enum ComparisonResultWriter {
  TP_WRITER("_TP.vcf"),
  UNKNOWN_TP_WRITER("_unknown_TP.vcf"),
  FP_WRITER("_FP.vcf"),
  UNKNOWN_FP_WRITER("_unknown_FP.vcf"),
  FN_WRITER("_FN.vcf"),
  JSON_WRITER("_report.json");

  private static final String encoding = "UTF-8";
  private final String suffix;
  private ComparisonResultWriter(String suffix) {
    this.suffix = suffix;
  }
  public PrintWriter getWriter(String prefix) throws FileNotFoundException, UnsupportedEncodingException{
    return new PrintWriter(prefix + suffix, encoding);
  }
}
