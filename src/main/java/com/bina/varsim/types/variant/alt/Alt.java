package com.bina.varsim.types.variant.alt;

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;

import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by guoy28 on 11/18/16.
 * class for ALT field in VCF format
 */
public class Alt {
  private SymbolicAllele symbolicAllele;
  private Breakend breakend;
  private FlexSeq seq;
  public Alt() {}

  public SymbolicAllele getSymbolicAllele() {
    return symbolicAllele;
  }

  public void setSymbolicAllele(SymbolicAllele symbolicAllele) {
    this.symbolicAllele = symbolicAllele;
  }

  public Breakend getBreakend() {
    return breakend;
  }

  public void setBreakend(Breakend breakend) {
    this.breakend = breakend;
  }

  public FlexSeq getSeq() {
    return seq;
  }

  public void setSeq(FlexSeq seq) {
    this.seq = seq;
  }

  public enum SVType {
    DEL, INS, DUP, INV, CNV,
    BND; //breakend
    public enum SVSubtype {
      TANDEM, //consecutive events
      ME, //mobile element
      TRA, //translocation (cut-and-paste)
      IPS; //interspersed
    }
  }
  public static class SymbolicAllele{
    /**
     * looks like <DUP:TANDEM>
     */
    private SVType major;
    private SVType.SVSubtype minor;

    @Override
    public String toString() {
      return "<" + major + (minor == null? ">" : (":" + minor + ">"));
    }

    public SymbolicAllele(String major, String minor) {
      this.major = SVType.valueOf(major);
      this.minor = SVType.SVSubtype.valueOf(minor);
    }
    public SymbolicAllele(String major) {
      this.major = SVType.valueOf(major);
      this.minor = null;
    }

    /**
     * parse <X:Y>
     * @param alt
     */
    public static SymbolicAllele SymbolicAlleleFactory(String alt) {
      Pattern r = Pattern.compile("<([^<:>]+):?([^<:>]*)>");
      Matcher m = r.matcher(alt);
      if (m.find()) {
        if( m.group(2).length() > 0) {
          return new SymbolicAllele(m.group(1), m.group(2));
        } else {
          return new SymbolicAllele(m.group(1));
        }
      }
      return null; //failed to parse as a symbolic allele
    }
  }
  public static class Breakend{
    /**
     * looks like ATGC[1:99[
     */
    private byte[] seq;
    private ChrString chr;
    private long pos;
    private boolean left; //ATGC[[, left=true, [[ATGC, left = false
    private boolean forward; //[1:99[ forward = true, ]1:99] forward = false
    public Breakend(byte[] seq, ChrString chr, long pos, boolean left, boolean forward) {
      this.seq = seq;
      this.chr = chr;
      this.pos = pos;
      this.left = left;
      this.forward = forward;
    }

    @Override
    public String toString() {
      if (left) {
        return new String(seq) + (forward?("[" + chr + ":" + pos + "["):("]" + chr + ":" + pos + "]"));
      } else {
        return  forward? ("[" + chr + ":" + pos + "["):("]" + chr + ":" + pos + "]") + new String(seq);
      }
    }
  }
}
