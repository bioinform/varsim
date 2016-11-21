package com.bina.varsim.types.variant.alt;

import com.bina.varsim.types.ChrString;
import com.bina.varsim.types.FlexSeq;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by guoy28 on 11/18/16.
 * class for ALT field in VCF format
 */
public final class Alt {
  //TODO convert this class to immutable
  private SymbolicAllele symbolicAllele;
  private Breakend breakend;
  private FlexSeq seq;

  /**
   * there are 3 possible ways to construct an ALT
   * 1) breakend
   * 2) symbolic allele
   * 3) alt sequence
   */
  public static Alt altFactory(String alt) {
    if (alt.startsWith("<")) {
      return new Alt(SymbolicAllele.symbolicAlleleFactory(alt));
    } else if (alt.indexOf("[") >= 0) {
      return new Alt(Breakend.breakendFactory(alt));
    } else {
      return new Alt(new FlexSeq(alt.getBytes()));
    }
  }

  public Alt() {}
  public Alt(SymbolicAllele s) {
    this.symbolicAllele = s;
  }
  public Alt(FlexSeq s) {
    this.seq = s;
  }
  public Alt(Breakend b) {
    this.breakend = b;
  }

  public Alt copy() {
    Alt a = new Alt();
    a.setBreakend(this.getBreakend()); //immutable
    a.setSeq(new FlexSeq(this.getSeq())); //mutable
    a.setSymbolicAllele(this.getSymbolicAllele()); //immutable
    return a;
  }

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
    TRA, //temporarily added for debugging
    BND; //breakend
    public enum SVSubtype {
      TANDEM, //consecutive events
      ME, //mobile element
      TRA, //translocation (cut-and-paste)
      IPS; //interspersed
    }
  }

  /**
   * immutable class for symbolic allele
   */
  public static final class SymbolicAllele{
    //TODO consider caching common symoblic alleles (e.g. <DEL>)
    /**
     * looks like <DUP:TANDEM>
     */
    private static Pattern r = Pattern.compile("<([^<:>]+):?([^<:>]*)>");
    private final SVType major;
    private final SVType.SVSubtype minor;

    public SVType getMajor() {
      return major;
    }

    public SVType.SVSubtype getMinor() {
      return minor;
    }

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
    public static SymbolicAllele symbolicAlleleFactory(String alt) {
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

  /**
   * immutable class for breakend
   */
  public static final class Breakend{
    /**
     * looks like ATGC[1:99[
     */
    private static Pattern r = Pattern.compile("(\\w+)([\\[\\]])([^:]+):(\\d+)([\\[\\]])|([\\[\\]])([^:]+):(\\d+)([\\[\\]])(\\w+)");
                                                //1     2       3       4     5           6         7       8        9     10
    private final byte[] seq;
    private final ChrString chr;
    private final long pos;
    private final boolean left; //ATGC[[, left=true, [[ATGC, left = false
    private final boolean forward; //[1:99[ forward = true, ]1:99] forward = false

    public Breakend(byte[] seq, ChrString chr, long pos, boolean left, boolean forward) {
      this.seq = seq;
      this.chr = chr;
      this.pos = pos;
      this.left = left;
      this.forward = forward;
    }
    public static Breakend breakendFactory(String alt) {
      Matcher m = r.matcher(alt);
      if (m.find()) {
        if (m.group(1) != null) {
          return new Breakend(m.group(1).getBytes(), new ChrString(m.group(3)), Long.parseLong(m.group(4)), true, m.group(2).indexOf("[") >= 0 );
        } else {
          return new Breakend(m.group(10).getBytes(), new ChrString(m.group(7)), Long.parseLong(m.group(8)), false, m.group(6).indexOf("[") >= 0 );
        }
      }
      return null; //failed to parse as a symbolic allele
    }

    @Override
    public String toString() {
      if (left) {
        return new String(seq) + (forward?("[" + chr + ":" + pos + "["):("]" + chr + ":" + pos + "]"));
      } else {
        return  forward? ("[" + chr + ":" + pos + "["):("]" + chr + ":" + pos + "]") + new String(seq);
      }
    }

    @Override
    public int hashCode() {
      return this.toString().hashCode();
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) return true;
      if (!(o instanceof Breakend)) {
        return false;
      }
      Breakend that = (Breakend) o;
      if (!this.toString().equals(that.toString())) return false;
      return true;
    }

    /**
     * return a clone copy of underlying sequence
     * the class is immutable so any change to
     * the sequence (array) is irrelevant, has no
     * effect on the object
     *
     * @return
     */
    public byte[] getSeq() { return seq.clone(); }

    public ChrString getChr() {
      return chr;
    }

    public long getPos() {
      return pos;
    }

    public boolean isLeft() {
      return left;
    }

    public boolean isForward() {
      return forward;
    }
  }
}
