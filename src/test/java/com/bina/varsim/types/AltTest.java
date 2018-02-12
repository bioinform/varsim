package com.bina.varsim.types;

import com.bina.varsim.types.variant.alt.Alt;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * Created by guoy28 on 11/18/16.
 */
public class AltTest {
  @Test
  public void breakendTest() {
    Alt.Breakend test = new Alt.Breakend("A".getBytes(),new ChrString("1") ,123,true, true);
    assertTrue(test.toString().equals("A[1:123["));
    assertTrue(Alt.Breakend.breakendFactory("A[1:1[").toString().equals("A[1:1["));
    assertFalse(Alt.Breakend.breakendFactory("]1:1]A").toString().equals("A[1:1["));
    assertTrue(Alt.Breakend.breakendFactory("]1:1999999]A").toString().equals("]1:1999999]A"));
  }
  @Test (expected = Exception.class)
  public void symbolicAlleleTest() {
    assertTrue(new Alt.SymbolicAllele("DUP").toString().equals("<DUP>"));
    assertTrue(new Alt.SymbolicAllele("DUP","TANDEM").toString().equals("<DUP:TANDEM>"));
    assertTrue(Alt.SymbolicAllele.symbolicAlleleFactory("GAT") == null);
    assertTrue(Alt.SymbolicAllele.symbolicAlleleFactory("<DUP:TANDEM").toString().equals("<DUP:TANDEM>"));
    new Alt.SymbolicAllele("TX");
  }
  @Test
  public void altTest() {
    assertTrue(Alt.altFactory("<DUP:TANDEM>").getSymbolicAllele().toString().equals("<DUP:TANDEM>"));
    assertTrue(Alt.altFactory("ANTCCC[X:123[").getBreakend().toString().equals("ANTCCC[X:123["));
    assertTrue(Alt.altFactory("ATCCG").getSeq().toString().equals("ATCCG"));
  }
  @Test(expected = IllegalArgumentException.class)
  public void illegalCharTest() {
    Alt.altFactory("ATCCGD");
  }
}
