package com.bina.varsim.types;

import com.bina.varsim.types.variant.alt.Alt;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * Created by guoy28 on 11/18/16.
 */
public class AltTest {
  @Test
  public void breakendTest() {
    Alt.Breakend test = new Alt.Breakend("A".getBytes(),new ChrString("1") ,123,true, true);
    assertTrue(test.toString().equals("A[1:123["));
  }
  @Test (expected = Exception.class)
  public void symbolicAlleleTest() {
    assertTrue(new Alt.SymbolicAllele("DUP").toString().equals("<DUP>"));
    assertTrue(new Alt.SymbolicAllele("DUP","TANDEM").toString().equals("<DUP:TANDEM>"));
    assertTrue(Alt.SymbolicAllele.SymbolicAlleleFactory("GAT") == null);
    assertTrue(Alt.SymbolicAllele.SymbolicAlleleFactory("<DUP:TANDEM").toString().equals("<DUP:TANDEM>"));
    new Alt.SymbolicAllele("TX");
  }
}
