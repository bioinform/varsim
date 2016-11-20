package com.bina.varsim.util;

import org.junit.Test;

import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 11/20/16.
 */
public class StringUtilitiesTest {
  @Test
  public void stringJoinerTest() {
    StringUtilities.StringJoiner joiner = new StringUtilities.StringJoiner(" ");
    assertTrue(joiner.toString().equals(""));
    joiner.add("x");
    assertTrue(joiner.toString().equals("x"));
    joiner.add("y");
    assertTrue(joiner.toString().equals("x y"));
  }
}
