package com.bina.varsim.util;

/**
 * Created by guoy28 on 6/26/17.
 */
public class VarSimMath {
  public static int max(Integer a, Integer b) {
      if (a == null && b == null) {
        throw new IllegalArgumentException("at least one Integer must be non-null");
      }
      if (a != null && b != null) {
        return Math.max(a, b);
      }
      return a == null ? b : a;
  }
}
