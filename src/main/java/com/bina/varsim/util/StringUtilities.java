package com.bina.varsim.util;

import com.bina.varsim.types.ChrString;

import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import java.util.StringTokenizer;

/**
 * Created by guoy28 on 10/10/16.
 */
public class StringUtilities {
    /**
     * convert integer array to delimited string
     * @param a
     * @return
     */
    public static String concatenateArray(final int[] a, final CharSequence delimiter) {
        if (a == null) {
            return "";
        }
        StringJoiner joiner = new StringJoiner(delimiter.toString());
        for (int i = 0; i < a.length; i++) {
            joiner.add(Integer.toString(a[i]));
        }
        return joiner.toString();
    }

    /**
     * convert String array to delimitted string
     * @param a
     * @param delimiter
     * @return
     */
    public static String concatenateArray(final String[] a, final CharSequence delimiter) {
        if (a == null) {
            return "";
        }
        StringJoiner joiner = new StringJoiner(delimiter.toString());
        for (int i = 0; i < a.length; i++) {
            joiner.add(a[i]);
        }
        return joiner.toString();
    }

    /**
     * convert ChrString array to delimited string
     * @param a
     * @param delimiter
     * @return
     */
    public static String concatenateArray(final ChrString[] a, final CharSequence delimiter) {
        if (a == null) {
            return "";
        }
        StringJoiner joiner = new StringJoiner(delimiter.toString());
        for (int i = 0; i < a.length; i++) {
            joiner.add(a[i].toString());
        }
        return joiner.toString();
    }

    /**
     * split string into String array
     */
    public static String[] fastSplit(final String s, final String delim) {
        StringTokenizer st = new StringTokenizer(s, delim);
        String[] result = new String[st.countTokens()];
        for (int i = 0; i < result.length; i++) {
            result[i] = st.nextToken();
        }
        return result;
    }

    /**
     * check if a string is a non-negative integer
     * empty string is considered illegal non-negative integer
     * reference: https://stackoverflow.com/questions/10575624/java-string-see-if-a-string-contains-only-numbers-and-not-letters
     * @param s
     * @return
     */
    public static boolean isNonNegativeInteger(final String s) {
      if (s == null || s.length() == 0) {
          return false;
      }
      for (int i = 0; i < s.length(); i++) {
          char c = s.charAt(i);
          if (c > '9' || c < '0') {
              return false;
          }
      }
      return true;
    }

    /**
     * parse an Integer, there is no sanity check for speed reasons!
     * @param s
     * @return
     */
    public static int parseInt(final String s) {
        int n = 0;
        int sign = -1;
        final int length = s.length();
        final char firstChar = s.charAt(0);
        if (firstChar == '-') {
            sign = 1;
        } else {
            n = '0' - firstChar;
        }
        for (int i = 1; i < length; i++) {
            n = n * 10 + '0' - s.charAt(i);
        }
        n = sign * n;
        return n;
    }
}
