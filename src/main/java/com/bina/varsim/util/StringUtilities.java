package com.bina.varsim.util;

import com.bina.varsim.types.ChrString;

import java.util.ArrayList;
import java.util.List;

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
     * partial implementation of Java 8 StringJoiner
     */
    public static class StringJoiner {
        private String delimiter;
        private List<String> fields;

      public StringJoiner(String delimiter) {
          this.delimiter = delimiter;
          this.fields = new ArrayList<String>();
      }

      public void add(String s) {
          fields.add(s);
      }

      @Override
      public String toString() {
        StringBuilder sb = new StringBuilder();
          for (int i = 0; i < fields.size() - 1; i++) {
            sb.append(fields.get(i));
            sb.append(delimiter);
          }
          if (fields.size() > 0)
            sb.append(fields.get(fields.size() - 1));
        return sb.toString();
      }
    }
}
