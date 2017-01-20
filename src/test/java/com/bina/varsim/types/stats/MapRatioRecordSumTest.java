package com.bina.varsim.types.stats;

import org.junit.Test;

import java.util.EnumSet;
import java.util.Set;
import com.bina.varsim.types.stats.MapRatioRecordSum.EventTypesForStats;

import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 1/18/17.
 */
public class MapRatioRecordSumTest {
  @Test
  public void emptyRecordTest() {
    MapRatioRecordSum x = new MapRatioRecordSum();
    assertTrue(x.toString().equals("---------\n" + "---------\n"));
  }
  @Test
  public void oneTPTest() {
    Set<EventTypesForStats> features = EnumSet.noneOf(EventTypesForStats.class);
    MapRatioRecordSum x = new MapRatioRecordSum();
            features.add(EventTypesForStats.All);
    x.incStat(features, 1, StatsNamespace.TP);
    assertTrue(x.toString().equals("0,Infinity,0.0000,1,0,0,NaN\n" +
            "1,Infinity,0.0000,1,0,0,NaN\n" +
            "---------\n" +
            "All\n" +
            "0,Infinity,0.0000,1,0,0,NaN\n" +
            "1,Infinity,0.0000,1,0,0,NaN\n" +
            "\n" +
            "---------\n"));
    x.incStat(features, -1, StatsNamespace.TP);
    features.add(EventTypesForStats.Sequence);
    x.incStat(features, 1, StatsNamespace.TP);
    assertTrue(x.toString().equals("0,Infinity,0.0000,2,0,0,NaN\n" +
            "1,Infinity,0.0000,2,0,0,NaN\n" +
            "---------\n" +
            "All\n" +
            "0,Infinity,0.0000,2,0,0,NaN\n" +
            "1,Infinity,0.0000,2,0,0,NaN\n" +
            "\n" +
            "Sequence\n" +
            "0,Infinity,0.0000,1,0,0,NaN\n" +
            "1,Infinity,0.0000,1,0,0,NaN\n" +
            "\n" +
            "---------\n"));
  }
}
