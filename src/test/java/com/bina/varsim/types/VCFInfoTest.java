package com.bina.varsim.types;

import org.junit.Test;

import java.util.List;

import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 10/7/16.
 */
public class VCFInfoTest {
    @Test
    public void VCFInfoParsingTest() {
        try {
            VCFInfo test = new VCFInfo("SVTYPE=TRA;TRASUBTYPE=CHIVALRY,SELFISHNESS;SVLEN=100,200;POS2=155000679,156000999;END2=155110000,155100999;CHR2=1,1");
            List<String> svtypeList = (List<String>) test.getValue("SVTYPE");
            assertTrue(svtypeList.get(0).equals("TRA"));
            List<Integer> posList = (List<Integer>) test.getValue("POS2");
            assertTrue(posList.get(0) == 155000679 && posList.get(1) == 156000999);
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
}
