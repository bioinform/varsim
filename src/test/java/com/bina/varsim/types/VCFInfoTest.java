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
            VCFInfo test = new VCFInfo("SVTYPE=TRA;END=156001999;TRASUBTYPE=ACCEPT,REJECT;SVLEN=100,200;POS2=155000679,156000999;END2=155110000,155100999;CHR2=1,1");
            String[] svtypeList = test.getValue("SVTYPE", String[].class);
            assertTrue(svtypeList[0].equals("TRA"));
            int[] posList = test.getValue("POS2", int[].class);
            assertTrue(posList[0] == 155000679 && posList[1] == 156000999);
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
}
