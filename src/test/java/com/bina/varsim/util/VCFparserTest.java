package com.bina.varsim.util;

import com.bina.varsim.types.variant.Variant;
import org.junit.Test;

import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 10/5/16.
 * some simple tests for VCF parsing
 */
public class VCFparserTest {
    @Test
    public void deletionParsingTest() {
        VCFparser runner = new VCFparser();
        Variant v = runner.process_line("12\t29557989\t.\tACAAAAGAAATGATCATGTTTGTAGGT\tAAAAAGAAATGATCATGTTTGTAGGT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        /*
        ref: A CAAAAGAAATGATCATGTTTGTAGGT
        alt: A A AAAGAAATGATCATGTTTGTAGGT
        */
        assertTrue(v.toString().equals("12\t29557989\t.\tAC\tA\t.\tPASS\tSVLEN=-1\tGT\t1|1"));
        assertTrue(v.isPhased());
        assertTrue(v.getRef().length == 1); //[C]
        assertTrue(v.getAlt(1).length() == 0); //""
    }
    @Test
    public void insertionParsingTest() {
        VCFparser runner = new VCFparser();
        Variant v = runner.process_line("12\t29557990\t.\tCTTT\tCGTTTT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        assertTrue(v.toString().equals("12\t29557990\t.\tC\tCGT\t.\tPASS\tSVLEN=2\tGT\t1|1"));
        assertTrue(v.isPhased());
        assertTrue(v.getRef().length == 0); //[]
        assertTrue(v.getAlt(1).length() == 2); //"GT"
    }
    @Test
    public void homopolymerInsertionParsingTest() {
        VCFparser runner = new VCFparser();
        Variant v = runner.process_line("12\t29557989\t.\tACT\tAAAACT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        assertTrue(v.toString().equals("12\t29557989\t.\tA\tAAAA\t.\tPASS\tSVLEN=3\tGT\t1|1"));
        assertTrue(v.isPhased());
        assertTrue(v.getRef().length == 0); //[]
        assertTrue(v.getAlt(1).length() == 3); //"GT"
    }
    @Test(expected=IllegalArgumentException.class)
    public void tandemDuplicationParsingTest() {
        VCFparser runner = new VCFparser();
        //genotype does not agree with copy number
        Variant v = runner.process_line("15\t85825565\tnssv534459\tT\t<DUP:TANDEM>\t.\tPASS\tSVTYPE=DUP;SVLEN=284016\tGT:CN\t0|1:2|2");
    }
}
