package com.bina.varsim.util;

import com.bina.varsim.types.variant.Variant;
import org.junit.Before;
import org.junit.Test;
import org.junit.Ignore;

import java.io.IOException;
import java.rmi.UnexpectedException;

import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 10/5/16.
 * some simple tests for VCF parsing
 */
public class VCFparserTest {
    private VCFparser parser;
    @Before
    public void setup() throws UnexpectedException{
        parser = new VCFparser();
        parser.processLine("##fileformat=VCFv4.1");
        parser.processLine("##reference=src/test/resources/DuplicationTest/oneDuplicationTest.fa");
        parser.processLine("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of variant\">");
        parser.processLine("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
        parser.processLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        parser.processLine("##ALT=<ID=DEL,Description=\"Deletion\">");
        parser.processLine("##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">");
        parser.processLine("##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">");
        parser.processLine("##ALT=<ID=DUP,Description=\"Duplication\">");
        parser.processLine("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">");
        parser.processLine("##ALT=<ID=INS,Description=\"Insertion of novel sequence\">");
        parser.processLine("##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">");
        parser.processLine("##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">");
        parser.processLine("##ALT=<ID=INV,Description=\"Inversion\">");
        parser.processLine("##ALT=<ID=CNV,Description=\"Copy number variable region\">");
        parser.processLine("##ALT=<ID=ITX,Description=\"Intra-chromosomal translocation\">");
        parser.processLine("##ALT=<ID=CTX,Description=\"Inter-chromosomal translocation\">");
        parser.processLine("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	test");
    }
    @Test
    public void deletionParsingTest() throws UnexpectedException {
        Variant v = parser.processLine("12\t29557989\t.\tACAAAAGAAATGATCATGTTTGTAGGT\tAAAAAGAAATGATCATGTTTGTAGGT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        /*
        ref: A CAAAAGAAATGATCATGTTTGTAGGT
        alt: A A AAAGAAATGATCATGTTTGTAGGT
        */
        assertTrue(v.toString().equals("12\t29557989\t.\tAC\tA\t.\tPASS\tSVTYPE=DEL;SVLEN=-1\tGT\t1|1"));
        assertTrue(v.isPhased());
        assertTrue(v.getReference().length == 1); //[C]
        assertTrue(v.getAlt(1).length() == 0); //""
    }
    @Test
    public void insertionParsingTest() throws UnexpectedException {
        Variant v = parser.processLine("12\t29557990\t.\tCTTT\tCGTTTT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        assertTrue(v.toString().equals("12\t29557990\t.\tC\tCGT\t.\tPASS\tSVLEN=2\tGT\t1|1"));
        assertTrue(v.isPhased());
        assertTrue(v.getReference().length == 0); //[]
        assertTrue(v.getAlt(1).length() == 2); //"GT"
    }
    @Test
    public void homopolymerInsertionParsingTest() throws UnexpectedException {
        Variant v = parser.processLine("12\t29557989\t.\tACT\tAAAACT\t.\tPASS\tSVLEN=-26\tGT\t1|1");
        assertTrue(v.toString().equals("12\t29557989\t.\tA\tAAAA\t.\tPASS\tSVLEN=3\tGT\t1|1"));
        assertTrue(v.isPhased());
        assertTrue(v.getReference().length == 0); //[]
        assertTrue(v.getAlt(1).length() == 3); //"GT"
    }
    @Test(expected = IllegalArgumentException.class)
    public void SymbolicAlleleSVLenTest() throws UnexpectedException {
        parser.processLine("1	3	.	T	<DUP:TANDEM>	.	PASS	SVTYPE=DUP;SVLEN=4,4	GT:CN	1|2:2|3");
    }
    @Test(expected = IllegalArgumentException.class)
    public void SymbolicAlleleMixWithNonSymbolicAlleleTest() throws IOException {
        parser.processLine("1	3	.	T	<DUP:TANDEM>,TTTT	.	PASS	SVTYPE=DUP;SVLEN=4,4	GT:CN	1|2:2|3");
    }
    @Test
    public void addExtraBaseForInsertion() throws UnexpectedException {
        Variant v = parser.processLine("1	3	.	T	TTTT	.	PASS	SVLEN=3	GT:CN	0|1:2|3");
        assertTrue(v.toString().equals("1	3	.	T	TTTT	.	PASS	SVLEN=3	GT	0|1"));
    }
    @Test
    public void addExtraBaseForDeletion() throws UnexpectedException {
        Variant v = parser.processLine("1	3	.	TTTT	T	.	PASS	SVLEN=-3	GT:CN	0|1:2|3");
        assertTrue(v.toString().equals("1	3	.	TTTT	T	.	PASS	SVTYPE=DEL;SVLEN=-3	GT	0|1"));
    }
    @Test
    public void vcfinfoTestDefaultBoolean() throws UnexpectedException {
        //for unrecognized flags, return boolen type
        Variant	v	=	parser.processLine("1	111	rs770821123	C	A	.	.	RS=770821123;RSPOS=10000111;VP=0x050000080005000002000100;dbSNPBuildID=144;SAO=0;SSR=0;WGT=1;VC=SNV;INT;ASP");
        assertTrue(v.toString().equals("1	111	rs770821123	C	A	.	.	SVLEN=1	GT	1|1"));
    }
}
