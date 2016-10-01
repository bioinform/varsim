package com.bina.varsim.tools.simulation;

import org.apache.commons.io.FileUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.fail;

/**
 * Created by guoy28 on 9/29/16.
 * test various methods in VCF2diploid class
 */
public class VCF2diploidTest {
    private int seed = 11;
    @Rule
    public TemporaryFolder tmpFolder = new TemporaryFolder();

    /*
    assume maven will automatically go to root directory of the project.
    based on http://stackoverflow.com/questions/28673651/how-to-get-the-path-of-src-test-resources-directory-in-junit
    it might not be a good idea to use ClassLoader (they are not java.io.File objects in a jar file).
     */
    @Test
    public void SNPTest() throws IOException {
        File wd = tmpFolder.newFolder("SNPTest");
        String reference = "src/test/resources/SNPTest/oneSNPTest.fa";
        String vcf = "src/test/resources/SNPTest/oneSNPTest.vcf";
        String map = "src/test/resources/SNPTest/oneSNPTest.map";
        String maternalReference = "src/test/resources/SNPTest/1_test_maternal.fa";
        String paternalReference = "src/test/resources/SNPTest/1_test_paternal.fa";

        Path outputVCFPath = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputMaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCFPath.toFile(), new File(vcf)));
        assertTrue(FileUtils.contentEquals(outputMaternalReferencePath.toFile(), new File(maternalReference)));
        assertTrue(FileUtils.contentEquals(outputPaternalReferencePath.toFile(), new File(paternalReference)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void MNPTest() throws IOException {
        File wd = tmpFolder.newFolder("MNPTest");
        String reference = "src/test/resources/MNPTest/oneMNPTest.fa";
        String vcf = "src/test/resources/MNPTest/oneMNPTest.vcf";
        String expectedVCF = "src/test/resources/MNPTest/oneMNPTestExpected.vcf";
        String map = "src/test/resources/MNPTest/oneMNPTest.map";
        String maternalReference = "src/test/resources/MNPTest/1_test_maternal.fa";
        String paternalReference = "src/test/resources/MNPTest/1_test_paternal.fa";

        Path outputVCFPath = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputMaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCFPath.toFile(), new File(expectedVCF)));
        assertTrue(FileUtils.contentEquals(outputMaternalReferencePath.toFile(), new File(maternalReference)));
        assertTrue(FileUtils.contentEquals(outputPaternalReferencePath.toFile(), new File(paternalReference)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void InsertionTest() throws IOException {
        File wd = tmpFolder.newFolder("InsertionTest");
        String reference = "src/test/resources/InsertionTest/oneInsertionTest.fa";
        String vcf = "src/test/resources/InsertionTest/oneInsertionTest.vcf";
        String map = "src/test/resources/InsertionTest/oneInsertionTest.map";
        String maternalReference = "src/test/resources/InsertionTest/1_test_maternal.fa";
        String paternalReference = "src/test/resources/InsertionTest/1_test_paternal.fa";

        Path outputVCFPath = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputMaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCFPath.toFile(), new File(vcf)));
        assertTrue(FileUtils.contentEquals(outputMaternalReferencePath.toFile(), new File(maternalReference)));
        assertTrue(FileUtils.contentEquals(outputPaternalReferencePath.toFile(), new File(paternalReference)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void DeletionTest() throws IOException {
        File wd = tmpFolder.newFolder("DeletionTest");
        String reference = "src/test/resources/DeletionTest/oneDeletionTest.fa";
        String vcf = "src/test/resources/DeletionTest/oneDeletionTest.vcf";
        String map = "src/test/resources/DeletionTest/oneDeletionTest.map";
        String maternalReference = "src/test/resources/DeletionTest/1_test_maternal.fa";
        String paternalReference = "src/test/resources/DeletionTest/1_test_paternal.fa";

        Path outputVCFPath = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputMaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCFPath.toFile(), new File(vcf)));
        assertTrue(FileUtils.contentEquals(outputMaternalReferencePath.toFile(), new File(maternalReference)));
        assertTrue(FileUtils.contentEquals(outputPaternalReferencePath.toFile(), new File(paternalReference)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void InversionTest() throws IOException {
        File wd = tmpFolder.newFolder("InversionTest");
        String reference = "src/test/resources/InversionTest/oneInversionTest.fa";
        String vcf = "src/test/resources/InversionTest/oneInversionTest.vcf";
        String map = "src/test/resources/InversionTest/oneInversionTest.map";
        String maternalReference = "src/test/resources/InversionTest/1_test_maternal.fa";
        String paternalReference = "src/test/resources/InversionTest/1_test_paternal.fa";

        Path outputVCFPath = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputMaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCFPath.toFile(), new File(vcf)));
        assertTrue(FileUtils.contentEquals(outputMaternalReferencePath.toFile(), new File(maternalReference)));
        assertTrue(FileUtils.contentEquals(outputPaternalReferencePath.toFile(), new File(paternalReference)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void DuplicationTest() throws IOException {
        File wd = tmpFolder.newFolder("DuplicationTest");
        String reference = "src/test/resources/DuplicationTest/oneDuplicationTest.fa";
        String vcf = "src/test/resources/DuplicationTest/oneDuplicationTest.vcf";
        String map = "src/test/resources/DuplicationTest/oneDuplicationTest.map";
        String maternalReference = "src/test/resources/DuplicationTest/1_test_maternal.fa";
        String paternalReference = "src/test/resources/DuplicationTest/1_test_paternal.fa";

        Path outputVCFPath = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputMaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReferencePath = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCFPath.toFile(), new File(vcf)));
        assertTrue(FileUtils.contentEquals(outputMaternalReferencePath.toFile(), new File(maternalReference)));
        assertTrue(FileUtils.contentEquals(outputPaternalReferencePath.toFile(), new File(paternalReference)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
}
