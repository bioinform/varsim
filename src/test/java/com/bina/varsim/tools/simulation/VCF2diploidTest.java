package com.bina.varsim.tools.simulation;

import org.apache.commons.io.FileUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
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
    @Test
    public void MixedTest() throws IOException {
        File wd = tmpFolder.newFolder("MixedTest");
        String reference = "src/test/resources/MixedTest/oneMixedTest.fa";
        String vcf = "src/test/resources/MixedTest/oneMixedTest.vcf";
        String map = "src/test/resources/MixedTest/oneMixedTest.map";
        String maternalReference = "src/test/resources/MixedTest/1_test_maternal.fa";
        String paternalReference = "src/test/resources/MixedTest/1_test_paternal.fa";

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
    public void balancedReciprocalTranslocationIntrachromosomalTest() throws IOException {
        File wd = tmpFolder.newFolder("BalancedReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal.vcf";
        String expectedVCF1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal.expected.2.vcf";
        String map = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal.map";
        String maternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedReciprocalTranslocationIntrachromosomalWithInversionTest() throws IOException {
        File wd = tmpFolder.newFolder("BalancedReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion.vcf";
        String expectedVCF1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion.expected.2.vcf";
        String map = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion.map";
        String maternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedReciprocalTranslocationInterchromosomalTest() throws IOException {
        File wd = tmpFolder.newFolder("BalancedReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal.vcf";
        String expectedVCF1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal.expected.2.vcf";
        String map = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal.map";
        String maternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedReciprocalTranslocationInterchromosomalWithInversionTest() throws IOException {
        File wd = tmpFolder.newFolder("BalancedReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion.vcf";
        String expectedVCF1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion.expected.2.vcf";
        String map = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion.map";
        String maternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedNonReciprocalTranslocationInterchromosomal() throws IOException {
        File wd = tmpFolder.newFolder("BalancedNonReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedNonReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal.vcf";
        String expectedVCF1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal.expected.2.vcf";
        String map = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal.map";
        String maternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedNonReciprocalTranslocationInterchromosomalWithInversionTest() throws IOException {
        File wd = tmpFolder.newFolder("BalancedNonReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedNonReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion.vcf";
        String expectedVCF1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion.expected.2.vcf";
        String map = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion.map";
        String maternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedNonReciprocalTranslocationInterHomologousChromosomal() throws IOException {
        File wd = tmpFolder.newFolder("BalancedNonReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedNonReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal.vcf";
        String expectedVCF1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal.expected.2.vcf";
        String map = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal.map";
        String maternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedNonReciprocalTranslocationIntrachromosomal() throws IOException {
        File wd = tmpFolder.newFolder("BalancedNonReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedNonReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal.vcf";
        String expectedVCF1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal.expected.2.vcf";
        String map = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal.map";
        String maternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Test
    public void balancedNonReciprocalTranslocationIntrachromosomalWithInversion() throws IOException {
        File wd = tmpFolder.newFolder("BalancedNonReciprocalTranslocationTest");
        String reference = "src/test/resources/BalancedNonReciprocalTranslocationTest/twelveTranslocationTest.fa";
        String vcf = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion.vcf";
        String expectedVCF1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion.expected.1.vcf";
        String expectedVCF2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion.expected.2.vcf";
        String map = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion.map";
        String maternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion_1_maternal.fa";
        String paternalReference1 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion_1_paternal.fa";
        String maternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion_2_maternal.fa";
        String paternalReference2 = "src/test/resources/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion_2_paternal.fa";

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputVCF2Path = Paths.get(wd.getCanonicalPath(), "2_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");
        Path outputMaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_maternal.fa");
        Path outputPaternalReference2Path = Paths.get(wd.getCanonicalPath(), "2_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
}
