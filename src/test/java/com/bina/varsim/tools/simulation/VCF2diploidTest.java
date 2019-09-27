package com.bina.varsim.tools.simulation;

import org.apache.commons.io.FileUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

import static com.bina.varsim.GlobalTestingOptions.updateVCF;
import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.fail;

/**
 * Created by guoy28 on 9/29/16.
 * test various methods in VCF2diploid class
 */
public class VCF2diploidTest {
    private int seed = 11;
    public void universalTestMethod(String directory) throws IOException {
        File wd = tmpFolder.newFolder("tmp");
        String reference = new File(directory, "reference.fa").toString();
        String vcf = new File(directory, "input.vcf").toString();
        String expectedVCF1 = new File(directory, "1.vcf").toString();
        String expectedVCF2 = new File(directory, "2.vcf").toString();
        String map = new File(directory, "expected.map").toString();
        String maternalReference1 = new File(directory, "maternal.1.fa").toString();
        String paternalReference1 = new File(directory, "paternal.1.fa").toString();
        String maternalReference2 = new File(directory, "maternal.2.fa").toString();
        String paternalReference2 = new File(directory, "paternal.2.fa").toString();

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
        if (updateVCF) {
            Files.copy(outputVCF1Path, Paths.get(expectedVCF1), StandardCopyOption.REPLACE_EXISTING);
            Files.copy(outputVCF2Path, Paths.get(expectedVCF2), StandardCopyOption.REPLACE_EXISTING);
        }
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputVCF2Path.toFile(), new File(expectedVCF2)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference2Path.toFile(), new File(maternalReference2)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference2Path.toFile(), new File(paternalReference2)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    public void universalTestMethod2(String directory) throws IOException {
        File wd = tmpFolder.newFolder("tmp");
        String reference = new File(directory, "reference.fa").toString();
        String vcf = new File(directory, "input.vcf").toString();
        String expectedVCF1 = new File(directory, "1.vcf").toString();
        String map = new File(directory, "expected.map").toString();
        String maternalReference1 = new File(directory, "maternal.1.fa").toString();
        String paternalReference1 = new File(directory, "paternal.1.fa").toString();

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "1_test.vcf");
        Path outputMaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_maternal.fa");
        Path outputPaternalReference1Path = Paths.get(wd.getCanonicalPath(), "1_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        if (updateVCF) {
            Files.copy(outputVCF1Path, Paths.get(expectedVCF1), StandardCopyOption.REPLACE_EXISTING);
        }
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReference1Path.toFile(), new File(maternalReference1)));
        assertTrue(FileUtils.contentEquals(outputPaternalReference1Path.toFile(), new File(paternalReference1)));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }

    /**
     * test male sample with variants on X chromosome
     * @param directory
     * @throws IOException
     */
    public void universalTestMethod3(String directory) throws IOException {
        File wd = tmpFolder.newFolder("tmp");
        String reference = new File(directory, "reference.fa").toString();
        String vcf = new File(directory, "input.vcf").toString();
        String expectedVCF1 = new File(directory, "expected.vcf").toString();
        String map = new File(directory, "expected.map").toString();
        String maternalReference1 = new File(directory, "chrX_maternal.fa").toString();

        Path outputVCF1Path = Paths.get(wd.getCanonicalPath(), "chrX_test.vcf");
        Path outputMaternalReferencePath = Paths.get(wd.getCanonicalPath(), "chrX_test_maternal.fa");
        Path outputPaternalReferencePath = Paths.get(wd.getCanonicalPath(), "chrX_test_paternal.fa");

        VCF2diploid runner = new VCF2diploid();
        String[] args = new String[]{
                "-chr", reference, "-outdir", wd.getCanonicalPath(),
                "-seed", Integer.toString(this.seed), "-id", "test",
                "-t", "MALE", "-vcf", vcf
        };
        runner.run(args);
        if (updateVCF) {
            Files.copy(outputVCF1Path, Paths.get(expectedVCF1), StandardCopyOption.REPLACE_EXISTING);
        }
        assertTrue(FileUtils.contentEquals(outputVCF1Path.toFile(), new File(expectedVCF1)));
        assertTrue(FileUtils.contentEquals(outputMaternalReferencePath.toFile(), new File(maternalReference1)));
        assertTrue(Files.notExists(outputPaternalReferencePath));
        assertTrue(FileUtils.contentEquals(runner.getOutputMap(), new File(map)));
    }
    @Rule
    public TemporaryFolder tmpFolder = new TemporaryFolder();

    /*
    assume maven will automatically go to root directory of the project.
    based on http://stackoverflow.com/questions/28673651/how-to-get-the-path-of-src-test-resources-directory-in-junit
    it might not be a good idea to use ClassLoader (they are not java.io.File objects in a jar file).
     */
    @Test
    public void SNPTest() throws IOException {
        universalTestMethod2("src/test/resources/SNPTest");
    }
    @Test
    public void MNPTest() throws IOException {
        universalTestMethod2("src/test/resources/MNPTest");
    }
    @Test
    public void InsertionTest() throws IOException {
        universalTestMethod2("src/test/resources/InsertionTest");
    }
    @Test
    public void DeletionTest() throws IOException {
        universalTestMethod2("src/test/resources/DeletionTest");
    }
    @Test
    public void InversionTest() throws IOException {
        universalTestMethod2("src/test/resources/InversionTest");
    }
    @Test
    public void DuplicationTest() throws IOException {
        universalTestMethod2("src/test/resources/DuplicationTest");
    }
    @Test
    public void MixedTest() throws IOException {
        universalTestMethod2("src/test/resources/MixedTest");
    }
    @Test
    public void balancedReciprocalTranslocationIntrachromosomalTest() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomal");
    }
    @Test
    public void balancedReciprocalTranslocationIntrachromosomalWithInversionTest() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationIntrachromosomalWithInversion");
    }
    @Test
    public void balancedReciprocalTranslocationInterchromosomalTest() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomal");
    }
    @Test
    public void balancedReciprocalTranslocationInterchromosomalWithInversionTest() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedReciprocalTranslocationTest/balancedReciprocalTranslocationInterchromosomalWithInversion");
    }
    @Test
    public void balancedNonReciprocalTranslocationInterchromosomal() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomal");
    }
    @Test
    public void balancedNonReciprocalTranslocationInterchromosomalWithInversionTest() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterchromosomalWithInversion");
    }
    @Test
    public void balancedNonReciprocalTranslocationInterHomologousChromosomal() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationInterHomologousChromosomal");
    }
    @Test
    public void balancedNonReciprocalTranslocationIntrachromosomal() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomal");
    }
    @Test
    public void balancedNonReciprocalTranslocationIntrachromosomalWithInversion() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/BalancedNonReciprocalTranslocationTest/balancedNonReciprocalTranslocationIntrachromosomalWithInversion");
    }
    @Test
    public void unbalancedTranslocationInterchromosomal() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedLossTranslocationTest/unbalancedTranslocationInterchromosomal");
    }
    @Test
    public void unbalancedTranslocationIntrachromosomal() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedLossTranslocationTest/unbalancedTranslocationIntrachromosomal");
    }

    /**
     * unbalanced loss
     * @throws IOException
     */
    @Test
    public void unbalancedTranslocationInterchromosomalWithInversion() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedLossTranslocationWithInversionTest/UnbalancedLossInterchromosomalTranslocationWithInversion");
    }
    @Test
    public void unbalancedTranslocationIntrachromosomalWithInversion() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedLossTranslocationWithInversionTest/UnbalancedLossIntrachromosomalTranslocationWithInversion");
    }
    @Test
    public void unbalancedGainTranslocationInterchromosomal() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedGainTranslocationTest/UnbalancedGainInterchromosomalTranslocation");
    }
    @Test
    public void unbalancedGainTranslocationIntrachromosomal() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedGainTranslocationTest/UnbalancedGainIntrachromosomalTranslocation");
    }
    @Test
    public void unbalancedGainTranslocationInterchromosomalWithInversion() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedGainTranslocationTest/UnbalancedGainInterchromosomalTranslocationWithInversion");
    }
    @Test
    public void unbalancedGainTranslocationIntrachromosomalWithInversion() throws IOException {
        universalTestMethod("src/test/resources/TranslocationTest/UnbalancedGainTranslocationTest/UnbalancedGainIntrachromosomalTranslocationWithInversion");
    }
    @Test
    public void homDelOverlapOthers() throws IOException {
        universalTestMethod2("src/test/resources/simulationTests/homDelOverlapOthers");
    }
    @Test
    public void hetDelOverlapOthers() throws IOException {
        universalTestMethod2("src/test/resources/simulationTests/hetDelOverlapOthers");
    }
    @Test
    public void maleChrXHet() throws IOException {
        universalTestMethod3("src/test/resources/simulationTests/maleChrXHet");
    }
    @Test
    public void cnvTest() throws IOException {
        universalTestMethod("src/test/resources/simulationTests/cnvTest");
    }
}
