package com.bina.varsim.tools.evaluation;

import com.bina.varsim.tools.simulation.VCF2diploid;
import org.apache.commons.io.FileUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import static htsjdk.samtools.SAMFileHeader.GroupOrder.reference;
import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 10/28/16.
 */
public class VCFCompareTest {
  @Rule
  public TemporaryFolder tmpFolder = new TemporaryFolder();

  public void universalTestMethod(String directory) throws IOException {
    File wd = tmpFolder.newFolder("tmp");
    String truthVcf = new File(directory, "truth.vcf").toString();
    String vcfForCompare = new File(directory, "compare.vcf").toString();
    String expectedFalseNegative = new File(directory, "test_FN.vcf").toString();
    String expectedFalsePositive = new File(directory, "test_FP.vcf").toString();
    String expectedTruePositive = new File(directory, "test_TP.vcf").toString();
    String expectedUnknownFalsePositive = new File(directory, "test_unknown_FP.vcf").toString();
    String expectedUnknownTruePositive = new File(directory, "test_unknown_TP.vcf").toString();

    Path outputFalseNegative = Paths.get(wd.getCanonicalPath(), "test_FN.vcf");
    Path outputFalsePositive = Paths.get(wd.getCanonicalPath(), "test_FP.vcf");
    Path outputTruePositive = Paths.get(wd.getCanonicalPath(), "test_TP.vcf");
    Path outputUnknownFalsePositive = Paths.get(wd.getCanonicalPath(), "test_unknown_FP.vcf");
    Path outputUnknownTruePositive = Paths.get(wd.getCanonicalPath(), "test_unknown_TP.vcf");

    //example: java -jar ~/Downloads/varsim_0.6.3/VarSim.jar vcfcompare -true_vcf simu.truth.vcf -prefix try small.lumpy.vcf
    String[] args = new String[]{
            "-true_vcf", truthVcf,
            "-prefix", Paths.get(wd.getCanonicalPath(),"test").toString(),
            vcfForCompare
    };
    VCFcompare.main(args);
    assertTrue(FileUtils.contentEquals(outputFalseNegative.toFile(), new File(expectedFalseNegative)));
    assertTrue(FileUtils.contentEquals(outputFalsePositive.toFile(), new File(expectedFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputTruePositive.toFile(), new File(expectedTruePositive)));
    assertTrue(FileUtils.contentEquals(outputUnknownFalsePositive.toFile(), new File(expectedUnknownFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputUnknownTruePositive.toFile(), new File(expectedUnknownTruePositive)));
  }
  @Test
  public void firstSimpleTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/identicalTest");
  }
}
