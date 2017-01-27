package com.bina.varsim.tools.evaluation;

import com.bina.varsim.tools.evaluation.*;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import static com.bina.varsim.VarSimToolNamespace.SAMCompare;
import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 1/12/17.
 */
public class SAMCompareTest {
  @Rule
  public TemporaryFolder tmpFolder = new TemporaryFolder();

  public void universalTestMethodShortRead(String directory) throws IOException {
    universalTestMethod(directory, new String[0], false);

  }
  public void universalTestMethodShortRead(String directory, String[] additionalArgs) throws IOException {
    universalTestMethod(directory, additionalArgs, false);
  }
  public void universalTestMethodLongRead(String directory) throws IOException {
    universalTestMethod(directory, new String[0], true);
  }
  public void universalTestMethodLongRead(String directory, String[] additionalArgs) throws IOException {
    universalTestMethod(directory, additionalArgs, true);
  }
  public void universalTestMethod(String directory, String[] additionalArgs, boolean isLongRead) throws IOException {
    File wd = tmpFolder.newFolder("tmp");
    String bam = new File(directory, "input.bam").toString();
    String readmap = new File(directory, "input.map.gz").toString();

    String expectedDeletionFalsePositive = new File(directory, "test_D_FP.bam").toString();
    String expectedInsertionFalsePositive = new File(directory, "test_I_FP.bam").toString();
    String expectedSubstitutionFalsePositive = new File(directory, "test_S_FP.bam").toString();
    String expectedTumFalsePositive = new File(directory, "test_TUM_FP.bam").toString();
    String expectedTandemDupFalsePositive = new File(directory, "test_T_FP.bam").toString();
    String expectedInversionFalsePositive = new File(directory, "test_V_FP.bam").toString();
    String expectedFalsePositive = new File(directory, "test_FP.bam").toString();

    Path outputDeletionFalsePositive = Paths.get(wd.getCanonicalPath(), "test_D_FP.bam");
    Path outputInsertionFalsePositive = Paths.get(wd.getCanonicalPath(), "test_I_FP.bam");
    Path outputSubstitutionFalsePositive = Paths.get(wd.getCanonicalPath(), "test_S_FP.bam");
    Path outputTumFalsePositive = Paths.get(wd.getCanonicalPath(), "test_TUM_FP.bam");
    Path outputTandemDupFalsePositive = Paths.get(wd.getCanonicalPath(), "test_T_FP.bam");
    Path outputInversionFalsePositive = Paths.get(wd.getCanonicalPath(), "test_V_FP.bam");
    Path outputFalsePositive = Paths.get(wd.getCanonicalPath(), "test_FP.bam");

    //example: java -jar ~/Downloads/varsim_0.6.3/VarSim.jar vcfcompare -true_vcf simu.truth.vcf -prefix try small.lumpy.vcf
    String[] args = null;
    if (isLongRead) {
      args = new String[]{
              "-prefix", Paths.get(wd.getCanonicalPath(),"test").toString(),
              "-read_map",
              readmap,
              bam
      };
    } else {
      args = new String[]{
              "-prefix", Paths.get(wd.getCanonicalPath(),"test").toString(),
              bam,
      };
    }
    SAMcompare.main(ArrayUtils.addAll(args, additionalArgs));

    assertTrue(FileUtils.contentEquals(outputDeletionFalsePositive.toFile(), new File(expectedDeletionFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputInsertionFalsePositive.toFile(), new File(expectedInsertionFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputSubstitutionFalsePositive.toFile(), new File(expectedSubstitutionFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputTumFalsePositive.toFile(), new File(expectedTumFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputTandemDupFalsePositive.toFile(), new File(expectedTandemDupFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputInversionFalsePositive.toFile(), new File(expectedInversionFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputFalsePositive.toFile(), new File(expectedFalsePositive)));
  }

  /***************Short read test***************/
  /**
   * make sure the original behavior is not altered
   * @throws IOException
   */
  @Test
  public void firstSAMCompareShortReadTest() throws IOException{
   universalTestMethodShortRead("src/test/resources/samcompareTests/first_samcompare_short_read_test");
  }

  /**
   * there should be one false positive for insertion
   * @throws IOException
   */
  @Test
  public void oneFalsePositiveSAMCompareShortReadTest() throws IOException{
    universalTestMethodShortRead("src/test/resources/samcompareTests/short_read_test_one_false_positive",new String[]{"-wig","5"});
  }

  /**
   * there are 2 copies of a gene in the reference, although the 2 copies carry different mutations
   * some reads may appear equally like to come from either one of the copies.
   *
   * there should be one true positive and one false positive evaluated at 95% identity threshold
   * i.e. one true multialignment one false multialignment
   *
   * @throws IOException if file cannot be found/read
   */
  @Test
  public void oneFalsePositiveOneTruePositiveMultiAlignmentSAMCompareShortReadTest() throws IOException{
    universalTestMethodShortRead("src/test/resources/samcompareTests/multiAlignmentTests/shortReadoneTPoneFPMultiAlignment",new String[]{"-wig","5","-identity_threshold","0.95"});
  }


  /***************Long read test***************/
  /**
   * make sure the original behavior is not altered
   * @throws IOException
   */
  @Test
  public void firstSAMCompareLongReadTest() throws IOException {
    universalTestMethodLongRead("src/test/resources/samcompareTests/first_samcompare_long_read_test");
  }
}
