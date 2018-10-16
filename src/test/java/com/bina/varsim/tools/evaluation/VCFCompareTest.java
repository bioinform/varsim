package com.bina.varsim.tools.evaluation;

import com.bina.varsim.tools.simulation.VCF2diploid;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Ignore;
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
import static htsjdk.samtools.SAMFileHeader.GroupOrder.reference;
import static junit.framework.TestCase.assertTrue;

/**
 * Created by guoy28 on 10/28/16.
 */
public class VCFCompareTest {
  @Rule
  public TemporaryFolder tmpFolder = new TemporaryFolder();

  public void universalTestMethod(String directory) throws IOException {
    universalTestMethod(directory, new String[0]);
  }
  public void universalTestMethod(String directory, String[] additionalArgs) throws IOException {
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
    Path outputJson = Paths.get(wd.getCanonicalPath(), "test_report.json");

    //example: java -jar ~/Downloads/varsim_0.6.3/VarSim.jar vcfcompare -true_vcf simu.truth.vcf -prefix try small.lumpy.vcf
    String[] args = new String[]{
            "-true_vcf", truthVcf,
            "-prefix", Paths.get(wd.getCanonicalPath(),"test").toString(),
            vcfForCompare
    };
    VCFcompare.main(ArrayUtils.addAll(args, additionalArgs));
    if (updateVCF) {
      Files.copy(outputFalseNegative, Paths.get(expectedFalseNegative), StandardCopyOption.REPLACE_EXISTING);
      Files.copy(outputFalsePositive, Paths.get(expectedFalsePositive), StandardCopyOption.REPLACE_EXISTING);
      Files.copy(outputTruePositive, Paths.get(expectedTruePositive), StandardCopyOption.REPLACE_EXISTING);
      Files.copy(outputUnknownTruePositive, Paths.get(expectedUnknownTruePositive), StandardCopyOption.REPLACE_EXISTING);
      Files.copy(outputUnknownFalsePositive, Paths.get(expectedUnknownFalsePositive), StandardCopyOption.REPLACE_EXISTING);
    }
    assertTrue(FileUtils.contentEquals(outputFalseNegative.toFile(), new File(expectedFalseNegative)));
    assertTrue(FileUtils.contentEquals(outputFalsePositive.toFile(), new File(expectedFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputTruePositive.toFile(), new File(expectedTruePositive)));
    assertTrue(FileUtils.contentEquals(outputUnknownFalsePositive.toFile(), new File(expectedUnknownFalsePositive)));
    assertTrue(FileUtils.contentEquals(outputUnknownTruePositive.toFile(), new File(expectedUnknownTruePositive)));
  }

  /**
   * comapre two identical VCFs
   * @throws IOException
   */
  @Test
  public void firstSimpleTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/identicalTest");
  }

  /**
   * test 2 variants with same mutations but different genotypes
   * -match_geno off
   * @throws IOException
   */
  @Test
  public void genotypeTest1() throws IOException {
    universalTestMethod("src/test/resources/validationTest/genotypeTest1");
  }
  /**
   * test 2 variants with same mutations but different genotypes
   * -match_geno on
   * @throws IOException
   */
  @Test
  public void genotypeTest2() throws IOException {
    universalTestMethod("src/test/resources/validationTest/genotypeTest2", new String[]{"-match_geno"});
  }

  /**
   * comparing two VCFs with two same variants but different reference padding at the beginning
   * @throws IOException
   */
  @Test
  public void differentPaddingTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/differentPaddingTest");
  }
  /**
   * compare two VCFs with two same left-normalized variants (more work than removing padding is required)
   * a variant in truth VCF will be split and compared against 2 variants in compare.vcf
   * @throws IOException
   */
  @Test
  public void canonicalizationTest1() throws IOException {
    universalTestMethod("src/test/resources/validationTest/canonicalizationTest1");
  }
  /**
   * compare two VCFs with two same left-normalized variants (more work than removing padding is required)
   * a variant in truth VCF will be split and compared against 1 variant in compare.vcf
   * @throws IOException
   */
  @Test
  public void canonicalizationTest2() throws IOException {
    universalTestMethod("src/test/resources/validationTest/canonicalizationTest2");
  }
  /**
   * compare two VCFs with two same left-normalized variants (more work than removing padding is required)
   * 2 variants in truth VCF will be compared against 1 variant in compare.vcf
   * @throws IOException
   */
  @Test
  public void canonicalizationTest3() throws IOException {
    universalTestMethod("src/test/resources/validationTest/canonicalizationTest3");
  }
  /**
   * truth is an 2bp MNP
   * compared against 2 SNPs are actually identical to the MNP (same phase, same genotype)
   * @throws IOException
   */
  @Test
  public void canonicalizationTest4SNPvsMNP() throws IOException {
    universalTestMethod("src/test/resources/validationTest/canonicalizationTest4SNPvsMNP");
  }
  /**
   * truth is an 2bp MNP
   * compared against 2 SNPs are actually identical to the MNP (same phase, same genotype)
   * @throws IOException
   */
  @Test
  public void canonicalizationTest5SNPvsMNPnoGT() throws IOException {
    universalTestMethod("src/test/resources/validationTest/canonicalizationTest5SNPvsMNPnoGT");
  }

  /**
   * test 2 SNPs when both maternal and paternal genotypes are missing
   * @throws IOException
   */
  @Test
  public void missingGTSNPTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/missingGT_SNP");
  }
  /**
   * test 2 SNPs with identical loci
   * @throws IOException
   */
  @Test
  public void duplicateLociSNPTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/sameLociSNP");
  }
  /**
   * compare two VCFs with two same balanced, nonreciprocal intrachromosomal translocations
   * (cut-and-paste). Validation will be done by evaluating breakend.
   * @throws IOException
   */
  @Test
  public void breakendEvaluationOnIdenticalBalancedNonreciprocalTranslocation() throws IOException {
    universalTestMethod("src/test/resources/validationTest/breakendTests/breakendEvaluationOnIdenticalBalancedNonreciprocalTranslocation");
  }
  /**
   * compare two VCFs with one insertion and one deletion
   * @throws IOException
   */
  @Test
  public void insertionVsInsertionTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/insertionVsInsertionTest", new String[]{"-wig", "20", "-over", "0.2"});
  }
  /**
   * compare two VCFs with 2 translocations, same source of duplication, different insertion positions
   * let it fail with a small wiggle setting
   * @throws IOException
   */
  @Test
  public void breakendEvaluationDifferentInsertionsFailWithSmallWiggle() throws IOException {
    universalTestMethod("src/test/resources/validationTest/breakendTests/nonreciprocalWithDifferentInsertionFailWithSmallWiggle", new String[]{"-wig", "1"});
  }
  /**
   * compare two VCFs with 2 translocations, same source of duplication, different insertion positions
   * let it succeed with a large wiggle
   * @throws IOException
   */
  @Test
  public void breakendEvaluationDifferentInsertionsSucceedWithLargeWiggle() throws IOException {
    universalTestMethod("src/test/resources/validationTest/breakendTests/nonreciprocalWithDifferentInsertionSucceedWithLargeWiggle", new String[]{"-wig", "15"});
  }
  /**
   * compare two VCFs with 2 translocations, different sources of duplication, same insertion positions
   * duplications come from two chromosomes, so must fail
   * @throws IOException
   */
  @Test
  public void breakendEvaluationDifferentDuplicationChromosomes() throws IOException {
    universalTestMethod("src/test/resources/validationTest/breakendTests/nonreciprocalWithDifferentDuplicatedSequences");
  }
  /**
   * compare two VCFs with 2 translocations, same sources of duplication, same insertion positions
   * duplications are in different orientations
   * the two translocations (aka interspersed duplications) should not match even if wiggle is very
   * large, because of different orientations in breakend representations.
   * @throws IOException
   */
  @Test
  public void breakendEvaluationDifferentDuplicationOrientations() throws IOException {
    universalTestMethod("src/test/resources/validationTest/breakendTests/nonreciprocalWithDifferentOrientationsAndLargeWiggle", new String[]{"-wig", "100"});
  }
  /**
   * validate one reciprocal balanced translocation with another
   *
   * there are 4 records
   *
   */
  @Test
  public void breakendEvaluationNonReciprocalBalancedWithLargeWiggle() throws IOException {
    universalTestMethod("src/test/resources/validationTest/breakendTests/nonreciprocalBalancedWithLargeWiggle", new String[]{"-wig", "10"});
  }
  /******tests for distance-based metrics**************/
  /**
   * 1 SNP in TP, 1 SNP in compare, with 1-bp distance
   *
   */
  @Test
  public void distanceSNPTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/snpDistance", new String[]{"-output_distance_metric"});
  }
  /**
   * 1 INS in TP, 1 INS in compare, with 100-bp 3' distance, 100-bp 5' distance, 100bp length difference
   *
   */
  @Test
  public void distanceINSTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/insDistance", new String[]{"-wig", "150","-over","0.5","-output_distance_metric"});
  }
  /**
   * 1 INS in TP, 1 INS in compare, with 100-bp 3' distance, 100-bp 5' distance, 100bp length difference
   *
   * predicted length is imprecise
   *
   */
  @Test
  public void impreciseDistanceINSTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/impreciseInsDistance", new String[]{"-wig", "150","-over","0.5","-output_distance_metric","-ignore_imprecise_length"});
  }
  /**
   * 1 INS in TP, 1 INS in compare, with 100-bp 3' distance, 100-bp 5' distance, 100bp length difference
   *
   * truth length is imprecise
   *
   */
  @Test
  public void impreciseDistanceINSTest2() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/impreciseInsDistance2", new String[]{"-wig", "150","-over","0.5","-output_distance_metric","-ignore_imprecise_length"});
  }
  /**
   * 1 INV in TP, 1 INV in compare, with 1000-bp 3' distance, 11000-bp 5' distance, 100bp length difference
   *
   */
  @Test
  public void distanceINVTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/invDistance", new String[]{"-wig", "2000","-over","0.55","-output_distance_metric"});
  }
  /**
   * 1 DEL in TP with 2 ALTs, 1 DEL in compare with 1 ALT, matching only one of the 2 ALTS in TP
   *
   */
  @Test
  public void distanceDELTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/delDistance", new String[]{"-wig", "10","-over","0.5","-output_distance_metric"});
  }
  /**
   * 1 DUP in TP, 1 DUP in compare, length the same, distance 10bp
   *
   */
  @Test
  public void distanceDUPTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/dupDistance", new String[]{"-wig", "20","-over","0.7","-output_distance_metric"});
  }
  /**
   * 10 DUP in TP, 10 DUP in compare, length the same, distance 1,2,3,4,5,6,7,8,9,10bp
   *
   */
  @Test
  public void distanceTenDUPTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/tenDupDistance", new String[]{"-wig", "20","-over","0.7","-output_distance_metric"});
  }

  /**
   * same loci, different inserted sequences
   * @throws IOException
   */
  @Test
  public void sameLociDifferentInsertSeq() throws IOException {
    universalTestMethod("src/test/resources/validationTest/sameLociDifferentInsertSeq");
  }
  /**
   * same loci, same inserted sequences
   * @throws IOException
   */
  @Test
  public void sameLociSameInsertSeq() throws IOException {
    universalTestMethod("src/test/resources/validationTest/sameLociSameInsertSeq");
  }
  /**
   * 1 DEL in TP, 2 DEL in compare, length the same, distance 30bp,10bp, respectively
   * wiggle=50bp, both DEL in compare can match, but we only match the closest one
   */
  @Test
  @Ignore public void distanceMultipleMatchingTest() throws IOException {
    //disabled until global matching is implemented
    universalTestMethod("src/test/resources/validationTest/distanceMetricTests/multipleMatchingDistance", new String[]{"-wig", "50","-over","0.7","-output_distance_metric"});
  }
  /**
   * complex variant
   * CAAAA -> CAAG 1/1
   * is called as 1 DEL + 1 SNV
   * partial matching disallowed (so that rtg can rescue this case)
   */
  @Test
  public void complexVariantNoPartialMatchingTest() throws IOException {
    universalTestMethod("src/test/resources/validationTest/complexVariantTests/noPartialMatchingTest", new String[]{"-disallow_partial_fp"});
  }

  /**
   * truth is
   * 1 DEL + 1 SNV
   * complex variant is
   * CAAAA -> CAAG 1/1
   * partial matching disallowed (so that rtg can rescue this case)
   *
   * basically reverse of previous test case
   */
  @Test
  public void complexVariantNoPartialMatchingTest2() throws IOException {
    universalTestMethod("src/test/resources/validationTest/complexVariantTests/noPartialMatchingTest2", new String[]{"-disallow_partial_fp"});
  }
}
