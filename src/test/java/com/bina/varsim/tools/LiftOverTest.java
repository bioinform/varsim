package com.bina.varsim.tools;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import static junit.framework.TestCase.assertTrue;

public class LiftOverTest {
  @Rule
  public TemporaryFolder tmpFolder = new TemporaryFolder();

  public void universalTestMethod(String directory) throws IOException {
    universalTestMethod(directory, new String[0]);
  }
  public void universalTestMethod(String directory, String[] additionalArgs) throws IOException {
    File wd = tmpFolder.newFolder("tmp");
    String inputMap = new File(directory, "input.map").toString();
    String inputBed = new File(directory, "input.bed").toString();
    String expectedBed = new File(directory, "expected.bed").toString();

    Path outputBed = Paths.get(wd.getCanonicalPath(), "test.bed");

    String[] args = new String[]{
            "-input", inputBed,
            "-map", inputMap,
            "-prefix", Paths.get(wd.getCanonicalPath(),"test").toString(),
    };
    LiftOver.main(ArrayUtils.addAll(args, additionalArgs));
    assertTrue(FileUtils.contentEquals(outputBed.toFile(), new File(expectedBed)));
  }

  /**
   * comapre two identical VCFs
   * @throws IOException
   */
  @Test
  public void firstSimpleTest() throws IOException {
    universalTestMethod("src/test/resources/LiftOverTest/simpleTest");
  }

  /**
   * lift over tandem duplications
   * @throws IOException
   */
  @Test
  public void tandemDupTest() throws IOException {
    universalTestMethod("src/test/resources/LiftOverTest/tandemDupTest");
  }

  /**
   * lift over large insertion
   * @throws IOException
   */
  @Test
  public void largeInsertionTest() throws IOException {
    universalTestMethod("src/test/resources/LiftOverTest/largeInsertionTest");
  }
}
