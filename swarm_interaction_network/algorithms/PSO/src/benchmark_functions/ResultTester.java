/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 * The tester for comparing the benchmark functions' results
 * 
 * @author Thomas Weise
 */
public class ResultTester {

  /** the test functions */
  public static final Function[] TESTS = new Function[] {// 
  // new F1(),
  // new F2(),
  // new F3(),
//   new F4(),
  // new F5(),
  // new F6(),
  // new F7(),
//   new F8(),
  // new F9(),
  // new F10(),
  // new F11(),
  // new F12(),
//   new F13(),
  // new F14(),
  // new F15(),
  // new F16(),
  // new F17(),
  // new F18(),
  // new F19(),
  new F20(),//
  };

  /**
   * The test function used to check whether the routines here have been
   * implemented correctly.
   * 
   * @param params
   *          the parameters
   */
  public static final void main(final String[] params) {
    final double[][] d;
    int i, j;

    d = loadData();

    for (i = 0; i < TESTS.length; i++) {
      System.out.println();
      System.out.println(TESTS[i].getShortName());
      for (j = 0; j < d.length; j++) {
        System.out.println();
        System.out.print('\t');
        System.out.print(TESTS[i].compute(d[j]));
      }
    }
  }

  /**
   * load the test data
   * 
   * @return the test data
   */
  private static final double[][] loadData() {
    InputStream is;
    BufferedReader r;
    String s, os;
    String[] ss;
    ArrayList/* <double[]> */ds;
    double[] d;
    int i;

    ds = null;

    try {
      is = ResultTester.class.getResourceAsStream("f20_opt.txt"); //$NON-NLS-1$
      if (is == null) {
        is = new FileInputStream(ResultTester.class.getPackage().getName()
            .replace('.', File.pathSeparatorChar)
            + File.pathSeparatorChar + "f20_opt.txt"); //$NON-NLS-1$
      }
      try {
        r = new BufferedReader(new InputStreamReader(is));
        try {

          ds = new ArrayList/* <double[]> */(100);

          while ((s = r.readLine()) != null) {
            s = s.trim();
            if (s.length() <= 0) {
              continue;
            }

            do {
              os = s;
              s = s.replace('\t', ' ').replace("  ", " "); //$NON-NLS-1$//$NON-NLS-2$
            } while (s != os);

            ss = s.split(" "); //$NON-NLS-1$
            if (ss.length != Defaults.DEFAULT_DIM) {
              continue;
            }

            d = new double[Defaults.DEFAULT_DIM];
            for (i = (d.length - 1); i >= 0; i--) {
              d[i] = Double.parseDouble(ss[i].trim());
            }

            ds.add(d);
          }

        } finally {
          r.close();
          is = null;
        }
      } finally {
        if (is != null) {
          is.close();
        }
      }

      if (ds != null) {
        return ((double[][]) (ds.toArray(new double[ds.size()][])));
      }
      return null;
    } catch (Throwable tt) {
      throw new RuntimeException(tt);
    }

  }

}
