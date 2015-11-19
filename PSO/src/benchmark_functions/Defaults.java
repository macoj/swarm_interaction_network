/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */

package benchmark_functions;
/**
 * This internal class holds the default values for everything
 * 
 * @author Thomas Weise
 */
final class Defaults {

  /** the default dimension */
  static final int DEFAULT_DIM = 1000;

  /** the default m */
  static int DEFAULT_M = 50;

  /**
   * Obtain the randomizer
   * 
   * @param c
   *          the class to get the randomizer for
   * @return the randomizer
   */
  static final Randomizer getRandomizer(final Class/* <?> */c) {
    long l;

    try {
      l = Long.parseLong(c.getSimpleName().substring(1));
    } catch (Throwable tt) {
      throw new RuntimeException(tt);
    }

    return new Randomizer(l);
  }

  // /** the full directory path of this package */
  // private static final String FULL_DIR = Kernel.class.getPackage()
  // .getName().replace('.', File.pathSeparatorChar);
  //
  // /** the rotation matrix suffix */
  //  private static final String ROT_MAT_SUFFIX = "_M.txt"; //$NON-NLS-1$
  //
  // /** the shift vector suffix */
  //  private static final String SHIFT_VECTOR_SUFFIX = "_o.txt"; //$NON-NLS-1$
  //
  // /** the permutation vector suffix */
  //  private static final String PERM_VECTOR_SUFFIX = "_P.txt"; //$NON-NLS-1$
  //
  // /**
  // * Obtain the default rotation matrix for class c
  // *
  // * @param c
  // * the class
  // * @return the default rotation matrix for class c
  // */
  // static final double[] defaultRotMatrix(final Class<?> c) {
  // final InputStream is;
  // final String s;
  // final double[] d;
  //
  // s = c.getSimpleName() + ROT_MAT_SUFFIX;
  //
  // try {
  // is = Kernel.class.getResourceAsStream(s);
  // if (is != null) {
  // d = Utils.loadRotMatrix(is);
  // } else {
  // d = Utils.loadRotMatrix(FULL_DIR + File.pathSeparatorChar + s);
  // }
  // DEFAULT_M = (int) (0.9d + Math.sqrt(d.length));
  // return d;
  // } catch (Throwable t) {
  // return Utils.defaultRandomRotMatrix(DEFAULT_M);
  // }
  // }
  //
  // /**
  // * Obtain the default shift vector for class c
  // *
  // * @param c
  // * the class
  // * @return the default shift vector for class c
  // */
  // static final double[] defaultShiftVec(final Class<?> c) {
  // final InputStream is;
  // final String s;
  // double min, max;
  // final double[] d;
  //
  // s = c.getSimpleName() + SHIFT_VECTOR_SUFFIX;
  //
  // try {
  // is = Kernel.class.getResourceAsStream(s);
  // if (is != null) {
  // d = Utils.loadShiftVector(is);
  // } else {
  // d = Utils.loadShiftVector(FULL_DIR + File.pathSeparatorChar + s);
  // }
  // return d;
  // } catch (Throwable t) {
  //
  // try {
  //        min = ((Number) (c.getDeclaredField("MIN").get(null))).doubleValue(); //$NON-NLS-1$      
  // } catch (Throwable tt) {
  // min = -100d;
  // }
  //
  // try {
  //        max = ((Number) (c.getDeclaredField("MAX").get(null))).doubleValue(); //$NON-NLS-1$      
  // } catch (Throwable tt) {
  // max = -100d;
  // }
  // return Utils.defaultRandomShiftVector(DEFAULT_DIM, min, max);
  // }
  // }
  //
  // /**
  // * Obtain the default permutation vector for class c
  // *
  // * @param c
  // * the class
  // * @return the default permutation vector for class c
  // */
  // static final int[] defaultPermVec(final Class<?> c) {
  // final InputStream is;
  // final String s;
  // final int[] d;
  //
  // s = c.getSimpleName() + PERM_VECTOR_SUFFIX;
  //
  // try {
  // is = Kernel.class.getResourceAsStream(s);
  // if (is != null) {
  // d = Utils.loadPermVector(is);
  // } else {
  // d = Utils.loadPermVector(FULL_DIR + File.pathSeparatorChar + s);
  // }
  // return d;
  // } catch (Throwable t) {
  // return Utils.defaultRandomPermVector(DEFAULT_DIM);
  // }
  // }

}
