/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
import java.io.File;

/**
 * A function which is shifted and permutated
 * 
 * @author Thomas Weise
 */
public abstract class ShiftedPermutatedFunction extends ShiftedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the permutation */
  protected final int[] m_p;

  /**
   * Create a new shifted function
   * 
   * @param o
   *          the shifted global optimum
   * @param p
   *          the permutation vector
   * @param min
   *          the minimum value a decision variable can take on
   * @param max
   *          the maximum value a decision variable can take on
   */
  protected ShiftedPermutatedFunction(final double[] o, final int[] p,
      final double min, final double max) {
    super(o, min, max);
    this.m_p = p;
  }

  /**
   * Obtain the permutation vector of the benchmark function
   * 
   * @return the permutation of the benchmark function
   */
  public final int[] getPermutationVector() {
    return this.m_p;
  }

  /**
   * Store the utility information
   * 
   * @param f
   *          the directory
   * @throws Throwable
   *           a possible io exception
   */
  // @Override
  void storeUtilityInfo(final File f) throws Throwable {
    Throwable s;

    s = null;

    try {
      super.storeUtilityInfo(f);
    } catch (Throwable t) {
      s = t;
    }

    try {
      Utils.storePermVector(new File(f, this.getShortName() + // 
          "_p.txt").toString(), this.m_p); //$NON-NLS-1$
    } catch (Throwable t) {
      if (s == null) {
        s = t;
      }
    }

    if (s != null) {
      throw s;
    }
  }

}
