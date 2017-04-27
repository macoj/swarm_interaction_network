/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
import java.io.BufferedWriter;
import java.io.File;

/**
 * A shifted, permutated, and rotated function
 * 
 * @author Thomas Weise
 */
public abstract class ShiftedPermutatedRotatedFunction extends
    ShiftedPermutatedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the rotation matrix */
  protected final double[] m_m;

  /** the rotation matrix' dimension */
  protected final int m_matDim;

  /** a temporary vector, at least of the same length as the matrix is wide */
  protected double[] m_tmp;

  /**
   * Create a new shifted function
   * 
   * @param o
   *          the shifted global optimum
   * @param p
   *          the permutation vector
   * @param m
   *          the rotation matrix
   * @param min
   *          the minimum value a decision variable can take on
   * @param max
   *          the maximum value a decision variable can take on
   */
  protected ShiftedPermutatedRotatedFunction(final double[] o,
      final int[] p, final double[] m, final double min, final double max) {
    super(o, p, min, max);
    this.m_m = m;
    this.m_matDim = getMatrixDim(m);
    this.m_tmp = this.createTmp(this.m_matDim);
  }

  /**
   * Compute the matrix dimension
   * 
   * @param m
   *          the matrix
   * @return the matrix dimension
   */
  protected static final int getMatrixDim(final double[] m) {
    return (int) (0.9d + Math.sqrt(m.length));
  }

  /**
   * Create the temporary vector
   * 
   * @param w
   *          the width (and height) of the rotation matrix
   * @return the temporary vector. must be at least of length w but can be
   *         longer
   */
  protected double[] createTmp(final int w) {
    return new double[w];
  }

  /**
   * Clone this rotated function: Ensures that each rotated function has an
   * own temporary variable. This is needed if multiple threads are used
   * -&gt; in this case, each thread has to use its own instance of the
   * function!
   * 
   * @return a clone of this function
   */
  // @Override
  public Object clone() {
    ShiftedPermutatedRotatedFunction f;

    f = ((ShiftedPermutatedRotatedFunction) (this.internalClone()));
    f.m_tmp = f.createTmp((int) (0.9d + Math.sqrt(f.m_m.length)));
    return f;
  }

  /**
   * Obtain the rotation matrix of the benchmark function
   * 
   * @return the rotation matrix of the benchmark function
   */
  public final double[] getRotationMatrix() {
    return this.m_m;
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
      Utils.storeRotMatrix(new File(f, this.getShortName() + // 
          "_m.txt").toString(), this.m_m); //$NON-NLS-1$
    } catch (Throwable t) {
      if (s == null) {
        s = t;
      }
    }

    if (s != null) {
      throw s;
    }
  }

  /**
   * Store the function information
   * 
   * @param w
   *          the writer to store
   * @throws Throwable
   *           a possible io exception
   */
  // @Override
  protected void storeFunctionInfo(final BufferedWriter w)
      throws Throwable {
    String s;

    super.storeFunctionInfo(w);

    w.newLine();
    w.write("Rotation : ");//$NON-NLS-1$
    s = Integer.toString(this.m_matDim);
    w.write(s + 'x' + s);
  }
}
