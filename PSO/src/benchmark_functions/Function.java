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
import java.io.FileWriter;
import java.io.Serializable;

/**
 * The base class for benchmark functions. All benchmark functions are for
 * minimization and have the same global optimum: <code>0</code>. Some
 * benchmark functions may use internal temporary variables and
 * <warning>not be threadsafe</warning>. Therefore, you should always use
 * each instance of a function for one single {#link java.lang.Thread}
 * only. You may clone or serialize function instances.
 * 
 * @author Thomas Weise 
 */
public abstract class Function implements Serializable, Cloneable {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the worst possible objective value */
  public static final double WORST = Double.POSITIVE_INFINITY;

  /** the best possible objective value */
  public static final double BEST = 0d;

  /** the dimension of the function */
  protected final int m_dimension;

  /** the minimum value which the decision variables can take on */
  protected final double m_min;

  /** the maximum value which the decision variables can take on */
  protected final double m_max;

  /**
   * Create the benchmark function
   * 
   * @param dimension
   *          the dimension of the function
   * @param min
   *          the minimum value which the decision variables can take on
   * @param max
   *          the maximum value which the decision variables can take on
   */
  protected Function(final int dimension, final double min,
      final double max) {
    super();
    this.m_dimension = dimension;
    this.m_min = min;
    this.m_max = max;
  }

  /**
   * Obtain the dimension of the function.
   * 
   * @return the dimension of the function
   */
  public final int getDimension() {
    return this.m_dimension;
  }

  /**
   * Obtain the minimum value which the decision variables can take on
   * 
   * @return the minimum value which the decision variables can take on
   */
  public final double getMin() {
    return this.m_min;
  }

  /**
   * Obtain the maximum value which the decision variables can take on
   * 
   * @return the maximum value which the decision variables can take on
   */
  public final double getMax() {
    return this.m_max;
  }

  /**
   * Compute the value of the benchmark function. This function takes into
   * consideration only the first {{@link #getDimension()} elements of the
   * candidate vector.
   * 
   * @param x
   *          the candidate solution vector
   * @return the value of the function
   */
  public abstract double compute(final double[] x);

  /**
   * Obtain the full name of the benchmark function (according to
   * &quot;Benchmark Functions for the CEC�2010 Special Session and
   * Competition on Large-Scale Global Optimization&quot; Ke Tang, Xiaodong
   * Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise CEC'2010)
   * 
   * @return the full name of the benchmark function
   */
  public abstract String getFullName();

  /**
   * Obtain the short name of the benchmark function (according to
   * &quot;Benchmark Functions for the CEC�2010 Special Session and
   * Competition on Large-Scale Global Optimization&quot; Ke Tang, Xiaodong
   * Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise CEC'2010)
   * 
   * @return the short name of the benchmark function
   */
  public String getShortName() {
    return this.getFullName();
  }

  /**
   * Obtain the optimum vector of the benchmark function
   * 
   * @return the optimum vector of the benchmark function
   */
  public abstract double[] getOptimum();

  /**
   * the internal clone routine
   * 
   * @return a clone of this object
   */
  final Function internalClone() {
    try {
      return ((Function) (super.clone()));
    } catch (Throwable t) {
      throw new RuntimeException(t); // will not happen since we implement
      // Cloneable
    }
  }

  /**
   * Clone this function
   * 
   * @return a clone of this function
   */
  // //@Override
  public Object clone() {
    return this; // default behavior: return this
  }

  /**
   * Store the information of this function in the following directory
   * 
   * @param dir
   *          the directory
   */
  public final void storeInfo(final String dir) {
    File f;
    BufferedWriter w;
    Throwable s;

    s = null;
    try {
      f = new File(dir).getCanonicalFile();

      try {
        w = new BufferedWriter(new FileWriter(//
            new File(f, this.getShortName() + ".txt"))); //$NON-NLS-1$
        try {
          this.storeFunctionInfo(w);
        } finally {
          w.close();
        }
      } catch (Throwable tt) {
        s = tt;
      }

      this.storeUtilityInfo(f);

    } catch (Throwable t) {
      if (s == null) {
        s = t;
      }
    }
    if (s != null) {
      if ((s instanceof RuntimeException)) {
        throw ((RuntimeException) s);
      }
      throw new RuntimeException(s);
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
  protected void storeFunctionInfo(final BufferedWriter w)
      throws Throwable {

    w.write("Short Name: ");//$NON-NLS-1$
    w.write(this.getShortName());
    w.newLine();
    w.write("Long Name: ");//$NON-NLS-1$
    w.write(this.getFullName());
    w.newLine();
    w.write("Dimension: ");//$NON-NLS-1$
    w.write(Integer.toString(this.m_dimension));
    w.newLine();
    w.write("x-min    : ");//$NON-NLS-1$
    w.write(Double.toString(this.m_min));
    w.newLine();
    w.write("x-max    : "); //$NON-NLS-1$
    w.write(Double.toString(this.m_min));
  }

  /**
   * Store the utility information
   * 
   * @param f
   *          the directory
   * @throws Throwable
   *           a possible io exception
   */
  void storeUtilityInfo(final File f) throws Throwable {
    //
  }

  /**
   * Obtain the string representation of this function's name
   * 
   * @return the string representation of this function's name
   */
  // //@Override
  public final String toString() {
    return this.getShortName() + ": " + this.getFullName(); //$NON-NLS-1$
  }
}
