/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * The Shifted Schwefel���s Problem 1.2: F19.
 * 
 * @author Thomas Weise
 */
public final class F19 extends ShiftedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 100d;

  /** the minimum value */
  public static final double MIN = (-MAX);

  /**
   * Create a new Shifted Schwefel���s Problem 1.2
   * 
   * @param o
   *          the shifted global optimum
   */
  public F19(final double[] o) {
    super(o, MIN, MAX);
  }
  
  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F19(int dimension) {
    this(Defaults.getRandomizer(F19.class).createShiftVector(dimension, MIN, MAX));
  }

  /**
   * Create a default instance of F19.
   * 
   * @param r
   *          the randomizer to use
   */
  public F19(final Randomizer r) {
    this(r.createShiftVector(Defaults.DEFAULT_DIM, MIN, MAX));//
  }

  /**
   * Create a default instance of F19.
   */
  public F19() {
    this(Defaults.getRandomizer(F19.class));
  }

  /**
   * Compute the value of the elliptic function. This function takes into
   * consideration only the first {{@link #getDimension()} elements of the
   * candidate vector.
   * 
   * @param x
   *          the candidate solution vector
   * @return the value of the function
   */
  // @Override
  public final double compute(final double[] x) {
    return Kernel.shiftedSchwefel12(x, this.m_o, 0, this.m_dimension);
  }

  /**
   * Obtain the full name of the benchmark function (according to
   * &quot;Benchmark Functions for the CEC���2010 Special Session and
   * Competition on Large-Scale Global Optimization&quot; , Ke Tang,
   * Xiaodong Li, P. N. Suganthan, and Zhenyu Yang, CEC'2010)
   * 
   * @return the full name of the benchmark function
   */
  // @Override
  public final String getFullName() {
    return "Shifted Schwefel���s Problem 1.2";//$NON-NLS-1$
  }

  /**
   * Obtain the short name of the benchmark function (according to
   * &quot;Benchmark Functions for the CEC���2010 Special Session and
   * Competition on Large-Scale Global Optimization&quot; , Ke Tang,
   * Xiaodong Li, P. N. Suganthan, and Zhenyu Yang, CEC'2010). If no short
   * name is defined, the full name will be used.
   * 
   * @return the short name of the benchmark function
   */
  // @Override
  public final String getShortName() {
    return "F19"; //$NON-NLS-1$
  }
}
