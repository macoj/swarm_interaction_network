/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * The shifted rastrigin's function: F2.
 * 
 * @author Thomas Weise
 */
public final class F2 extends ShiftedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 5d;

  /** the minimum value */
  public static final double MIN = (-MAX);

  /**
   * Create a new shifted rastrigin's function
   * 
   * @param o
   *          the shifted global optimum
   */
  public F2(final double[] o) {
    super(o, MIN, MAX);
  }
  
  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F2(int dimension) {
	    this(Defaults.getRandomizer(F2.class).createShiftVector(dimension, MIN, MAX));
  }

  /**
   * Create a default instance of F2.
   * 
   * @param r
   *          the randomizer to use
   */
  public F2(final Randomizer r) {
    this(r.createShiftVector(Defaults.DEFAULT_DIM, MIN, MAX));//
  }

  /**
   * Create a default instance of F2.
   */
  public F2() {
    this(Defaults.getRandomizer(F2.class));
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
    return Kernel.shiftedRastrigin(x, this.m_o, 0, this.m_dimension);
  }

  /**
   * Obtain the full name of the benchmark function (according to
   * &quot;Benchmark Functions for the CEC���2010 Special Session and
   * Competition on Large-Scale Global Optimization&quot; Ke Tang, Xiaodong
   * Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise CEC'2010)
   * 
   * @return the full name of the benchmark function
   */
  // @Override
  public final String getFullName() {
    return "Shifted Rastrigin���s Function";//$NON-NLS-1$
  }

  /**
   * Obtain the short name of the benchmark function (according to
   * &quot;Benchmark Functions for the CEC���2010 Special Session and
   * Competition on Large-Scale Global Optimization&quot; Ke Tang, Xiaodong
   * Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise CEC'2010)
   * 
   * @return the short name of the benchmark function
   */
  // @Override
  public final String getShortName() {
    return "F2"; //$NON-NLS-1$
  }
}
