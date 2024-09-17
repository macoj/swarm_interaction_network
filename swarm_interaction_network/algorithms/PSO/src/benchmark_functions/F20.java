/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * The Shifted Rosenbrock���s Function: F20.
 * 
 * @author Thomas Weise
 */
public final class F20 extends ShiftedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 100d;

  /** the minimum value */
  public static final double MIN = (-MAX);

  /** the optimum */
  private transient double[] m_opt;
  
  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F20(int dimension) {
	    this(Defaults.getRandomizer(F20.class).createShiftVector(dimension, MIN, MAX - 1d));
  }

  /**
   * Create a new Shifted Rosenbrock���s Function
   * 
   * @param o
   *          the shifted global optimum
   */
  public F20(final double[] o) {
    super(o, MIN, MAX);
  }

  /**
   * Create a default instance of F20.
   * 
   * @param r
   *          the randomizer to use
   */
  public F20(final Randomizer r) {
    this(r.createShiftVector(Defaults.DEFAULT_DIM, MIN, MAX - 1d));//
  }

  /**
   * Create a default instance of F20.
   */
  public F20() {
    this(Defaults.getRandomizer(F20.class));
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
    return Kernel.shiftedRosenbrock(x, this.m_o, 0, this.m_dimension);
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
    return "Shifted Rosenbrock���s Function";//$NON-NLS-1$
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
    return "F20"; //$NON-NLS-1$
  }

  /**
   * Obtain the optimum vector of the benchmark function
   * 
   * @return the optimum vector of the benchmark function
   */
  // @Override
  public final double[] getOptimum() {
    double[] d;
    int i;

    d = this.m_opt;
    if (d != null) {
      return d;
    }

    d = ((double[]) (this.m_o.clone()));
    for (i = (d.length - 1); i >= 0; i--) {
      d[i] += 1d;
    }

    return (this.m_opt = d);
  }
}
