/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * The D/2m-group Shifted m-dimensional Rosenbrock���s Function: F13.
 * 
 * @author Thomas Weise
 */
public final class F13 extends ShiftedPermutatedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 100d;

  /** the minimum value */
  public static final double MIN = (-MAX);

  /** the m-value */
  private final int m_m;
  
  /** the optimum vector */
  private transient double[] m_opt;
  
  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F13(int dimension) {
	    this(Defaults.getRandomizer(F13.class), dimension);
  }

  /**
   * Create a new D/2m--group v
   * 
   * @param o
   *          the shifted global optimum
   * @param p
   *          the permutation vector
   * @param m
   *          the fraction of nonseparability
   */
  public F13(final double[] o, final int[] p, final int m) {
    super(o, p, MIN, MAX);
    this.m_m = m;
  }

  /**
   * Create a default instance of F13.
   * 
   * @param r
   *          the randomizer to use
   */
  public F13(final Randomizer r) {
    this(r, Defaults.DEFAULT_DIM);//
  }
  
  public F13(final Randomizer r, int dimension) {
	    this(r.createShiftVector(dimension, MIN, MAX - 1d),//
	        r.createPermVector(dimension),//
	        Defaults.DEFAULT_M);//
  }

  /**
   * Create a default instance of F13.
   */
  public F13() {
    this(Defaults.getRandomizer(F13.class));
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
    final int max, gs, d;
    double s;
    int i, e;

    gs = this.m_m;
    d = this.m_dimension;
    max = (d / (gs << 1));

    s = 0d;
    e = 0;
    for (i = 0; i < max; i++) {
      s += Kernel.shiftedPermRosenbrock(x, this.m_o, this.m_p, e, gs); //
      e += gs;
    }

    return (s + Kernel.shiftedPermSphere(x, this.m_o, this.m_p, e,
        this.m_dimension - e));
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
    return "D/2m-group Shifted m-dimensional Rosenbrock���s Function";//$NON-NLS-1$
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
    return "F13"; //$NON-NLS-1$
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
    final int[] perm;
    final int max;

    d = this.m_opt;
    if (d != null) {
      return d;
    }

    d = ((double[]) (this.m_o.clone()));
    perm = this.m_p;
    max = (this.m_m*((this.m_dimension / (this.m_m << 1))));

    for (i = 0; i < max; i++) {
      d[perm[i]] += 1d;
    }
 

    return (this.m_opt = d);
  }
}
