/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * <p>
 * The Single-group Shifted and m-rotated Rastrigin���s Function: F5.
 * </p>
 * <p>
 * This function is not <warning>not threadsafe</warning> because it uses
 * internal temporary variables. Therefore, you should always use each
 * instance of this function for one single {#link java.lang.Thread} only.
 * You may clone or serialize function instances to use multiple threads.
 * <p>
 * 
 * @author Thomas Weise
 */
public final class F5 extends ShiftedPermutatedRotatedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 5d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F5(int dimension) {
	    this(Defaults.getRandomizer(F5.class), dimension);
  }

  /**
   * Create a new Single-group Shifted and m-rotated Rastrigin���s Function
   * 
   * @param o
   *          the shifted global optimum
   * @param p
   *          the permutation vector
   * @param m
   *          the rotation matrix
   */
  public F5(final double[] o, final int[] p, final double[] m) {
    super(o, p, m, MIN, MAX); 
  }

  /**
   * Create a default instance of F5.
   * 
   * @param r
   *          the randomizer to use
   */
  public F5(final Randomizer r) {
    this(r, Defaults.DEFAULT_DIM);//
  }

  public F5(final Randomizer r, int dimension) {
	    this(r.createShiftVector(dimension, MIN, MAX),//
	        r.createPermVector(dimension),//
	        r.createRotMatrix1D(Defaults.DEFAULT_M));//
	  }
  
  /**
   * Create a default instance of F5.
   */
  public F5() {
    this(Defaults.getRandomizer(F5.class));
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
    return (Kernel.shiftedPermRotRastrigin(x, this.m_o, this.m_p,//
        this.m_m, 0, this.m_matDim, this.m_tmp) * 1e6) + //
        Kernel.shiftedPermRastrigin(x, this.m_o, this.m_p, this.m_matDim,//
            this.m_dimension - this.m_matDim);
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
    return "Single-group Shifted and m-rotated Rastrigin���s Function";//$NON-NLS-1$
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
    return "F5"; //$NON-NLS-1$
  }
}
