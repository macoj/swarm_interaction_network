/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * The shifted ackley's function: F3.
 * 
 * @author Thomas Weise
 */
public final class F3 extends ShiftedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 32d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F3(int dimension) {
	    this(Defaults.getRandomizer(F3.class).createShiftVector(dimension, MIN, MAX));
  }

  /**
   * Create a new shifted ackley's function
   * 
   * @param o
   *          the shifted global optimum
   */
  public F3(final double[] o) {
    super(o, MIN, MAX);
  }

  /**
   * Create a default instance of F3.
   * 
   * @param r
   *          the randomizer to use
   */
  public F3(final Randomizer r) {
    this(r.createShiftVector(Defaults.DEFAULT_DIM, MIN, MAX));//
  }

  /**
   * Create a default instance of F3.
   */
  public F3() {
    this(Defaults.getRandomizer(F3.class));
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
    return Kernel.shiftedAckley(x, this.m_o, 0, this.m_dimension);
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
    return "Shifted Ackley���s Function";//$NON-NLS-1$
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
    return "F3"; //$NON-NLS-1$
  }
  
  /**
   * The test function used to check whether the routines here have been
   * implemented correctly.
   * 
   * @param params
   *          the parameters
   */
  public static final void main(final String[] params) {
        
    System.out.println(
    -20d * 1 + 20d - Math.exp(1)  + Math.E);
    
//    System.out.println(Kernel.protAckley(0d,100,1000));
  }
}
