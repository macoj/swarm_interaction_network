/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * The shifted elliptic function: F1.
 * 
 * @author Thomas Weise
 */
public final class F1 extends ShiftedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 100d;

  /** the minimum value */
  public static final double MIN = (-MAX);

  /** the lookup table */
  private final double[] m_lookup;


  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F1(int dimension) {
	    this(Defaults.getRandomizer(F1.class).createShiftVector(dimension, MIN, MAX));
  }  
	  
  /**
   * Create a new shifted elliptic function
   * 
   * @param o
   *          the shifted global optimum
   */
  public F1(final double[] o) {
    super(o, MIN, MAX);
    this.m_lookup = Kernel.createPowLookup(this.m_dimension);
  }

  /**
   * Create a default instance of F1.
   * 
   * @param r
   *          the randomizer to use
   */
  public F1(final Randomizer r) {
    this(r.createShiftVector(Defaults.DEFAULT_DIM, MIN, MAX));
  }

  /**
   * Create a default instance of F1.
   */
  public F1() {
    this(Defaults.getRandomizer(F1.class));
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
    return Kernel.shiftedElliptic(x, this.m_o, 0, this.m_dimension,
        this.m_lookup);
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
    return "Shifted Elliptic Function";//$NON-NLS-1$
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
    return "F1"; //$NON-NLS-1$
  }
}
