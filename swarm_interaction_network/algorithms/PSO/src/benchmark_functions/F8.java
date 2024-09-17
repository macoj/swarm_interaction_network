/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
import java.io.BufferedWriter;

/**
 * The Single-group Shifted m-dimensional Rosenbrock���s Function: F8.
 *  
 * @author Thomas Weise
 */
public final class F8 extends ShiftedPermutatedFunction {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the maximum value */
  public static final double MAX = 100d;

  /** the minimum value */
  public static final double MIN = (-MAX);

  /** the m-value */
  private final int m_m;
  
  /**
   * Create a new function with a given dimensionality
   * 
   * @param dimension
   *          the dimension
   */
  public F8(int dimension) {
	    this(Defaults.getRandomizer(F8.class), dimension);
  }

  /** the optimum vector */
  private transient double[] m_opt;

  /**
   * Create a new Single-group Shifted m-dimensional Rosenbrock���s Function
   * 
   * @param o
   *          the shifted global optimum
   * @param p
   *          the permutation vector
   * @param m
   *          the fraction of nonseparability
   */
  public F8(final double[] o, final int[] p, final int m) {
    super(o, p, MIN, MAX);
    this.m_m = m;
  }

  // /**
  // * @return x
  // */
  // private static final double[] fakeo(){
  // double[] d;
  // int i;
  //    
  // d = new double[1000];
  // for(i=999;i>=0;i--){
  // d[i]=99;
  // }
  // return d;
  // }

  /**
   * Create a default instance of F8.
   * 
   * @param r
   *          the randomizer to use
   */
  public F8(final Randomizer r) {
    this(r, Defaults.DEFAULT_DIM);//
  }

  public F8(final Randomizer r, int dimension) {
    this(r.createShiftVector(dimension, MIN, MAX - 1d),//
        r.createPermVector(dimension),//
        Defaults.DEFAULT_M);//
  }
  /**
   * Create a default instance of F8.
   */
  public F8() {
    this(Defaults.getRandomizer(F8.class));
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
    return (Kernel.shiftedPermRosenbrock(x, this.m_o, this.m_p, //
        0, this.m_m) * 1e6) + //
        Kernel.shiftedPermSphere(x, this.m_o, this.m_p, this.m_m,//
            this.m_dimension - this.m_m);
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
    return "Single-group Shifted m-dimensional Rosenbrock���s Function";//$NON-NLS-1$
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
    return "F8"; //$NON-NLS-1$
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

    super.storeFunctionInfo(w);

    w.newLine();
    w.write("m-value  : ");//$NON-NLS-1$
    w.write(Integer.toString(this.m_m));
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

    d = this.m_opt;
    if (d != null) {
      return d;
    }

    d = ((double[]) (this.m_o.clone()));
    perm = this.m_p;
    for (i = (this.m_m - 1); i >= 0; i--) {
      d[perm[i]] += 1d;
    }

    return (this.m_opt = d);
  }
}