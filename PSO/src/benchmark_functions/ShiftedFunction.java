/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
import java.io.File;

/**
 * A shifted function
 * 
 * @author Thomas Weise
 */
public abstract class ShiftedFunction extends Function {

  /** the serial version id */
  private static final long serialVersionUID = 1;

  /** the shifted global optimum */
  protected final double[] m_o;

  /**
   * Create a new shifted function
   * 
   * @param o
   *          the shifted global optimum
   * @param min
   *          the minimum value a decision variable can take on
   * @param max
   *          the maximum value a decision variable can take on
   */
  protected ShiftedFunction(final double[] o, final double min,
      final double max) {
    super(o.length, min, max);
    this.m_o = o;
  }

  /**
   * Obtain the optimum vector of the benchmark function
   * 
   * @return the optimum vector of the benchmark function
   */
  // @Override
  public double[] getOptimum() {
    return this.m_o;
  }

  /**
   * Obtain the shift vector of the benchmark function
   * 
   * @return the shift vector of the benchmark function
   */
  public final double[] getShiftVector() {
    return this.m_o;
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
      Utils.storeShiftVector(new File(f, this.getShortName() + //
          "_o.txt").toString(), this.m_o); //$NON-NLS-1$
    } catch (Throwable t) {
      if (s == null) {
        s = t;
      }
    }

    if (s != null) {
      throw s;
    }
  }
}
