package benchmark_functions;

import java.util.Random;

public final class F_Random extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 1;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  private Random random_number_generator;
  
  public F_Random() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Random(int dimension) {
	    super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 0.0;
		}
		random_number_generator = new Random(System.currentTimeMillis());
	}
  
  public final double compute(final double[] x) {
	return random_number_generator.nextDouble();  
  }

  public final String getFullName() {
    return "Random Function";
  }

  public final String getShortName() {
    return "Random"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
