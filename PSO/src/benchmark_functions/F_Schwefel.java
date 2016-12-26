package benchmark_functions;

public final class F_Schwefel extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 500d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  public F_Schwefel() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Schwefel(int dimension) {
	    super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 420.96;
		}
	}
  
  public final double compute(final double[] x) {
	double sum = 0.0;
	for (int i = 0; i < this.m_dimension; i++) {
		double xi = x[i];
		sum += xi * Math.sin(Math.sqrt(Math.abs(xi)));
	}
	return 418.9829*((double) (this.m_dimension)) - sum; 
  }

  public final String getFullName() {
    return "Schwefel Function";
  }

  public final String getShortName() {
    return "Schwefel"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
