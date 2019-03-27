package benchmark_functions;

public final class F_Griewank extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 600;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  public F_Griewank() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Griewank(int dimension) {
	    super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 0.0;
		}
	}
  
  public final double compute(final double[] x) {
	double sum = 0.0;
	double prod = 1.0;
	for (int i = 1; i <= this.m_dimension; i++) {
		double xi = x[(i - 1)];
		sum += (xi * xi);
		prod *= Math.cos(xi / Math.sqrt(i));
	}
	return (sum / 4000.0) - prod + 1.0;
  }

  public final String getFullName() {
    return "Griewank Function";
  }

  public final String getShortName() {
    return "Griewank"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
