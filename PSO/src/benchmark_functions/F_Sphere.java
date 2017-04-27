package benchmark_functions;

public final class F_Sphere extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 100d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  public F_Sphere() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Sphere(int dimension) {
	    super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 0.0;
		}
	}
  
  public final double compute(final double[] x) {
    double sum_of_squares = 0.0;
	for (int i = 0; i < this.m_dimension; i++) {
		sum_of_squares += (x[i]*x[i]);
	}
	return sum_of_squares;  
  }

  public final String getFullName() {
    return "Sphere Function";
  }

  public final String getShortName() {
    return "Sphere"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
