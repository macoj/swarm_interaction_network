package benchmark_functions;

public final class F_Rosenbrock extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 2.048d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  public F_Rosenbrock() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Rosenbrock(int dimension) {
		super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 1.0;
		}
	}
  
  public final double compute(final double[] x) {
	double result_value = 0;
	for (int i=0; i<this.m_dimension-1; i++) 
	{
		double xi = x[i];
		double xi_1 = x[i+1];
		result_value += 100.0*Math.pow(xi_1 - xi*xi, 2.0) + Math.pow(1.0 - xi, 2.0);
	}	
	return result_value;
  }

  public final String getFullName() {
    return "Rosenbrock Function";
  }

  public final String getShortName() {
    return "Rosenbrock"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
