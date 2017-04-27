package benchmark_functions;

public final class F_Rastrigin extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 5.12d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  public F_Rastrigin() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Rastrigin(int dimension) {
	    super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 0.0;
		}
	}
  
  public final double compute(final double[] x) {
	double result = 0;
	for (int i=0; i<this.m_dimension; i++) 
	{
		double xi = x[i];
		result += xi*xi - 10.0*Math.cos(2.0*Math.PI*xi);
	}
	return result + 10.0*((double) (this.m_dimension));  
  }

  public final String getFullName() {
    return "Rastrigin Function";
  }

  public final String getShortName() {
    return "Rastrigin"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
