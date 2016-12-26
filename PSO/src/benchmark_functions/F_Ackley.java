package benchmark_functions;

public final class F_Ackley extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 32.768d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  public F_Ackley() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Ackley(int dimension) {
	  	super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 0.0;
		}
	}
  
  public final double compute(final double[] x) {
	double sum_1 = 0.0;
	double sum_2 = 0.0;
	for (int i=0; i<this.m_dimension; i++) 
	{
		double xi = x[i];
		sum_1 += xi * xi;
		sum_2 += Math.cos(2.0*Math.PI*xi);

	}
	return -20.0*Math.exp(-0.2*Math.sqrt(sum_1/this.m_dimension))-Math.exp(sum_2/this.m_dimension)+20.0+Math.E; 
  }

  public final String getFullName() {
    return "Ackley Function";
  }

  public final String getShortName() {
    return "Ackley"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
