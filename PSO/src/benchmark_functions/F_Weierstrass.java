package benchmark_functions;

public final class F_Weierstrass extends Function {

  private static final long serialVersionUID = 1L;

  /** the maximum value */
  public static final double MAX = 0.5d;

  /** the minimum value */
  public static final double MIN = (-MAX);
  
  public double[] optimum;

  public F_Weierstrass() {
	  this(Defaults.DEFAULT_DIM);
  }
  
  public F_Weierstrass(int dimension) {
	    super(dimension, MIN, MAX);
		optimum = new double[dimension];
		for (int i = 0; i < optimum.length; i++) {
			optimum[i] = 0.0;
		}
	}
  
  public final double compute(final double[] x) {
	int k_max = 20;
	double a = 0.5;
	double b = 3.0;
	double sum1 = 0.0;
	for (int i = 0 ; i < this.m_dimension ; i ++) {
		for (int k = 0 ; k <= k_max ; k ++) {
			sum1 += Math.pow(a, k) * Math.cos(Math.PI * 2.0 * Math.pow(b, k) * (x[i] + 0.5));
		}
	}

	double sum2 = 0.0;
	for (int k = 0 ; k <= k_max ; k ++) {
		sum2 += Math.pow(a, k) * Math.cos(Math.PI * 2.0 * Math.pow(b, k) * (0.5));
	}

	return (sum1 - sum2*((double )(this.m_dimension)));  
  }

  public final String getFullName() {
    return "Weierstrass Function";
  }

  public final String getShortName() {
    return "Weierstrass"; 
  }

  public double[] getOptimum() {
	return optimum;
  }
}
