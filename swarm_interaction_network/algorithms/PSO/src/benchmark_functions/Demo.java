package benchmark_functions;
import java.util.Random;

/*
 * Copyright (c) 2009 Thomas Weise
 * http://www.it-weise.de/
 * tweise@gmx.de
 * 
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */

/**
 * This is a demo class illustrating the usage of the benchmark functions.
 * For this purpose, we apply a very simple hill climbing strategy.
 */
public class Demo {

  /**
   * This is a demo class illustrating the usage of the benchmark
   * functions.
   * 
   * @param params
   *          the parameters which are completely ignored ^^
   */
  public static final void main(final String[] params) {
    final Function f;
    final double[] best, cur;
    final int dim;
    final double min, max;
    double b, c, old, neu;
    int i, it;
    Random r;

    // Create the default instance of F12.
    // Each objective function has been implemented in a class of a
    // corrsponding name.
    // Please be aware that you need one separate instance per thread!
    f = new F12();

    // Create a random number generator.
    r = new Random();

    // Obtain the dimension of the objective function: 1000 in the default
    // case.
    dim = f.getDimension();

    // Obtain the lower and upper bound of the search space.
    // This boundary is the same for all dimensions.
    min = f.getMin();
    max = f.getMax();

    best = new double[dim];
    cur = new double[dim];

    // All functions are subject to minimization, hence this is the worst
    // objective value possible.
    b = Double.POSITIVE_INFINITY;

    // Initialize the hill climber with a random point from the search
    // space.
    for (i = (dim - 1); i >= 0; i--) {
      cur[i] = (min + (r.nextDouble() * (max - min)));
    }

    // Run the hill climber infinitely long.
    for (it = 0;; it++) {
      // Compute the objective value.
      c = f.compute(cur);

      // Did we find an improvement?
      if (c < b) {
        b = c;
        System.arraycopy(cur, 0, best, 0, dim);
        System.out.println("new best: " + b + //$NON-NLS-1$
            " in iteration " + it);//$NON-NLS-1$
      }

      // A very simple mutation strategy.
      do {
        i = r.nextInt(dim);
        old = best[i];
        do {
          neu = (old + (r.nextGaussian() * Math.exp(r.nextInt(40) - 35)));
        } while ((neu < min) || (neu > max));
        cur[i] = neu;
      } while (r.nextInt(10) > 0);
    }

  }
}
