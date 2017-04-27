/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
/**
 * The internal mathematical kernel of the benchmark suite.
 * 
 * @author Thomas Weise
 */
final class Kernel {

  /** two pi */
  private static final double TPI = (Math.PI + Math.PI);

  /**
   * Multiply a vector with a matrix.
   * 
   * @param vec
   *          the vector
   * @param mat
   *          the matrix
   * @param dest
   *          the destination vector
   */
  final static void multiply(final double[] vec, final double[] mat,
      final double[] dest) {
    final int count, upper, max;
    int i, j, k;
    double s;

    count = dest.length;
    upper = (count - 1);
    max = (count * upper);
    for (i = upper; i >= 0; i--) {
      s = 0d;
      for (k = upper, j = max + i; k >= 0; k--, j -= count) {
        s += (mat[j] * vec[k]);
      }

      dest[i] = s;
    }
  }

  /**
   * Compute the shifted sphere function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedSphere(final double[] x, final double[] o,
      final int start, final int count) {
    int j;
    double s, z;

    s = 0d;
    j = (count - 1 + start);
    for (; j >= start; j--) {
      z = (x[j] - o[j]);
      s += (z * z);
    }

    return s;
  }

  /**
   * Compute the shifted permutated sphere function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedPermSphere(final double[] x,
      final double[] o, final int[] P, final int start, final int count) {
    int j, k;
    double s, z;

    s = 0d;
    j = (count - 1 + start);
    for (; j >= start; j--) {
      k = P[j];
      z = (x[k] - o[k]);
      s += (z * z);
    }

    return s;
  }

  // /**
  // * Compute the shifted elliptic function
  // *
  // * @param x
  // * the input vector
  // * @param o
  // * the global optimum
  // * @param start
  // * the start index
  // * @param count
  // * the number of elements to consider in the computation
  // * @return the result
  // */
  // static final double shiftedElliptic(final double[] x, final double[]
  // o,
  // final int start, final int count) {
  // int i, j;
  // final double pow;
  // double s, z;
  //
  // s = 0d;
  // i = (count - 1);
  // j = (i + start);
  // pow = (1d / i);
  // for (; i >= 0; i--, j--) {
  // z = (x[j] - o[j]);
  // s += (Math.pow(1e6, i * pow) * z * z);
  // }
  //
  // return s;
  // }

  /**
   * Create a lookup table for the pow function
   * 
   * @param count
   *          the counter
   * @return the table
   */
  static final double[] createPowLookup(final int count) {
    int i;
    final double[] d;
    final double pow;

    d = new double[count];
    i = (count - 1);
    pow = (1d / i);
    d[i] = 1e6d;
    d[0] = 1d;

    for (--i; i > 0; i--) {
      d[i] = Math.pow(1e6d, i * pow);
    }

    return d;
  }

  /**
   * Compute the shifted elliptic function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @param lookup
   *          the lookup table
   * @return the result
   */
  static final double shiftedElliptic(final double[] x, final double[] o,
      final int start, final int count, final double[] lookup) {
    int i, j;
    // final double pow;
    double s, z;

    s = 0d;
    i = (count - 1);
    j = (i + start);
    // pow = (1d / i);
    for (; i >= 0; i--, j--) {
      z = (x[j] - o[j]);
      s += (/* Math.pow(1e6, i pow) */lookup[i] * z * z);
    }

    return s;
  }

  /**
   * Compute the shifted permutated elliptic function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @param lookup
   *          the lookup table
   * @return the result
   */
  static final double shiftedPermElliptic(final double[] x,
      final double[] o, final int[] P, final int start, final int count,
      final double[] lookup) {
    int i, j, k;
    // final double pow;
    double s, z;

    s = 0d;
    i = (count - 1);
    j = (i + start);
    // pow = (1d / i);
    for (; i >= 0; i--, j--) {
      k = P[j];
      z = (x[k] - o[k]);
      s += (/* Math.pow(1e6, i pow) */lookup[i] * z * z);
    }

    return s;
  }

  /**
   * Compute the shifted, permutated, and rotated elliptic function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param M
   *          the rotation matrix
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @param z
   *          a temporary array
   * @param lookup
   *          the lookup table
   * @return the result
   */
  static final double shiftedPermRotElliptic(final double[] x,
      final double[] o, final int[] P, final double[] M, final int start,
      final int count, final double[] z, final double[] lookup) {
    final int upper, max;
    int i, j, k;
    // final double pow;
    double rz, s;

    // compute z
    i = upper = (count - 1);
    j = (i + start);
    for (; i >= 0; i--, j--) {
      k = P[j];
      z[i] = (x[k] - o[k]);
    }

    // rotate and compute function at the same time:
    max = (count * upper);
    // pow = (1d / upper);
    s = 0d;
    for (i = upper; i >= 0; i--) {

      // rotate
      rz = 0d;
      for (k = upper, j = max + i; k >= 0; k--, j -= count) {
        rz += (M[j] * z[k]);
      }

      // compute function
      s += (/* Math.pow(1e6, i pow) */lookup[i] * rz * rz);
    }

    return s;
  }

  /**
   * Compute the shifted rastrigin's function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedRastrigin(final double[] x, final double[] o,
      final int start, final int count) {
    int j;
    double s, z;

    s = 0d;
    j = (count - 1 + start);
    for (; j >= start; j--) {
      z = (x[j] - o[j]);
      s += ((z * z) - (10d * Math.cos(TPI * z)) + 10d);
    }

    return s;
  }

  /**
   * Compute the shifted permutated rastrigin's function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedPermRastrigin(final double[] x,
      final double[] o, final int[] P, final int start, final int count) {
    int j, k;
    double s, z;

    s = 0d;
    j = (count - 1 + start);
    for (; j >= start; j--) {
      k = P[j];
      z = (x[k] - o[k]);
      s += ((z * z) - (10d * Math.cos(TPI * z)) + 10d);
    }

    return s;
  }

  /**
   * Compute the shifted, permutated, and rotated rastrigin function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param M
   *          the rotation matrix
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @param z
   *          a temporary array
   * @return the result
   */
  static final double shiftedPermRotRastrigin(final double[] x,
      final double[] o, final int[] P, final double[] M, final int start,
      final int count, final double[] z) {
    final int upper, max;
    int i, j, k;
    double rz, s;

    // compute z
    i = upper = (count - 1);
    j = (i + start);
    for (; i >= 0; i--, j--) {
      k = P[j];
      z[i] = (x[k] - o[k]);
    }

    // rotate and compute function at the same time:
    max = (count * upper);
    s = 0d;
    for (i = upper; i >= 0; i--) {

      // rotate
      rz = 0d;
      for (k = upper, j = max + i; k >= 0; k--, j -= count) {
        rz += (M[j] * z[k]);
      }

      // compute function
      s += ((rz * rz) - (10d * Math.cos(TPI * rz)) + 10d);
    }

    return s;
  }

  /**
   * Compute the shifted ackley's function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedAckley(final double[] x, final double[] o,
      final int start, final int count) {
    int j;
    double s1, s2, z;

    s1 = 0d;
    s2 = 0d;
    j = (count - 1 + start);
    for (; j >= start; j--) {
      z = (x[j] - o[j]);
      s1 += (z * z);
      s2 += Math.cos(TPI * z);
    }

    return /*
            * (-20d Math.exp(-0.2d Math.sqrt(s1 / count))) - // Math.exp(s2
            * / count) + 20d + Math.E;
            */
    protAckley(s1, s2, count);
  }

  /**
   * The protected ackley's function
   * 
   * @param s1
   *          the first sum
   * @param s2
   *          the second sum
   * @param count
   *          the counter
   * @return the result
   */
  private static final double protAckley(final double s1, final double s2,
      final double count) {
    final double e1, e2;

    if (s1 == 0d) {
      e1 = 1d;
    } else {
      e1 = Math.exp(-0.2d * Math.sqrt(s1 / count));
    }

    if (s2 == count) {
      e2 = Math.E;
    } else {
      e2 = Math.exp(s2 / count);
    }

    return ((20d - (20d * e1)) + (Math.E - e2));
  }

  /**
   * Compute the shifted permutated ackley's function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedPermAckley(final double[] x,
      final double[] o, final int[] P, final int start, final int count) {
    int j, k;
    double s1, s2, z;

    s1 = 0d;
    s2 = 0d;
    j = (count - 1 + start);
    for (; j >= start; j--) {
      k = P[j];
      z = (x[k] - o[k]);
      s1 += (z * z);
      s2 += Math.cos(TPI * z);
    }

    return /*
            * (-20d Math.exp(-0.2d Math.sqrt(s1 / count))) - // Math.exp(s2
            * / count) + 20d + Math.E;
            */
    protAckley(s1, s2, count);
  }

  /**
   * Compute the shifted, permutated, and rotated ackley's function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param M
   *          the rotation matrix
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @param z
   *          a temporary array
   * @return the result
   */
  static final double shiftedPermRotAckley(final double[] x,
      final double[] o, final int[] P, final double[] M, final int start,
      final int count, final double[] z) {
    final int upper, max;
    int i, j, k;
    double rz, s1, s2;

    // compute z
    i = upper = (count - 1);
    j = (i + start);
    for (; i >= 0; i--, j--) {
      k = P[j];
      z[i] = (x[k] - o[k]);
    }

    // rotate and compute function at the same time:
    max = (count * upper);
    s1 = 0d;
    s2 = 0d;
    for (i = upper; i >= 0; i--) {

      // rotate
      rz = 0d;
      for (k = upper, j = max + i; k >= 0; k--, j -= count) {
        rz += (M[j] * z[k]);
      }

      // compute function
      s1 += (rz * rz);
      s2 += Math.cos(TPI * rz);
    }

    return /*
            * (-20d Math.exp(-0.2d Math.sqrt(s1 / count))) - // Math.exp(s2
            * / count) + 20d + Math.E;
            */
    protAckley(s1, s2, count);
  }

  /**
   * Compute the shifted version of schwefel's problem 1.2
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedSchwefel12(final double[] x,
      final double[] o, final int start, final int count) {
    int j;
    final int end;
    double s1, s2, z;

    s1 = 0d;
    s2 = 0d;
    end = (start + count);
    for (j = start; j < end; j++) {
      z = (x[j] - o[j]);
      s1 += z;
      s2 += (s1 * s1);
    }

    return s2;
  }

  /**
   * Compute the shifted permutated schwefel's problem 1.2
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedPermSchwefel12(final double[] x,
      final double[] o, final int[] P, final int start, final int count) {
    int j, k;
    final int end;
    double s1, s2, z;

    s1 = 0d;
    s2 = 0d;
    end = (start + count);
    for (j = start; j < end; j++) {
      k = P[j];
      z = (x[k] - o[k]);
      s1 += z;
      s2 += (s1 * s1);
    }

    return s2;
  }

  // upsi, schwefel's problem 1.2 is already nonseparable, so we do not
  // need to rotate it :-D
  //
  // /**
  // * Compute the shifted, permutated, and rotated schwefel's problem 1.2
  // *
  // * @param x
  // * the input vector
  // * @param o
  // * the global optimum
  // * @param P
  // * the permutation
  // * @param M
  // * the rotation matrix
  // * @param start
  // * the start index
  // * @param count
  // * the number of elements to consider in the computation
  // * @param z
  // * a temporary array
  // * @return the result
  // */
  // static final double shiftedPermRotSchwefel12(final double[] x,
  // final double[] o, final int[] P, final double[] M, final int start,
  // final int count, final double[] z) {
  // final int upper, max;
  // int i, j, k;
  // double rz, s1, s2;
  //
  // // compute z
  // i = upper = (count - 1);
  // j = (i + start);
  // for (; i >= 0; i--, j--) {
  // k = P[j];
  // z[i] = (x[k] - o[k]);
  // }
  //
  // // rotate and compute function at the same time:
  // max = (count * upper);
  // s1 = 0d;
  // s2 = 0d;
  // for (i = 0; i < count; i++) {
  //
  // // rotate
  // rz = 0d;
  // for (k = upper, j = max + i; k >= 0; k--, j -= count) {
  // rz += (M[j] * z[k]);
  // }
  //
  // // compute function
  // s1 += rz;
  // s2 += (s1 * s1);
  // }
  //
  // return s2;
  // }

  /**
   * Compute the shifted version of rosenbrock's function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedRosenbrock(final double[] x,
      final double[] o, final int start, final int count) {
    int j;
    double s, z, oz, t;

    j = (count + start - 1);
    z = (x[j] - o[j]);
    s = 0d;

    for (--j; j >= start; j--) {
      oz = z;
      z = (x[j] - o[j]);

      t = ((z * z) - oz);
      s += (100d * t * t);

      t = (z - 1d);
      s += (t * t);
    }

    return s;
  }

  /**
   * Compute the shifted permutated rosenbrock's function
   * 
   * @param x
   *          the input vector
   * @param o
   *          the global optimum
   * @param P
   *          the permutation
   * @param start
   *          the start index
   * @param count
   *          the number of elements to consider in the computation
   * @return the result
   */
  static final double shiftedPermRosenbrock(final double[] x,
      final double[] o, final int[] P, final int start, final int count) {
    int j, k;
    double s, z, oz, t;

    s = 0d;
    j = (count + start - 1);

    k = P[j];
    z = (x[k] - o[k]);

    for (--j; j >= start; j--) {
      oz = z;
      k = P[j];
      z = (x[k] - o[k]);

      t = ((z * z) - oz);
      s += (100d * t * t);
      t = (z - 1d);
      s += (t * t);
    }

    return s;
  }

  // /**
  // * Compute the bth power of a
  // *
  // * @param a
  // * the base
  // * @param b
  // * the exponent
  // * @return a^b
  // */
  // private static final double pow(final double a, final double b) {
  // return StrictMath.exp(b * StrictMath.log(a));
  // }

  /**
   * The test function used to check whether the routines here have been
   * implemented correctly.
   * 
   * @param params
   *          the parameters
   */
  public static final void main(final String[] params) {
    // double[] x, m, o, z;
    // int[] p;
    //
    // x = new double[] { 1d, 2d, 3d };
    // m = new double[] { 3, -4, 1, 5, 3, -7, -9, 2, 6 };
    // o = new double[] { 7d, 5d, 2d };
    // p = new int[] { 0, 1, 2 };
    // z = new double[3];
    //
    // System.out.println(shiftedPermRotElliptic(x, o, p, m, 0, 3, z));

    long t;
    int i;

    // t = -System.nanoTime();
    // for(i=10000;i>=0;i--){
    // pow(2,3);
    // }
    // System.out.println(System.nanoTime()-t);

    t = -System.nanoTime();
    for (i = 10000; i >= 0; i--) {
      Math.pow(2, 3);
    }
    System.out.println(System.nanoTime() - t);
  }

}
