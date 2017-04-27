/*
 * Copyright (c) 2009 Thomas Weise for NICAL
 * http://www.it-weise.de/
 * tweise@gmx.de
 *
 * GNU LESSER GENERAL PUBLIC LICENSE (Version 2.1, February 1999)
 */
package benchmark_functions;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URL;
import java.util.regex.Pattern;

/**
 * With this class, we provide utility functions such as the creation,
 * loading, or storing of vectors and matrices.
 * 
 * @author Thomas Weise
 */
public final class Utils {

  /**
   * Load a shifting vector from a file
   * 
   * @param file
   *          the path to the file to load
   * @return the vector loaded
   */
  public static final double[] loadShiftVector(final String file) {
    FileInputStream r;
    try {
      r = new FileInputStream(file);
      try {
        return loadShiftVector(r);
      } catch (Throwable t1) {
        try {
          r.close();
        } catch (Throwable t2) {
          // not interesting here
        }
        throw t1;
      }

    } catch (RuntimeException t1) {
      throw t1;
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a shifting vector from an uri
   * 
   * @param uri
   *          the uri to load from
   * @return the vector loaded
   */
  public static final double[] loadShiftVector(final URI uri) {
    URL u;
    InputStream r;
    try {
      u = uri.toURL();
      r = u.openStream();
      try {
        return loadShiftVector(r);
      } catch (Throwable t1) {
        try {
          r.close();
        } catch (Throwable t2) {
          // not interesting here
        }
        throw t1;
      }

    } catch (RuntimeException t1) {
      throw t1;
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a shifting vector from an input stream
   * 
   * @param in
   *          the input stream to load from
   * @return the vector loaded
   */
  public static final double[] loadShiftVector(final InputStream in) {
    BufferedReader r;
    double[] buf, tmp;
    int fill;
    String s;
    double d;

    try {
      buf = new double[1000];
      fill = 0;
      r = new BufferedReader(new InputStreamReader(in));

      try {

        while ((s = r.readLine()) != null) {
          s = s.trim();
          if (s.length() <= 0) {
            continue;
          }

          d = Double.parseDouble(s);
          if (fill >= buf.length) {
            tmp = new double[fill << 1];
            System.arraycopy(buf, 0, tmp, 0, fill);
            buf = tmp;
          }

          buf[fill++] = d;
        }

      } finally {
        r.close();
      }

      if (fill == buf.length) {
        return buf;
      }

      tmp = new double[fill];
      System.arraycopy(buf, 0, tmp, 0, fill);
      return tmp;

    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Store a shifting vector to a file
   * 
   * @param file
   *          the path to the file to load
   * @param o
   *          the vector to be stored
   */
  public static final void storeShiftVector(final String file,
      final double[] o) {
    BufferedWriter w;
    int i;

    try {
      w = new BufferedWriter(new FileWriter(file));
      try {

        for (i = 0; i < o.length; i++) {
          if (i > 0) {
            w.newLine();
          }
          w.write(Double.toString(o[i]));
          if ((i % 100) <= 0) {
            w.flush();
          }
        }

      } finally {
        w.close();
      }
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a permutation vector from a file
   * 
   * @param file
   *          the path to the file to load
   * @return the vector loaded
   */
  public static final int[] loadPermVector(final String file) {
    FileInputStream r;
    try {
      r = new FileInputStream(file);
      try {
        return loadPermVector(r);
      } catch (Throwable t1) {
        try {
          r.close();
        } catch (Throwable t2) {
          // not interesting here
        }
        throw t1;
      }

    } catch (RuntimeException t1) {
      throw t1;
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a permutation vector from an uri
   * 
   * @param uri
   *          the uri to load from
   * @return the vector loaded
   */
  public static final int[] loadPermVector(final URI uri) {
    URL u;
    InputStream r;
    try {
      u = uri.toURL();
      r = u.openStream();
      try {
        return loadPermVector(r);
      } catch (Throwable t1) {
        try {
          r.close();
        } catch (Throwable t2) {
          // not interesting here
        }
        throw t1;
      }

    } catch (RuntimeException t1) {
      throw t1;
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a permutation vector from an input stream
   * 
   * @param in
   *          the input stream
   * @return the vector loaded
   */
  public static final int[] loadPermVector(final InputStream in) {
    BufferedReader r;
    int[] buf, tmp;
    int fill, d;
    String s;

    try {
      buf = new int[1000];
      fill = 0;
      r = new BufferedReader(new InputStreamReader(in));

      try {

        while ((s = r.readLine()) != null) {
          s = s.trim();
          if (s.length() <= 0) {
            continue;
          }

          d = Integer.parseInt(s);
          if (fill >= buf.length) {
            tmp = new int[fill << 1];
            System.arraycopy(buf, 0, tmp, 0, fill);
            buf = tmp;
          }

          buf[fill++] = (d - 1);
        }

      } finally {
        r.close();
      }

      if (fill == buf.length) {
        return buf;
      }

      tmp = new int[fill];
      System.arraycopy(buf, 0, tmp, 0, fill);
      return tmp;

    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Store a permutation vector to a file
   * 
   * @param file
   *          the path to the file to load
   * @param P
   *          the vector to be stored
   */
  public static final void storePermVector(final String file, final int[] P) {
    BufferedWriter w;
    int i;

    try {
      w = new BufferedWriter(new FileWriter(file));
      try {

        for (i = 0; i < P.length; i++) {
          if (i > 0) {
            w.newLine();
          }
          w.write(Integer.toString(P[i] + 1));
          if ((i % 100) <= 0) {
            w.flush();
          }
        }

      } finally {
        w.close();
      }
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a rotation matrix from a file
   * 
   * @param file
   *          the path to the file to load
   * @return the matrix loaded
   */
  public static final double[] loadRotMatrix(final String file) {
    FileInputStream r;
    try {
      r = new FileInputStream(file);
      try {
        return loadRotMatrix(r);
      } catch (Throwable t1) {
        try {
          r.close();
        } catch (Throwable t2) {
          // not interesting here
        }
        throw t1;
      }

    } catch (RuntimeException t1) {
      throw t1;
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a rotation matrix from an uri
   * 
   * @param uri
   *          the uri to load from
   * @return the matrix loaded
   */
  public static final double[] loadRotMatrix(final URI uri) {
    URL u;
    InputStream r;
    try {
      u = uri.toURL();
      r = u.openStream();
      try {
        return loadRotMatrix(r);
      } catch (Throwable t1) {
        try {
          r.close();
        } catch (Throwable t2) {
          // not interesting here
        }
        throw t1;
      }

    } catch (RuntimeException t1) {
      throw t1;
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Load a rotation matrix from an input stream
   * 
   * @param in
   *          the input stream
   * @return the matrix loaded
   */
  public static final double[] loadRotMatrix(final InputStream in) {
    BufferedReader r;
    double[] buf, tmp;
    int fill, i, l;
    double d;
    String[] ss;
    String s;
    final Pattern p;

    p = Pattern.compile(" "); //$NON-NLS-1$

    try {
      buf = new double[1000];
      fill = 0;
      r = new BufferedReader(new InputStreamReader(in));

      try {

        outer: while ((s = r.readLine()) != null) {
          s = s.trim();

          if (s.length() <= 0) {
            continue outer;
          }

          s = s.replace('\t', ' ').replace('\f', ' ').replace('\b', ' ');
          ss = p.split(s);

          l = ss.length;
          if (l <= 0) {
            continue outer;
          }

          inner: for (i = 0; i < l; i++) {
            s = ss[i].trim();
            if (s.length() <= 0) {
              continue inner;
            }

            d = Double.parseDouble(s);
            if (fill >= buf.length) {
              tmp = new double[fill << 1];
              System.arraycopy(buf, 0, tmp, 0, fill);
              buf = tmp;
            }

            buf[fill++] = d;
          }
        }

      } finally {
        r.close();
      }

      if (fill == buf.length) {
        return buf;
      }

      tmp = new double[fill];
      System.arraycopy(buf, 0, tmp, 0, fill);
      return tmp;

    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  /**
   * Store a rotation matrix to a file
   * 
   * @param file
   *          the path to the file to load
   * @param o
   *          the matrix to be stored
   */
  public static final void storeRotMatrix(final String file,
      final double[] o) {
    BufferedWriter w;
    final int len;
    int i;

    try {
      w = new BufferedWriter(new FileWriter(file));
      try {
        len = ((int) (0.9d + Math.sqrt(o.length)));

        for (i = 0; i < o.length; i++) {
          if (i > 0) {
            if ((i % len) <= 0) {
              w.newLine();
            } else {
              w.write(' ');
            }
          }
          w.write(Double.toString(o[i]));
          if ((i % 100) <= 0) {
            w.flush();
          }
        }

      } finally {
        w.close();
      }
    } catch (Throwable t) {
      throw new RuntimeException(t); // bypass exception checking
    }
  }

  // /** the default seed for randomized function parameters */
  // protected static final long DEFAULT_SEED = 1;

  // /**
  // * Create a random permutation vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @return the permutation vector
  // */
  // public static final int[] defaultRandomPermVector(final int dim) {
  // return createRandomPermVector(dim, DEFAULT_SEED);
  // }
  //
  // /**
  // * Create a random permutation vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @return the permutation vector
  // */
  // public static final int[] createRandomPermVector(final int dim) {
  // return createRandomPermVector(dim, new Random());
  // }
  //
  // /**
  // * Create a random permutation vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @param seed
  // * the random seed
  // * @return the permutation vector
  // */
  // public static final int[] createRandomPermVector(final int dim,
  // final long seed) {
  // return createRandomPermVector(dim, new Random(seed));
  // }

  // /**
  // * Create a random permutation vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @param r
  // * the randomizer
  // * @return the permutation vector
  // */
  // public static final int[] createRandomPermVector(final int dim,
  // final Random r) {
  // final int[] d;
  // int i, j, k, t;
  //
  // d = new int[dim];
  //
  // for (i = (dim - 1); i >= 0; i--) {
  // d[i] = i;
  // }
  //
  // synchronized (r) {
  // for (i = (dim << 3); i >= 0; i--) {
  // j = r.nextInt(dim);
  // do {
  // k = r.nextInt(dim);
  // } while (k == j);
  //
  // t = d[j];
  // d[j] = d[k];
  // d[k] = t;
  // }
  // }
  //
  // return d;
  // }
  //
  // /**
  // * Create a random rotation matrix
  // *
  // * @param dim
  // * the number of dimensions
  // * @return the rotation matrix
  // */
  // public static final double[] createRandomRotMatrix(final int dim) {
  // return createRandomRotMatrix(dim, new Random());
  // }
  //
  // /**
  // * Obtain the default rotation matrix
  // *
  // * @param dim
  // * the number of dimensions
  // * @return the default rotation matrix
  // */
  // public static final double[] defaultRandomRotMatrix(final int dim) {
  // return createRandomRotMatrix(dim, DEFAULT_SEED);
  // }
  //
  // /**
  // * Create a random rotation matrix
  // *
  // * @param dim
  // * the number of dimensions
  // * @param seed
  // * the random seed
  // * @return the rotation matrix
  // */
  // public static final double[] createRandomRotMatrix(final int dim,
  // final long seed) {
  // return createRandomRotMatrix(dim, new Random(seed));
  // }
  //
  // /**
  // * Create a random rotation matrix
  // *
  // * @param dim
  // * the number of dimensions
  // * @param r
  // * the randomizer
  // * @return the rotation matrix
  // */
  // public static final double[] createRandomRotMatrix(final int dim,
  // final Random r) {
  // final double[] m;
  // int i, j, k;
  // double dp, t;
  //
  // m = new double[dim * dim];
  //
  // synchronized (r) {
  // outer: for (;;) {
  //
  // // initialize
  // for (i = (m.length - 1); i >= 0; i--) {
  // m[i] = r.nextDouble();
  // }
  //
  // // main loop
  // for (i = (dim - 1); i >= 0; i--) {
  //
  // //
  // for (j = (dim - 1); j > i; j--) {
  //
  // // dot product
  // dp = 0d;
  // for (k = (dim - 1); k >= 0; k--) {
  // dp += (m[(i * dim) + k] * m[(j * dim) + k]);
  // }
  //
  // // subtract
  // for (k = (dim - 1); k >= 0; k--) {
  // m[(i * dim) + k] -= (dp * m[(j * dim) + k]);
  // }
  // }
  //
  // // normalize
  // dp = 0d;// subtract
  // for (k = (dim - 1); k >= 0; k--) {
  // t = m[(i * dim) + k];
  // dp += (t * t);
  // }
  //
  // // linear dependency -> restart
  // if (dp <= 0d) {
  // continue outer;
  // }
  // dp = (1d / Math.sqrt(dp));
  //
  // for (k = (dim - 1); k >= 0; k--) {
  // m[(i * dim) + k] *= dp;
  // }
  // }
  //
  // return m;
  // }
  // }
  // }
  //
  // /**
  // * Create the default shift vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @param min
  // * the minimum value a decision variable can take on
  // * @param max
  // * the maximum value a decision variable can take on
  // * @return the shift vector
  // */
  // public static final double[] defaultRandomShiftVector(final int dim,
  // final double min, final double max) {
  // return createRandomShiftVector(dim, min, max, DEFAULT_SEED);
  // }
  //
  // /**
  // * Create a random shift vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @param min
  // * the minimum value a decision variable can take on
  // * @param max
  // * the maximum value a decision variable can take on
  // * @param seed
  // * the randomizer seed
  // * @return the shift vector
  // */
  // public static final double[] createRandomShiftVector(final int dim,
  // final double min, final double max, final long seed) {
  // return createRandomShiftVector(dim, min, max, new Random(seed));
  // }
  //
  // /**
  // * Create a random shift vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @param min
  // * the minimum value a decision variable can take on
  // * @param max
  // * the maximum value a decision variable can take on
  // * @return the shift vector
  // */
  // public static final double[] createRandomShiftVector(final int dim,
  // final double min, final double max) {
  // return createRandomShiftVector(dim, min, max, new Random());
  // }
  //
  // /**
  // * Create a random shift vector
  // *
  // * @param dim
  // * the number of dimensions
  // * @param min
  // * the minimum value a decision variable can take on
  // * @param max
  // * the maximum value a decision variable can take on
  // * @param r
  // * the randomizer
  // * @return the shift vector
  // */
  // public static final double[] createRandomShiftVector(final int dim,
  // final double min, final double max, final Random r) {
  // final double[] d;
  // final double w, hw, middle;
  // int i;
  //
  // w = (max - min);
  // hw = (0.5d * w);
  // middle = (min + hw);
  // d = new double[dim];
  //
  // synchronized (r) {
  // for (i = (dim - 1); i >= 0; i--) {
  // d[i] = (middle + (r.nextDouble() * w) - hw);
  // }
  // }
  //
  // return d;
  // }

  /**
   * The test function used to check whether the routines here have been
   * implemented correctly.
   * 
   * @param params
   *          the parameters
   */
  public static final void main(final String[] params) {
    double[] m;
    int n, i, j;

    n = 50;
    //storeRotMatrix("E:\\1.txt", createRandomRotMatrix(n)); //$NON-NLS-1$
    //m = loadRotMatrix("E:\\1.txt"); //$NON-NLS-1$
    m = loadRotMatrix("e:/ackley_M_D50.txt"); //$NON-NLS-1$
    // m = createRandomRotMatrix(n);
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        if (j > 0) {
          System.out.print('\t');
        }
        System.out.print(m[i * n + j]);
      }
      System.out.println();
    }
  }
}
