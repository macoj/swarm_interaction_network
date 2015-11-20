package pso;

// "Adaptive Particle Swarm Optimization"
// Zhi-Hui Zhan, Jun Zhang, Yun Li, and Henry Shu-Hung Chung
/*
Step 1: At the current position, calculate the mean distance
of each particle i to all the other particles. For
example, this mean distance can be measured using
an Euclidian metric

d_i = 1/(N-1) x sum_{j=1, i!=j}^N distance(i, j)

where N and D are the population size and the
number of dimensions, respectively.

Step 2: Denote d_i of the globally best particle as d_g. Compare all 
di's, and determine the maximum and minimum distances dmax and dmin. 
Compute an "evolutionary factor" f as defined by
		f = (d_g - d_min)/(d_max - d_min)
*/
public class AnalysisZhan {

}
