package pso;
/**
*
* @author Bruno Andrade
*/

import java.io.PrintWriter;
/*
CHECKS:
	NONE 
*/

// "Measuring exploration/exploitation in particle swarms using swarm diversity"
// Evolutionary Computation, 2008. CEC 2008. (IEEE World Congress on Computational Intelligence). IEEE Congress on
// Olorunda, O. ; Engelbrecht A.P.

/*
 (1) The diameter of a swarm is defined as the maximum distance 
     between any 2 particles in the swarm.
 (2) the radius can be calculated as the distance between the 
     swarm center and the particle in the swarm which is furthest 
     away from it.
 */
public class AnalysisOlorunda {

	private PSO pso;
	private PrintWriter printer;
	private double diameter;
	private double radius;
	private double average_distance_around_center;
	private double normalized_average_distance_around_center;
	private double average_distance_among_particles;
	private double[] spatial_swarm_center; 
	private double[] distances_to_center;
	private double[][] distances_among_particles;

	public AnalysisOlorunda(PSO pso, PrintWriter printer) {
		this.pso = pso;
		this.printer = printer;        
		this.spatial_swarm_center = new double[pso.DIMENSION];
		this.distances_to_center = new double[pso.NUMBER_OF_PARTICLES];
		this.distances_among_particles = new double[pso.NUMBER_OF_PARTICLES][pso.NUMBER_OF_PARTICLES];
	}

	private double euclidianDistance(double[] v1, double[] v2) {
		double sum = 0.0d;
		for (int d = 0; d < v2.length; d++) {
			sum += Math.pow(v1[d] - v2[d], 2);
		}
		sum = Math.sqrt(sum);
		return sum;
	}

	private void evaluateDiameter(){
		diameter = -1.0;
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			for (int j = i; j < pso.NUMBER_OF_PARTICLES; j++) {
				if (i != j) {		
					distances_among_particles[i][j] = euclidianDistance(pso.particle_position[i], pso.particle_position[j]);
					if (distances_among_particles[i][j]>diameter) {
						diameter = distances_among_particles[i][j];
					}
				}
			}
		}
		printer.println("diameter:#"+ this.pso.current_iteration + " " + diameter);
	}

	private void evaluateRadius(){
		radius = -1;
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			if(distances_to_center[i]>radius){
				radius = distances_to_center[i];
			}
		}
		printer.println("radius:#"+ this.pso.current_iteration + " " + radius);
	}

	private void evaluateAverageDistanceAroundCenter(){
		double sum = 0;
		for (int particle = 0; particle < pso.NUMBER_OF_PARTICLES; particle++) {
			sum += distances_to_center[particle];
		}
		average_distance_around_center = sum/((float) pso.NUMBER_OF_PARTICLES);
		normalized_average_distance_around_center = average_distance_around_center/diameter; 
		printer.println("average_around_center:#"+ this.pso.current_iteration +  " " + average_distance_around_center);
		printer.println("normalized_average_around_center:#"+ this.pso.current_iteration +  " " + normalized_average_distance_around_center);
	}

	private void evaluateDistancesToCenter() {
		for (int particle = 0; particle < pso.NUMBER_OF_PARTICLES; particle++) {
			distances_to_center[particle] = euclidianDistance(pso.particle_position[particle], spatial_swarm_center);
		}
	}

	private void evaluateAverageOfAverage(){
		double sum;
		double[] average_distance_around_particle = new double[pso.NUMBER_OF_PARTICLES];
		for (int particle = 0; particle < pso.NUMBER_OF_PARTICLES; particle++) {
			sum = 0.0;
			for (int other = 0; other < pso.NUMBER_OF_PARTICLES; other++) {
				if(particle != other) {
					sum += distances_among_particles[particle][other];
				}
			}
			average_distance_around_particle[particle] = sum/((float) pso.NUMBER_OF_PARTICLES);
		}
		sum = 0.0;
		for (int particle = 0; particle < pso.NUMBER_OF_PARTICLES; particle++) {
			sum += average_distance_around_particle[particle];
		}

		average_distance_among_particles = sum/((float) pso.NUMBER_OF_PARTICLES);

		printer.println("average_of_average_around_all_particles:#"+ this.pso.current_iteration + " " + average_distance_among_particles);

	}
	
	private void evaluateSpatialSwarmCenter(){
		double sum;
		for(int d = 0; d < pso.DIMENSION; d++){
			sum = 0.0;
			for(int particle = 0; particle < pso.NUMBER_OF_PARTICLES; particle++){
				sum += pso.particle_position[particle][d];
			}
			spatial_swarm_center[d] = sum/((float) pso.NUMBER_OF_PARTICLES);
		}
	}

	public void iterate(){
		evaluateSpatialSwarmCenter();
		evaluateDistancesToCenter();
		evaluateDiameter();
		evaluateRadius();
		evaluateAverageDistanceAroundCenter();
		evaluateAverageOfAverage();
	}
}



