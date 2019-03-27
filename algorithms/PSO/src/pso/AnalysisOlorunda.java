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
public class AnalysisOlorunda implements Analysis {

	private PSO pso;
	private PrintWriter printer;
	private double diameter;
	private double radius;
	private double average_distance_around_center;
	private double normalized_average_distance_around_center;
	private double average_distance_among_particles;
	private double[] position_swarm_center; 
	private double[] velocity_swarm_center;
	private double[] distances_to_center;
	private double[][] distances_among_particles;

	public AnalysisOlorunda(PSO pso, PrintWriter printer) {
		this.pso = pso;
		this.printer = printer;        
		this.position_swarm_center = new double[pso.DIMENSION];
		this.velocity_swarm_center = new double[pso.DIMENSION];
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
			distances_to_center[particle] = euclidianDistance(pso.particle_position[particle], position_swarm_center);
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
	
	private void evaluateCoherence() {
		double coherence;
		double swarm_center_speed = 0d;
		double average_particle_speed = 0d;
		
		// speed of the center of the swarm:
		for (int d = 0; d < pso.DIMENSION; d++) {
			swarm_center_speed += Math.pow(velocity_swarm_center[d], 2);  // apparently we should multiply by NUMBER_OF_PARTICLES
		}
		swarm_center_speed = Math.sqrt(swarm_center_speed) / ((float) pso.NUMBER_OF_PARTICLES);
		
		// average speed of the particles
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			double speed_i = 0d;
			for (int d = 0; d < pso.DIMENSION; d++) {
				speed_i += Math.pow(pso.particle_velocity[i][d], 2);  
			}
			speed_i = Math.sqrt(speed_i) ;
			average_particle_speed += speed_i; 
		}
		average_particle_speed = average_particle_speed / ((float) pso.NUMBER_OF_PARTICLES); 
		coherence = swarm_center_speed / average_particle_speed; 
		// if particles are moving, but the center is not, coherence tends to zero. 
		printer.println("coherence:#"+ this.pso.current_iteration + " " + coherence);
	}
	
	private void evaluateSwarmCenter(){
		double position_dimension;
		double velocity_dimension;
		for(int d = 0; d < pso.DIMENSION; d++){
			position_dimension = 0.0;
			velocity_dimension = 0.0;
			for(int particle = 0; particle < pso.NUMBER_OF_PARTICLES; particle++){
				position_dimension += pso.particle_position[particle][d];
				velocity_dimension += pso.particle_velocity[particle][d];
			}
			position_swarm_center[d] = position_dimension/((float) pso.NUMBER_OF_PARTICLES);
			velocity_swarm_center[d] = velocity_dimension/((float) pso.NUMBER_OF_PARTICLES);
		}
	}

	public void iterate(){
		evaluateSwarmCenter();
		evaluateDistancesToCenter();
		evaluateDiameter();
		evaluateRadius();
		evaluateAverageDistanceAroundCenter();
		evaluateAverageOfAverage();
		evaluateCoherence();
	}
}



