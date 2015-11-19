package pso;

import java.io.PrintWriter;

/**
 *
 * @author Bruno Andrade
 */
public class Diversity {

	private PSO pso;
	private PrintWriter printer;
	private double diameter;
	private double radius;
	private double avgAroundCenter;
	private double avgOfAvg;
	private double[] center; // To define what is the swarm center

	public Diversity(PSO pso, PrintWriter printer) {
		this.pso = pso;
		this.printer = printer;        
		center = new double[pso.DIMENSION];
	}

	private double euclidianDistance(double[] ds, double[] ds2) {
		double sum = 0.0f;
		for (int i = 0; i < ds2.length; i++) {
			sum += Math.pow(ds[i] - ds2[i], 2);
		}
		sum = Math.sqrt(sum);
		return sum;
	}

	private void calculateDiameter(){

		double[][] particlesDistance = new double[pso.NUMBER_OF_PARTICLES][pso.NUMBER_OF_PARTICLES];
		diameter = -1;

		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			for (int j = i; j < pso.NUMBER_OF_PARTICLES; j++) {
				if(i!=j)		
					particlesDistance[i][j] = euclidianDistance(pso.particle_position[i], pso.particle_position[j]);
					if(particlesDistance[i][j]>diameter){
						diameter = particlesDistance[i][j];
					}
			}
		}

//		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
//			for (int j = i; j < pso.NUMBER_OF_PARTICLES; j++) {
//
//				if(i!=j && particlesDistance[i][j]>diameter){
//					diameter = particlesDistance[i][j];
//				}
//			}
//		}

		printer.println("diameter:#"+ this.pso.current_iteration + " " + diameter);

	}

	private void calculateRadius(){

		double[] particlesDistance = new double[pso.NUMBER_OF_PARTICLES];
		radius = -1;
//		calculateCenterOfSwarm();
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {

			particlesDistance[i] = euclidianDistance(pso.particle_position[i], center);
		}

		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {

			if(particlesDistance[i]>radius){
				radius = particlesDistance[i];
			}
		}

		printer.println("radius:#"+ this.pso.current_iteration + " " + radius);
	}

	private void calculateAvgAroundCenter(){
		double[] particlesDistance = new double[pso.NUMBER_OF_PARTICLES];
//		calculateCenterOfSwarm();
		double sum = 0;
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {

			particlesDistance[i] = euclidianDistance(pso.particle_position[i], center);
		}

		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {

			sum += particlesDistance[i];
		}

		avgAroundCenter = sum/particlesDistance.length; 

		printer.println("average_around_center:#"+ this.pso.current_iteration +  " " + avgAroundCenter);
	}



	private void calculateAvgOfAvg(){
		double[][] particlesDistance = new double[pso.NUMBER_OF_PARTICLES][pso.NUMBER_OF_PARTICLES];
		double sum;
		double[] avg_by_particle = new double[pso.NUMBER_OF_PARTICLES];

		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			for (int j = i; j < pso.NUMBER_OF_PARTICLES; j++) {
				if(i!=j)		
					particlesDistance[i][j] = euclidianDistance(pso.particle_position[i], pso.particle_position[j]);

			}
		}

		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {

			sum=0;
			for (int j = i; j < pso.NUMBER_OF_PARTICLES; j++) {
				if(i!=j)
					sum += particlesDistance[i][j];		
			}

			avg_by_particle[i] = sum/pso.NUMBER_OF_PARTICLES;
		}

		sum=0;

		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			sum += avg_by_particle[i];
		}

		avgOfAvg = sum/pso.NUMBER_OF_PARTICLES;

		printer.println("average_of_average_around_all_particles:#"+ this.pso.current_iteration + " " + avgOfAvg);

	}

	private void calculateCenterOfSwarm(){
		double sum;

		for(int i=0; i<pso.DIMENSION; i++){
			sum=0;
			for(int j=0; j<pso.NUMBER_OF_PARTICLES; j++){
				sum += pso.particle_position[j][i];
			}
			center[i] = sum/pso.NUMBER_OF_PARTICLES;
		}
	}

	public void iterate(){
		calculateCenterOfSwarm();
		calculateDiameter();
		calculateRadius();
		calculateAvgAroundCenter();
		calculateAvgOfAvg();
	}
}



