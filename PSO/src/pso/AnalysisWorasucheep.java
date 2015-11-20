package pso;
/**
*
* @author Bruno Andrade
*/

import java.io.PrintWriter;
import java.util.LinkedList;

// "A Particle Swarm Optimization with stagnation detection and dispersion"
// Evolutionary Computation, 2008. CEC 2008. (IEEE World Congress on Computational Intelligence). IEEE Congress on
// Worasucheep, C.

/*
 NOT CHECKED 
*/
public class AnalysisWorasucheep {
	
	private PSO pso;
	private PrintWriter printer;
	private LinkedList<Double> window = new LinkedList<Double>();
	private double previousFbest;
	private double currentFbest;
	private double avg_particles_velocity;
	private double previous_avg_particles_velocity;
	private double ratio;
	private final int WINDOW_SIZE = 500;
	private final int CONVERGENCE_CHECK = 1;
	
        public AnalysisWorasucheep(PSO pso, PrintWriter printer){
		this.pso = pso;
		this.printer = printer;
		ratio = Double.NaN;
	}
        
	public double getAvg_particles_velocity() {
		return avg_particles_velocity;
	}
	public double getPrevious_avg_particles_velocity() {
		return previous_avg_particles_velocity;
	}

	public double getR() {
		return ratio;
	}

	public double getCurrentFbest() {
		return currentFbest;
	}
	
	public double getPreviousFbest(){
		return this.previousFbest;
	}

	public void detectStagnation(){
		
		if(pso.current_iteration>=WINDOW_SIZE){
			double fc = currentFbest;
			double fp = previousFbest;
			double vc = avg_particles_velocity;
			double vp = getAverageVelocityByWindow();
			
			
			ratio = Math.abs(1-(fc/fp)/1-(vc/vp));
			if(ratio<Math.pow(10, -5)){
				//runDispersion();
			}
		}
		
	}

	private void runDispersion() {
		// In order to perform this action, it would be necessary to change the "updateVelocity" method in the PSO class,
		// making it to receive parameters.
	}
	
	private double getAvgVelocityOfAllParticles(){
		
		double sum;
		double[] avg_by_particle = new double[this.pso.DIMENSION];
		
		
		// average velocity of all particles
		for (int i=0; i<this.pso.particle_velocity.length; i++){   
			
//			sum = 0;
//			for(int j=0; j<this.pso.DIMENSION; j++){
//				sum += this.pso.particle_velocity[i][j];
//			}
//			
//			avg_by_particle[i] = sum/this.pso.DIMENSION;
                    
                    sum = 0;
                    for(int j=0; j<this.pso.DIMENSION; j++){
			sum += Math.pow(this.pso.particle_velocity[i][j], 2);
                    }
			
                    avg_by_particle[i] = Math.sqrt(sum);
			
		}
		
		sum=0;
		for(int i=0; i<this.pso.particle_velocity.length; i++){
			sum += avg_by_particle[i];
		}
		
		previous_avg_particles_velocity = avg_particles_velocity;
		
		avg_particles_velocity = sum/this.pso.particle_velocity.length;
		
		if(pso.current_iteration<WINDOW_SIZE)
			window.add(avg_particles_velocity);
		else{
			window.poll();
			window.add(avg_particles_velocity);
		}
		
		if(this.pso.current_iteration==0)
			previous_avg_particles_velocity = Double.NaN;
		
		return avg_particles_velocity;
	}	
	
	private void updateBestSwarmFitness(){
		previousFbest = currentFbest;
		currentFbest = pso.particle_best_value[pso.swarm_gbest];
		if(pso.current_iteration==0)
			previousFbest = Double.NaN;
	}
	
	private void updatePreviousBestSwarmFitness(){
		this.previousFbest = this.currentFbest;
	}
	
	private double getAverageVelocityByWindow(){
		double sum = 0;
		double averageVelocityByWindow;
		
		for(int i=0; i<window.size(); i++){
			sum += window.get(i);
		}
		
		averageVelocityByWindow = sum/window.size();
		
		return averageVelocityByWindow;
	}
	
	public void iterate(){
		this.getAvgVelocityOfAllParticles();
		updateBestSwarmFitness();
		String rLine = "R:#" + this.pso.current_iteration + " " + ratio;
		printer.println(rLine + "\n");
		if(pso.current_iteration%CONVERGENCE_CHECK==0)
			this.detectStagnation();
				
	}
}