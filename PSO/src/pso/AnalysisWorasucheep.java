package pso;
/**
*
* @author Bruno Andrade
*/

import java.io.PrintWriter;
import java.util.LinkedList;
/*
 CHECKS:
	OK (Marcos)
*/


// "A Particle Swarm Optimization with stagnation detection and dispersion"
// Evolutionary Computation, 2008. CEC 2008. (IEEE World Congress on Computational Intelligence). IEEE Congress on
// Worasucheep, C.

/*
To measure the stagnation in PSO, [5] takes into consideration also the velocities
and defines the improvement ratio as:
	R = (1 - fc/fp)/(1 - vc/vp),
where fc is the current fitness value of the best particle, fp is the previous and vc
is the current average velocity of all particles, while vp is the previously recorded
one. Stagnation is detected when R drops under a preset value 'e'.
 */
public class AnalysisWorasucheep {
	
	private PSO pso;
	private PrintWriter printer;
	private LinkedList<Double> window = new LinkedList<Double>();
	private double previous_best_value;
	private double current_best_value;
	private double current_average_particles_velocity;
	private double previous_average_particles_velocity;
	private double ratio;
	private final int WINDOW_SIZE = 500;
	private final int CONVERGENCE_CHECK = 1;
	
    public AnalysisWorasucheep(PSO pso, PrintWriter printer) {
		this.pso = pso;
		this.printer = printer;
		ratio = Double.NaN;
	}
        
	public double getAvg_particles_velocity() {
		return current_average_particles_velocity;
	}
	public double getPrevious_avg_particles_velocity() {
		return previous_average_particles_velocity;
	}

	public double getR() {
		return ratio;
	}

	public double getCurrentFbest() {
		return current_best_value;
	}
	
	public double getPreviousFbest(){
		return this.previous_best_value;
	}

	public void evaluateRatio(){
		if(pso.current_iteration >= WINDOW_SIZE){
			double fc = current_best_value;
			double fp = previous_best_value;
			double vc = current_average_particles_velocity;
			double vp = getAverageVelocityFromWindow();
			ratio = (1 - fc/fp)/(1 - vc/vp);
			ratio = Math.abs(ratio);
			String rLine = "R:#" + this.pso.current_iteration + " " + ratio;
			printer.println(rLine + "\n");
		}
	}

	private double evaluateParticlesAverageVelocity(){
		double sum;
		double[] norms_of_particles_velocities = new double[this.pso.DIMENSION];
		
		previous_average_particles_velocity = current_average_particles_velocity;
		
		for (int particle = 0; particle < this.pso.NUMBER_OF_PARTICLES; particle++) {   
            sum = 0d;
            for(int dimension = 0; dimension < this.pso.DIMENSION; dimension++) {
            	sum += Math.pow(this.pso.particle_velocity[particle][dimension], 2);
            }
            norms_of_particles_velocities[particle] = Math.sqrt(sum);
		}
		
		sum = 0d;
		
		for (int particle = 0; particle < this.pso.NUMBER_OF_PARTICLES; particle++) {
			sum += norms_of_particles_velocities[particle];
		}
		
		current_average_particles_velocity = sum/((double) this.pso.NUMBER_OF_PARTICLES);
		
		if(pso.current_iteration < WINDOW_SIZE) {
			window.add(current_average_particles_velocity);
		} else {
			window.poll();  // is the running time constant?  
			window.add(current_average_particles_velocity);
		}
		
		if (this.pso.current_iteration == 0) {
			previous_average_particles_velocity = Double.NaN;
		}
			
		return current_average_particles_velocity;
	}	
	
	private void updateBestSwarmFitness(){
		previous_best_value = current_best_value;
		current_best_value = pso.particle_best_value[pso.swarm_gbest];
		if (pso.current_iteration == 0) {
			previous_best_value = Double.NaN;
		}
	}
	
	private double getAverageVelocityFromWindow(){
		double sum = 0d;
		double averageVelocityByWindow;
		
		for (int i = 0; i < this.window.size(); i++) {
			sum += window.get(i);
		}
		
		averageVelocityByWindow = sum/((double) this.WINDOW_SIZE);
		
		return averageVelocityByWindow;
	}
	
	public void iterate() {
		this.evaluateParticlesAverageVelocity();
		this.updateBestSwarmFitness();
		if (pso.current_iteration%CONVERGENCE_CHECK == 0) {
			this.evaluateRatio();
		}
	}
}