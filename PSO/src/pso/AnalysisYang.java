package pso;
/**
*
* @author Bruno Andrade
*/

import java.io.PrintWriter;
import java.util.Random;
/*
CHECKS:
	OK (Marcos)
*/

// "A modified particle swarm optimizer with dynamic adaptation" 
// Applied Mathematics and Computation, Volume 189, Issue 2, 15 June 2007, Pages 1205-1213 
// Xueming Yang, Jinsha Yuan, Jiangye Yuan, Huina Mao 

/*
The evolution speed factor and aggregation degree factor: this paper defines the evolution 
speed and aggregation degree using the following formula:
	(1) Evolutionary speed factor
	    
	    h_i^t = |[min(F(pbest_i^(t-1)), F(pbest_i^t)) / max(F(pbest_i^(t-1)), F(pbest_i^t))]|
	
	where F(pbest_i^t) is the fitness value of pbest_i^t. Under the assumption and definition 
	above, it can be obtained that 0 < h <= 1. This parameter takes the run history of each 
	particle into account, and reflects the evolutionary speed of each particle, that is, 
	the smaller the value of h, the faster the speed.

	(2) Aggregation degree
		
		s = | min(F_tbest, Fbar_t) / max(F_tbest, Fbar_t) | 
		 
	where Fbar_t is the mean fitness of all particles in the swarm at the tth iteration. 
	Note that F(gbest_t) cannot be substituted for F_tbest, since F_tbest represents the optimal value 
	found in this iteration, while F(gbest_t) the optimal value that the whole swarm has found 
	up to the tth iteration. Compared with that in [8], the formula of aggregation degree factor 
	in this paper can faster response to evolving state.
 */
public class AnalysisYang {
	
	private PSO pso;
    private PrintWriter printer;
	private double[] previousPBest;
	private double[] currentPBest;
	private double fTBest;
	private double averageFitness;
	private double[] evoFactor;
	private double aggregationFactor;
	private double[] inertiaWeights;
	
	public AnalysisYang(PSO pso, PrintWriter printer) {
		this.pso = pso;
        this.printer = printer;
		previousPBest = pso.particle_best_value; 
		currentPBest = pso.particle_best_value;
		evoFactor = new double[pso.NUMBER_OF_PARTICLES];
		inertiaWeights = new double[pso.NUMBER_OF_PARTICLES];
		for (int i = 0; i < inertiaWeights.length; i++) {
			inertiaWeights[i] = 1;
		}
	}

	public PSO getPso() {
		return pso;
	}

	public double[] getPreviousPBest() {
		return previousPBest;
	}

	public double[] getCurrentPBest() {
		return currentPBest;
	}

	public double getfTBest() {
		return fTBest;
	}

	public double getAverageFitnessOfAll() {
		return averageFitness;
	}

	public double[] getEvoFactor() {
		return evoFactor;
	}

	public double getAggregFactor() {
		return aggregationFactor;
	}

	public double[] getInertiaWeights() {
		return inertiaWeights;
	}

	private void calculateEvolutonarySpeedFactor(){
		// we do not use this and will not. 
		for (int i=0; i < pso.NUMBER_OF_PARTICLES; i++) {
			evoFactor[i] = Math.min(previousPBest[i], currentPBest[i])/ Math.max(previousPBest[i], currentPBest[i]);
			evoFactor[i] = Math.abs(evoFactor[i]);
		}
	}
	
	private void calculateAggregationDegree(){
		// s = | min(F_tbest, Fbar_t) / max(F_tbest, Fbar_t) | 
		aggregationFactor = Math.min(fTBest, averageFitness)/ Math.max(fTBest, averageFitness);
		aggregationFactor = Math.abs(aggregationFactor);
        printer.println("Aggregation Factor:#"+ this.pso.current_iteration + " " + aggregationFactor);
	}
	
	private void updatePBest(){
		previousPBest = currentPBest;
		currentPBest = pso.particle_best_value;
	}
	
	private void updateAverageFitness(){
		// Fbar_t is the mean fitness of all particles in the swarm at the tth iteration
		double sum = 0.0;
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			sum += pso.particle_current_fitness[i];
		}
		averageFitness = sum/pso.NUMBER_OF_PARTICLES;
	}
	
	private void updateFTBest(){
		// F_tbest represents the optimal value found in THIS iteration
		fTBest = Double.MAX_VALUE;
		for (int i=0; i<pso.NUMBER_OF_PARTICLES; i++) {
			if (pso.particle_current_fitness[i] < fTBest) {
				fTBest = pso.particle_current_fitness[i];
			}
		}
	}
	
	public void iterate(){
		updatePBest();
		updateAverageFitness();
		updateFTBest();
		calculateEvolutonarySpeedFactor();
		calculateAggregationDegree();
	}
}
