package pso;
/**
*
* @author Bruno Andrade
*/

import java.io.PrintWriter;
import java.util.Random;

// "A modified particle swarm optimizer with dynamic adaptation" 
// Applied Mathematics and Computation, Volume 189, Issue 2, 15 June 2007, Pages 1205-1213 
// Xueming Yang, Jinsha Yuan, Jiangye Yuan, Huina Mao 

/*
 NOT CHECKED 
*/
public class AnalysisYang {
	
	private PSO pso;
    private PrintWriter printer;
	private double[] previousPBest;
	private double[] currentPBest;
	private double fTBest;
	private double averageFitnessOfAll;
	private double[] evoFactor;
	private double aggregFactor;
	private double[] inertiaWeights;
	
	public AnalysisYang(PSO pso, PrintWriter printer) {
		this.pso = pso;
        this.printer = printer;
		previousPBest = pso.particle_best_value;
		currentPBest = pso.particle_best_value;
		evoFactor = new double[pso.NUMBER_OF_PARTICLES];
		inertiaWeights = new double[pso.NUMBER_OF_PARTICLES];
		
		for (double w : inertiaWeights) {
			w = 1;
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
		return averageFitnessOfAll;
	}

	public double[] getEvoFactor() {
		return evoFactor;
	}

	public double getAggregFactor() {
		return aggregFactor;
	}

	public double[] getInertiaWeights() {
		return inertiaWeights;
	}

	private void calculateEvoFactor(){
				
		for(int i=0; i<pso.NUMBER_OF_PARTICLES; i++){
			evoFactor[i] = Math.min(previousPBest[i], currentPBest[i])/ Math.max(previousPBest[i], currentPBest[i]);
		}
		
	}
	
	private void calculateAggregFactor(){
		aggregFactor = Math.min(fTBest, averageFitnessOfAll)/ Math.max(fTBest, averageFitnessOfAll);
        printer.println("Aggregation Factor:#"+ this.pso.current_iteration + " " + aggregFactor);
	}
	
	private void calculateInertiaWeight(){
		Random r = new Random();
		double a = r.nextDouble();
		double b = r.nextDouble();
		
		for(int i=0; i<pso.NUMBER_OF_PARTICLES; i++){
			inertiaWeights[i] = inertiaWeights[i] - a * (1 - evoFactor[i]) + b*aggregFactor;
		}
		
	}
			
	private void updatePBest(){
		previousPBest = currentPBest;
		currentPBest = pso.particle_best_value;
	}
	
	private void updateAverageFitnessofAllAndFTBest(){
		
		double sum = 0;
		double[] current_fitness = new double[pso.NUMBER_OF_PARTICLES];
		fTBest = Double.MAX_VALUE;
		
		for(int i=0; i<pso.NUMBER_OF_PARTICLES; i++){
			
			current_fitness[i] = pso.FUNCTION.compute(pso.particle_position[i]);
			sum += current_fitness[i];
		}
		
		for(int i=0; i<pso.NUMBER_OF_PARTICLES; i++){
			if(current_fitness[i]<fTBest)
				fTBest = current_fitness[i];
		}
		
		averageFitnessOfAll = sum/pso.NUMBER_OF_PARTICLES;
	}
	
	public void iterate(){
		
		updatePBest();
		updateAverageFitnessofAllAndFTBest();
		calculateEvoFactor();
		calculateAggregFactor();
		calculateInertiaWeight();
	}
}
