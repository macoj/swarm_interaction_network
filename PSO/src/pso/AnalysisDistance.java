package pso;

import java.io.PrintWriter;

public class AnalysisDistance implements Analysis {
	
	private PSO pso;
	private PrintWriter printer;
	
	public AnalysisDistance(PSO pso, PrintWriter printer) {
		this.pso = pso;
		this.printer = printer;        
	}

	private double euclidianDistance(double[] v1, double[] v2) {
		double sum = 0.0d;
		for (int d = 0; d < v2.length; d++) {
			sum += Math.pow(v1[d] - v2[d], 2);
		}
		sum = Math.sqrt(sum);
		return sum;
	}
	
	public void iterate() {
		//TODO: use the concept of vectors to speed it up.
		double[][] distances = new double[pso.NUMBER_OF_PARTICLES][pso.NUMBER_OF_PARTICLES];

		this.printer.print("distance:# ");
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			for (int j = i; j < pso.NUMBER_OF_PARTICLES; j++) {
				if (i != j) {		
					distances[i][j] = euclidianDistance(pso.particle_position[i], pso.particle_position[j]);
					distances[j][i] = distances[i][j];
				}
			}
		}
		for (int i = 0; i < pso.NUMBER_OF_PARTICLES; i++) {
			for (int j = 0; j < pso.NUMBER_OF_PARTICLES; j++) {
				this.printer.print(distances[i][j] + " ");
			}
		}
		this.printer.print("\n");
	}

}

