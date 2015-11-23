package pso;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;
import benchmark_functions.*;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;

public class PSO implements Runnable {
	
    private PrintWriter printWriter;
	
	double swarm_evolutionary_state_1;
	double swarm_evolutionary_state_2;
	double swarm_evolutionary_state_2_reduced;
	
	public int swarm_gbest;
	
	int current_iteration;
	
	boolean swarm_neighborhood_graph[][];
	int swarm_influence_graph[][] = null;
	int swarm_influence_history[][] = null;
	int swarm_influence_history_maximum_length;
	boolean swarm_influence_graph_weighted;
	
	double swarm_neighborhood_average_path_length;

	Graph<Integer, String> swarm_neigbourhood_graph_jung;
	
	TOPOLOGY rr;
	TOPOLOGY swarm_initial_topology;
	TOPOLOGY_MECHANISM swarm_topology_mechanism;
	
	double swarm_random_topology_p; 
	
	double particle_position[][];
	double particle_velocity[][];
	double particle_current_fitness[];
	double particle_best_value_position[][];
	double particle_best_value[];
	int particle_best_neighbour[];
	int particle_failures[];
	double particle_failures_power[];
	
	double swarm_distance_particles[][];
	double[][] swarm_centroids = null;
	double[][] swarm_clusters_kohonen = null;
	int[] swarm_clusters = null;
	
	
	double PARTICLE_MAXX, PARTICLE_MAXV, PARTICLE_MINX, PARTICLE_INITIAL_RANGE_L, PARTICLE_INITIAL_RANGE_R;
	//it: 9999 gbest: 40 : 19.71921474112517
	int DIMENSION, NUMBER_OF_PARTICLES, MAXITER, RUNS;
	int swarm_number_of_clusters;
	Random random_number_generator;
	long random_number_generator_seed;
	
	Random random_number_generator_independent;
	long random_number_generator_seed_independent;
	
	
	int swarm_initial_maximum_neighbors;
	Function FUNCTION;
	
	enum TOPOLOGY { GLOBAL, RING, CLAN, RANDOM, VON_NEUMANN, BARABASI, PARTNERS, THREESOME_PARTNERS, ORGY_PARTNERS, NSOME_PARTNERS, K_MEANS_BASED, NO_TOPOLOGY, KOHONEN_BASED };
	
	enum TOPOLOGY_MECHANISM { NONE, BARABASI_BASED, BARABASI_BASED_DYNAMIC_A, BARABASI_BASED_DYNAMIC_B, BARABASI_BASED_DYNAMIC_C, BARABASI_BASED_DYNAMIC_D, DYNAMIC_2011, BARABASI_BASED_DYNAMIC_LOG };
	
	static final int swarm_topology_sleep = 1;
	
	int particles_failures_threshold;
	int[] particles_failures_threshold_particular = null;

/*
	F1  : Shifted Elliptic Function 									UNIMODAL
	F2  : Shifted Rastrigin's Function									MULTIMODAL
	F3  : Shifted Ackley's Function										MULTIMODAL
	F4  : Single-group Shifted and m-rotated Elliptic Function 			UNIMODAL
	F5  : Single-group Shifted and m-rotated Rastrigin's Function		MULTIMODAL
	F6  : Single-group Shifted and m-rotated Ackley's Function			MULTIMODAL		FUCKIN' BETTER! (a inicializacao deve ter muitas cores! -- Mulligan?)
	F7  : Single-group Shifted m-dimensional Schwefel's Problem 1.2 	UNIMODAL		WORSE
	F8  : Single-group Shifted m-dimensional Rosenbrock's Function		MULTIMODAL
	F9  : D/2m-group Shifted and m-rotated Elliptic Function 			UNIMODAL
	F10 : D/2m-group Shifted and m-rotated Rastrigin's Function			MULTIMODAL
	F11 : D/2m-group Shifted and m-rotated Ackley's Function			MULTIMODAL		SLIGHTLY BETTER
	F12 : D/2m-group Shifted m-dimensional Schwefel's Problem 1.2 		UNIMODAL
	F13 : D/2m-group Shifted m-dimensional Rosenbrock's Function		MULTIMODAL
	F14 : D/m-group Shifted and m-rotated Elliptic Function 			UNIMODAL
	F15 : D/m-group Shifted and m-rotated Rastrigin's Function			MULTIMODAL
	F16 : D/m-group Shifted and m-rotated Ackley's Function				MULTIMODAL		SLIGHTLY BETTER
	F17 : D/m-group Shifted m-dimensional Schwefel's Problem 1.2 		UNIMODAL
	F18 : D/m-group Shifted m-dimensional Rosenbrock's Function			MULTIMODAL
	F19 : D/m-group Shifted Schwefel's Problem 1.2 						UNIMODAL
	F20 : D/m-group Shifted Rosenbrock's Function						MULTIMODAL
*/
	public PSO(String[] args) {
		
		int evaluation = 1000000;
		// args[0] is the function
		if (args.length > 0) {
			FUNCTION = functionFromInt(Integer.parseInt(args[0]));
			
		} else {
			FUNCTION = new F6();	// it uses F6 otherwise
		}
		DIMENSION = 1000;
		NUMBER_OF_PARTICLES = 100; 
		MAXITER = evaluation / NUMBER_OF_PARTICLES;
		//MAXITER=50;
		RUNS = 1;
		PARTICLE_MAXX = FUNCTION.getMax();
		PARTICLE_MINX = FUNCTION.getMin();
		
		PARTICLE_MAXV = 0.5f*(PARTICLE_MAXX - PARTICLE_MINX);
		   
		PARTICLE_INITIAL_RANGE_L = PARTICLE_MINX;
		PARTICLE_INITIAL_RANGE_R = PARTICLE_MAXX;
		
		particles_failures_threshold = 50;
		current_iteration = 0;
		
		swarm_neighborhood_graph = new boolean[NUMBER_OF_PARTICLES][NUMBER_OF_PARTICLES];

		particle_position = new double[NUMBER_OF_PARTICLES][DIMENSION];
		particle_velocity = new double[NUMBER_OF_PARTICLES][DIMENSION];
		particle_best_value_position = new double[NUMBER_OF_PARTICLES][DIMENSION];
		particle_best_value = new double[NUMBER_OF_PARTICLES];
		particle_current_fitness = new double[NUMBER_OF_PARTICLES];
		particle_best_neighbour = new int[NUMBER_OF_PARTICLES];
		particle_failures = new int[NUMBER_OF_PARTICLES];
		
		particle_failures_power = new double[NUMBER_OF_PARTICLES];
		particles_failures_threshold_particular = new int[NUMBER_OF_PARTICLES];
		initializeRNG();
		
		swarm_influence_graph_weighted = false; // TODO consertar isso!
		
		swarm_influence_history_maximum_length = 5;
		swarm_influence_graph = new int[NUMBER_OF_PARTICLES][NUMBER_OF_PARTICLES];
		swarm_influence_history = new int[NUMBER_OF_PARTICLES][swarm_influence_history_maximum_length];
		for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
			for (int j = 0; j < NUMBER_OF_PARTICLES; j++) {
				swarm_influence_graph[i][j] = 0;
			}
			for (int j = 0; j < swarm_influence_history_maximum_length; j++) {
				swarm_influence_history[i][j] = 0;
			}
		}

		//swarm_initial_topology = TOPOLOGY.K_MEANS_BASED;
		//swarm_initial_topology = TOPOLOGY.KOHONEN_BASED;
		//swarm_topology_mechanism = TOPOLOGY_MECHANISM.BARABASI_BASED_DYNAMIC_D;
		//swarm_topology_mechanism = TOPOLOGY_MECHANISM.BARABASI_BASED_DYNAMIC_A;
		//swarm_topology_mechanism = TOPOLOGY_MECHANISM.BARABASI_BASED_DYNAMIC_C;
		//swarm_topology_mechanism = TOPOLOGY_MECHANISM.BARABASI_BASED_DYNAMIC_DISSERTACAO;
		//swarm_topology_mechanism = TOPOLOGY_MECHANISM.BARABASI_BASED_DYNAMIC_LOG;
		//it: 1499 gbest: 72 : 5.866975361772289E8
		
		//swarm_topology_mechanism = TOPOLOGY_MECHANISM.NONE;
		//swarm_initial_topology = TOPOLOGY.RING;
		//swarm_initial_topology = TOPOLOGY.RANDOM;
		//swarm_random_topology_p = 0.2;
		
		//swarm_initial_topology = TOPOLOGY.VON_NEUMANN;
		swarm_initial_topology = TOPOLOGY.GLOBAL;
		//swarm_initial_topology = TOPOLOGY.CLAN;
		//swarm_initial_topology = TOPOLOGY.NO_TOPOLOGY;
		swarm_gbest = -1;
		
		swarm_number_of_clusters = 3;
		
		swarm_initial_maximum_neighbors = 3;
		
		swarm_distance_particles = new double[NUMBER_OF_PARTICLES][NUMBER_OF_PARTICLES];
		this.swarm_neigbourhood_graph_jung = new SparseMultigraph<Integer, String>();
		for (int vertex = 0; vertex < NUMBER_OF_PARTICLES; vertex++) {
			this.swarm_neigbourhood_graph_jung.addVertex((Integer) vertex);
		}
		
	}
	
	private Function functionFromInt(int parseInt) {
		Function func_return = null;
		switch (parseInt) {
		case 1:
			func_return = new F1();
			break;
		case 2:
			func_return = new F2();
			break;
		case 3:
			func_return = new F3();
			break;
		case 4:
			func_return = new F4();
			break;
		case 5:
			func_return = new F5();
			break;
		case 6:
			func_return = new F6();
			break;
		case 7:
			func_return = new F7();
			break;
		case 8:
			func_return = new F8();
			break;
		case 9:
			func_return = new F9();
			break;
		case 10:
			func_return = new F10();
			break;
		case 11:
			func_return = new F11();
			break;
		case 12:
			func_return = new F12();
			break;
		case 13:
			func_return = new F13();
			break;
		case 14:
			func_return = new F14();
			break;
		case 15:
			func_return = new F15();
			break;
		case 16:
			func_return = new F16();
			break;
		case 17:
			func_return = new F17();
			break;
		case 18:
			func_return = new F18();
			break;
		case 19:
			func_return = new F19();
			break;
		case 20:
		default:
			func_return = new F20();
			break;				
		}
		return func_return;
	}
	
	private void initializeRNG() {
		//it:#0 1500 2013336.6304360947
		//random_number_generator_seed: 1339596216600
		//it:#0 1500 1482251.7006799742
		//random_number_generator_seed: 1339596216600
		//it:#0 1500 1482251.7006799742
		//random_number_generator_seed: 1339596216600
		random_number_generator_seed = (long) 1353096711862L;
		//1353096711862
		//it:#0 1278 2358.9153
		random_number_generator_seed = System.currentTimeMillis();
		random_number_generator = new Random(random_number_generator_seed);
		random_number_generator_seed_independent = System.currentTimeMillis();
		//random_number_generator_seed_independent = (long) 1353375901830L; // random_number_generator_seed;
		random_number_generator_independent = new Random(random_number_generator_seed_independent);
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		PSO pso = new PSO(args);
		pso.run();
	}
        
    public void setTopology(TOPOLOGY topology){
        this.swarm_initial_topology = topology;
    }
    
    public void setTopologyMechanism(TOPOLOGY_MECHANISM topologyMechanism){
        this.swarm_topology_mechanism = topologyMechanism;
    }

    public void setPrinter(PrintWriter printer){
        this.printWriter = printer;
    }
    
    /**
     *
     * @throws FileNotFoundException
     */
    @Override
	public void run(){
    	DecimalFormat decimal_format = new DecimalFormat("#.####");
        AnalysisOlorunda diversity = new AnalysisOlorunda(this, printWriter);
        AnalysisYang dapso = new AnalysisYang(this, printWriter);
        AnalysisWorasucheep psoDD = new AnalysisWorasucheep(this, printWriter);
        
        double run_final_values[] = new double[RUNS];
        for (int run = 0; run < RUNS; run++) {
            this.initializeRNG();
            this.initializePSO();

            do {
                this.evaluatePositionAndUpdatePersonal();
                this.updateGBest();
                this.updateNBest();
                this.updateNeighbourhood();
                this.updateVelocity();
                this.updatePosition();
                this.updateInfluenceGraph();
                
                String it = "it:#" + this.current_iteration + " " + particle_best_value[this.swarm_gbest];
                printWriter.println(it);

                String igLine = "ig:#" + this.current_iteration + " ";
                printWriter.print(igLine);
                int[][] ig = this.getInfluenceGraph();
                for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
                    for (int j = 0; j < NUMBER_OF_PARTICLES; j++) {
//						System.out.print(ig[i][j]+" ");
                        printWriter.print(ig[i][j]+" ");
                    }                            
                }
                printWriter.println();
                
                printWriter.print("pbest:#" + this.current_iteration + " ");
                for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
                    printWriter.print(this.particle_best_value[i] + " ");
                }
                printWriter.println();
                
                //diversity.iterate();
                //dapso.iterate();
                //psoDD.iterate();
                printWriter.println();
                
                this.current_iteration = this.current_iteration + 1;
                /*
                try {
                    Thread.sleep(1);
                } catch (InterruptedException ex) {
                    Logger.getLogger(PSO.class.getName()).log(Level.SEVERE, null, ex);
                }
                */
            } while (this.current_iteration < MAXITER);
            run_final_values[run] = particle_best_value[this.swarm_gbest];
        }
        printWriter.close();
	}

	private void updateInfluenceGraph() {
		// updates the history overall
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			// the particle received information from these particles in the last iterations...
			swarm_influence_history[particle][this.current_iteration%swarm_influence_history_maximum_length] = particle_best_neighbour[particle];
		}
		
		for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
			for (int j = 0; j < NUMBER_OF_PARTICLES; j++) {
				swarm_influence_graph[i][j] = 0;
			}
		}
		
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			swarm_influence_graph[particle_best_neighbour[particle]][particle]++;// [i][j] = i sends information to j
		}
		
	}
	
	private void updateNeighbourhood() {
		if (swarm_topology_mechanism != null) {
			switch (swarm_topology_mechanism) {
			case DYNAMIC_2011:
				this.updateNeighbourhood_BA_DYNAMIC_DISSERTACAO();
				break;
			default:
				break;
			}
		}
		
	}
	
	private void updateNeighbourhood_BA_DYNAMIC_DISSERTACAO() {
		// we do not want to update the neighborhood every time following
		// the same sequence of particles, so we randomize the index and
		// update following this random sequence.
		Random random_number_generator = null;
		random_number_generator = this.random_number_generator_independent;
		//random_number_generator = this.random_number_generator;   
		int ranking[] = new int[NUMBER_OF_PARTICLES];
		int ranking_of_particle[] = new int[NUMBER_OF_PARTICLES];
		int particle_index_shuffled[] = new int[NUMBER_OF_PARTICLES];
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			ranking[particle] = particle;
			particle_index_shuffled[particle] = particle;
		}
		
		// shuffle the index shuffled:
		for (int particle = 0; particle < 2*NUMBER_OF_PARTICLES; particle++) {
			int particle_a = Math.abs(random_number_generator.nextInt()%NUMBER_OF_PARTICLES);
			int particle_b = Math.abs(random_number_generator.nextInt()%NUMBER_OF_PARTICLES);
			int aux = particle_index_shuffled[particle_a];
			particle_index_shuffled[particle_a] = particle_index_shuffled[particle_b];
			particle_index_shuffled[particle_b] = aux; 
		}
		
		// a dirty bubble sort
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
				if (particle_best_value[ranking[particle]] > particle_best_value[ranking[neighbour]]) {
					int swap = ranking[neighbour];
					ranking[neighbour] = ranking[particle];
					ranking[particle] = swap;
				}
			}
		}
	
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			ranking_of_particle[ranking[particle]] = NUMBER_OF_PARTICLES - particle - 1;
		}
		// ranking[0] is the index of the worst particle
		// ranking[NUMBER_OF_AGENTS - 1] is the index of the best particle
		// ranking_of_particle[i] is the position of the particle i in the rank; if ranking_of_particle[i] = 0, it means that i is the worst particle 
		
		//random_number_generator = this.random_number_generator_independent;
		// in the same order of particle updating, the swarm behaviour is different
		// when different seeds are used. 
		
		//// -> esse laco pega uma particula aleatoriamente no enxame!
		for (int particle_ordered = 0; particle_ordered < ranking_of_particle.length; particle_ordered++) {
			
			int particle = particle_index_shuffled[particle_ordered];
			
			if (particle_failures[particle] < particles_failures_threshold_particular[particle]) {
				continue;
			}
	
			double r;  
			int who = 0;
			float sum = 0;
			double power = 1; 
						
			float pa_sum = (float) pa(NUMBER_OF_PARTICLES, power);
			r = random_number_generator.nextDouble();
			sum = 0;
			who = NUMBER_OF_PARTICLES - 1; 	// who = -1; se passar da roulette wheel, significa que o power esta muito grande
											//			 e a precisao nao eh suficiente, entao a divisao com pa_sum nao eh valida :/  
			int wheel;
			for (wheel = 0; wheel < NUMBER_OF_PARTICLES ; wheel++) {
				// ESCOLHE PELA ROULETTE WHEEL UM PARA SE CONECTAR
				sum += Math.pow(((double) wheel + 1), power) / (pa_sum);
				if (ranking[wheel] == particle) {
					continue;
				}
				if (sum >= r) {
					who = ranking[wheel];
					break;
				}
			}
			
			if (who == -1) {
				System.out.println(r+" and sum ="+sum+" and power: "+power);
				System.exit(1);
			}
			
			//who = chooseByRouletteWheelRankingBased_exponencial(); //chooseByRouletteWheelFitnessBased();
			// se vai conectar a um melhor, tudo bem.
			//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] == 0)) {
			//if ((ranking_of_particle[who] < ranking_of_particle[particle]) && (particle_failures[who] < particle_failures[particle])) {
			
			if ((particle_failures[who] == 0)/* && (ranking_of_particle[who] > ranking_of_particle[particle])*/) {
				if (!swarm_neighborhood_graph[who][particle]) {
							particle_failures[particle] = 0; 
				}
			} else {
				if (swarm_neighborhood_graph[who][particle]) {
					removeNeighbour(particle, who);
				}
			}		
					
			if (particle_failures[particle] == 0) {
				addNeighbour(who, particle);
			}
		}
	}

	private int numberOfNeighbours(int i) {
		int sum = 0;
		for (int neighbour = 0; neighbour < swarm_neighborhood_graph[i].length; neighbour++) {
			if (swarm_neighborhood_graph[i][neighbour]) {
				sum++;
			}
		}
		return sum;
	}

	private void updatePosition() {
		for (int particle = 0; particle < particle_position.length; particle++) {
			for (int dimension = 0; dimension < particle_position[particle].length; dimension++) {
				particle_position[particle][dimension] = particle_position[particle][dimension] + particle_velocity[particle][dimension];
				if (particle_position[particle][dimension] > PARTICLE_MAXX) {
					particle_position[particle][dimension] = PARTICLE_MAXX;
					particle_velocity[particle][dimension] = -particle_velocity[particle][dimension];
				}
				if (particle_position[particle][dimension] < PARTICLE_MINX) {
					particle_position[particle][dimension] = PARTICLE_MINX;
					particle_velocity[particle][dimension] = -particle_velocity[particle][dimension];
				}		
			}
		}
	}

	private void updateVelocity() {
		double c1 = 2.05f; 
		double c2 = 2.05f;
		double phi = c1 + c2;
		double factor = 2 / Math.abs(2 - phi - Math.sqrt(phi*phi - 4*phi));
		for (int particle = 0; particle < particle_velocity.length; particle++) {
			int best_neighbour = particle_best_neighbour[particle];
			// update particle with constriction factor
			for (int dimension = 0; dimension < particle_velocity[particle].length; dimension++) {
				double cognitive = c1*randomDouble()*(particle_best_value_position[particle][dimension] - particle_position[particle][dimension]);
				double social = c2*randomDouble()*(particle_best_value_position[best_neighbour][dimension] - particle_position[particle][dimension]);
				particle_velocity[particle][dimension] = factor * (particle_velocity[particle][dimension] + cognitive + social);
				// clamp the particles
				if (particle_velocity[particle][dimension] > PARTICLE_MAXV) {
					particle_velocity[particle][dimension] = PARTICLE_MAXV;
				} else if (particle_velocity[particle][dimension] < -PARTICLE_MAXV) {
					particle_velocity[particle][dimension] = -PARTICLE_MAXV;
				}
			}
		}
	}
	
	private double randomDouble() {
		return this.random_number_generator.nextDouble();
	}
	
	private double randomDoubleIndependent() {
		return this.random_number_generator_independent.nextDouble();
	}
	
	public int updateGBest() {
		// assuming here that the problem is minimization!
		this.swarm_gbest = this.findGBest();
		return this.swarm_gbest;
	}
	
	public int findGBest() {
		// assuming here that the problem is minimization! 
		double minval = Double.MAX_VALUE;
		int gbest = 0;
		for (int particle = 0; particle < particle_best_value.length; particle++) {
			if (particle_best_value[particle] < minval) {
				minval = particle_best_value[particle];
				gbest = particle;
			}
		}
		return gbest;
	}
	
	private void updateNBest() {
		// assuming here that the problem is minimization! 
		for (int particle = 0; particle < swarm_neighborhood_graph.length; particle++) {
			
			double minval = Double.MAX_VALUE;
			
			for (int neighbour = 0; neighbour < swarm_neighborhood_graph[particle].length; neighbour++) {
				
				if (neighbour == particle) {
					continue;
				}
				
				if (swarm_neighborhood_graph[particle][neighbour]) {
					if (particle_best_value[neighbour] < minval) {
						minval = particle_best_value[neighbour];
						particle_best_neighbour[particle] = neighbour;
					}
				}
				
			}
		}
		
	}

	private void evaluatePositionAndUpdatePersonal() {
		double zero = 0.000001D;
		for (int particle = 0; particle < particle_position.length; particle++) {
			particle_current_fitness[particle] = FUNCTION.compute(particle_position[particle]); 
			if (particle_current_fitness[particle] < particle_best_value[particle]) {
				double delta = Math.abs(particle_current_fitness[particle] -  particle_best_value[particle]);
				if (delta < zero) {
					particle_failures[particle]++;
				} else {
					particle_failures[particle] = 0;
					particle_best_value[particle] = particle_current_fitness[particle];
					for (int dimension = 0; dimension < particle_best_value_position[particle].length; dimension++) {
						particle_best_value_position[particle][dimension] = particle_position[particle][dimension];
					}
				} 
			} else {
				particle_failures[particle]++;
			}
		}	
	}

	private void initializePSO() {
		for (int particle = 0; particle < particle_position.length; particle++) {
			
			initializeParticle(particle);
			
		}
		initializeTopology();
		
	}
	
	private void initializeParticle(int particle) {
		for (int dimension = 0; dimension < particle_position[particle].length; dimension++) {
			particle_position[particle][dimension] =  (PARTICLE_INITIAL_RANGE_R - PARTICLE_INITIAL_RANGE_L) * randomDouble() + PARTICLE_INITIAL_RANGE_L;
			particle_best_value_position[particle][dimension] = particle_position[particle][dimension]; 
			particle_velocity[particle][dimension] = PARTICLE_MAXV*randomDouble();
			
			if (randomDouble() > 0.5) {
				particle_velocity[particle][dimension] = -particle_velocity[particle][dimension];
			}
		}
		particle_best_value[particle] = Double.MAX_VALUE;
		particle_current_fitness[particle] = particle_best_value[particle];
		particle_failures[particle] = 0;
		particle_failures_power[particle] = -0.9;
		particles_failures_threshold_particular[particle] = particles_failures_threshold;
	}
	
	private void initializeTopology() {
		for (int i = 0; i < swarm_neighborhood_graph.length; i++) {
			for (int j = 0; j < swarm_neighborhood_graph[i].length; j++) {
				removeNeighbour(i, j);
			}
		}	
		
		switch (swarm_initial_topology) {
		case GLOBAL:
			for (int i = 0; i < swarm_neighborhood_graph.length; i++) {
				for (int j = 0; j < swarm_neighborhood_graph[i].length; j++) {
					if (i != j) {
						addNeighbour(i, j);
					}
				}
			}		
			break;
		case RING:
			for (int i = 0; i < swarm_neighborhood_graph.length; i++) {
				addNeighbour(i, (i + 1) % swarm_neighborhood_graph.length);
			}		
			break;
		case VON_NEUMANN:
			int columns = swarm_neighborhood_graph.length / 2;	
			//     |   |   |   |
			//   - a - b - c - d -
			//     |   |   |   | 
			//   - e - f - g - h -
			//     |   |   |   |
			for (int j = 0; j < swarm_neighborhood_graph.length; j++) {
				int line = j / columns;
				addNeighbour(j, (j+1)%columns + (j/columns)*columns);
				addNeighbour(j, ((j+1)%columns + (j/columns + 1)*columns)%swarm_neighborhood_graph.length);
			}
			break;
		case RANDOM:
			for (int i = 0; i < swarm_neighborhood_graph.length; i++) {
				for (int j = 0; j < swarm_neighborhood_graph.length; j++) {
					if (i != j) {
						if (randomDouble() <= swarm_random_topology_p) {
							addNeighbour(i, j);
						}
					}
				}	
			}
			break;
		case THREESOME_PARTNERS:
			for (int i = 0; i < swarm_neighborhood_graph.length; i += 3) {
				addNeighbour(i, (i+1)%swarm_neighborhood_graph.length);
				addNeighbour(i, (i+2)%swarm_neighborhood_graph.length);
				addNeighbour((i+1)%swarm_neighborhood_graph.length, (i+2)%swarm_neighborhood_graph.length);
			}
			break;
		case NSOME_PARTNERS:
			int n = swarm_number_of_clusters;
			for (int i = 0; i < swarm_neighborhood_graph.length; i += n) {
				for (int p_i = 0; p_i < n; p_i++) {
					for (int p_j = p_i; p_j < n; p_j++) {
						if ((i+p_i)%swarm_neighborhood_graph.length != (i+p_j)%swarm_neighborhood_graph.length) {
							addNeighbour((i+p_i)%swarm_neighborhood_graph.length, (i+p_j)%swarm_neighborhood_graph.length);
						}
					}
					
					//addNeighbour(i, (i+p_i)%swarm_neighborhood_graph.length);
				}
				//addNeighbour(i, (i+1)%swarm_neighborhood_graph.length);
				//addNeighbour(i, (i+2)%swarm_neighborhood_graph.length);
				//addNeighbour((i+1)%swarm_neighborhood_graph.length, (i+2)%swarm_neighborhood_graph.length);
			}
			break;
		default:
			break;
		}	
	}

	private void addNeighbour(int i, int j) {
		if (i < NUMBER_OF_PARTICLES && j < NUMBER_OF_PARTICLES) {
			swarm_neighborhood_graph[i][j] = true;
			swarm_neighborhood_graph[j][i] = true;
			//swarm_viewer.addEdge(i, j, false);
			if (!(this.swarm_neigbourhood_graph_jung.containsEdge(i+""+j) || this.swarm_neigbourhood_graph_jung.containsEdge(j+""+i))) {
				this.swarm_neigbourhood_graph_jung.addEdge(i+"-"+j, (Integer) i, (Integer) j);
			}
		}
	}
	
	private void removeNeighbour(int i, int j) {
		if (i < NUMBER_OF_PARTICLES && j < NUMBER_OF_PARTICLES) {
			swarm_neighborhood_graph[i][j] = false;
			swarm_neighborhood_graph[j][i] = false;
                        //swarm_viewer.removeEdge(i, j, false);
			this.swarm_neigbourhood_graph_jung.removeEdge(i+"-"+j);
			this.swarm_neigbourhood_graph_jung.removeEdge(j+"-"+i);
		}
	}
	
	
	private double pa(int n, double p) {
		double sum = 0;
		int i;
		for (i = 1; i <= n; i++) {
			sum += Math.pow((double) i, p);
		}
		return sum;
	}
	
	int[][] getInfluenceGraph() {
		return getInfluenceGraph(1);
	}
	
	
	int[][] getInfluenceGraph(int k) {
		int[][] influence = null;
		if (this.swarm_influence_graph == null) {
			influence = new int[NUMBER_OF_PARTICLES][NUMBER_OF_PARTICLES];
			for (int i = 0; i < influence.length; i++) {
				for (int j = 0; j < influence[i].length; j++) {
					influence[i][j] = 0;
				}
			}
			for (int i = 0; i < influence.length; i++) {
				influence[particle_best_neighbour[i]][i] = 1;
			}
		} else {
			if (k == 1) {
				influence = this.swarm_influence_graph.clone();
			} else {
				influence = new int[NUMBER_OF_PARTICLES][NUMBER_OF_PARTICLES];
				for (int i = 0; i < influence.length; i++) {
					for (int j = 0; j < influence[i].length; j++) {
						influence[i][j] = 0;
					}
				}

				int steps_back = Math.abs(k); // 0 is the current iteration, just a trick to the loop above works right
				
				if (steps_back > this.current_iteration) {
					steps_back = this.current_iteration + 1;
				}
				
				for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle ++) {
					for (int history = 0; history < steps_back; history++) {
						int history_step = (this.current_iteration - history)%swarm_influence_history_maximum_length;
						influence[swarm_influence_history[particle][history_step]][particle]++;
					}
				}
			}
		}
		return influence;
	}
}