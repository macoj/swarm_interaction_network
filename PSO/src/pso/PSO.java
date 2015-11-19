package pso;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import benchmark_functions.*;

//import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;
import java.util.logging.Level;
import java.util.logging.Logger;

public class PSO implements Runnable{
	
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
	double particle_best_value_position[][];
	double particle_best_value[];
	int particle_best_neighbour[];
	int particle_failures[];
	double particle_failures_power[];
	
	double particle_fev_1[];
	double particle_fev_2[];
	double particle_fev_2_reduced[];
	double swarm_distance_particles[][];
	double[][] swarm_centroids = null;
	double[][] swarm_clusters_kohonen = null;
	int[] swarm_clusters = null;
	
	boolean swarm_viewer_enabled = false;
	boolean network_science_metrics = false;
	
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
	
	enum TOPOLOGY_MECHANISM { NONE, BARABASI_BASED, BARABASI_BASED_DYNAMIC_A, BARABASI_BASED_DYNAMIC_B, BARABASI_BASED_DYNAMIC_C, BARABASI_BASED_DYNAMIC_D, BARABASI_BASED_DYNAMIC_DISSERTACAO, BARABASI_BASED_DYNAMIC_LOG };
	
	static final int swarm_topology_sleep = 1;
	
	int particles_failures_threshold;
	int[] particles_failures_threshold_particular = null;
	
	double sphere (double x[]) {
		double sum = 0.0f;
		for (int i = 0; i < x.length; i++) {
			sum += x[i]*x[i];
		}
		sum = (double) Math.sqrt(sum);
		return sum;
	}
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
		
		int evaluation = 300000;
		
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
		particle_best_neighbour = new int[NUMBER_OF_PARTICLES];
		particle_failures = new int[NUMBER_OF_PARTICLES];
		
		particle_failures_power = new double[NUMBER_OF_PARTICLES];
		particles_failures_threshold_particular = new int[NUMBER_OF_PARTICLES];
		network_science_metrics = false;
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

		
		//		swarm_initial_topology = TOPOLOGY.CLAN;
		//swarm_initial_topology = TOPOLOGY.NO_TOPOLOGY;
		swarm_gbest = -1;
		
//		swarm_viewer = new SwarmViewer_graphstream();
		//swarm_viewer = new GraphView();
		//swarm_fitness_plotter = new Plotter_prefuse();
		//swarm_fev_plotter = new Plotter_prefuse();
		swarm_number_of_clusters = 3;
		
		swarm_viewer_enabled = false;
		
		swarm_initial_maximum_neighbors = 3;
		
		swarm_distance_particles = new double[NUMBER_OF_PARTICLES][NUMBER_OF_PARTICLES];
		this.particle_fev_1 = new double[NUMBER_OF_PARTICLES];
		this.particle_fev_2 = new double[NUMBER_OF_PARTICLES];
		this.particle_fev_2_reduced = new double[NUMBER_OF_PARTICLES];
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
	
	
	public int[] k_means(int k) {
		
		if (swarm_clusters == null) {
			swarm_clusters = new int[NUMBER_OF_PARTICLES];
			if (swarm_centroids == null) {
				swarm_centroids = new double[k][DIMENSION];
				for (int cluster = 0; cluster < k; cluster++) {
					for (int d = 0; d < DIMENSION; d++) {
						swarm_centroids[cluster][d] = (PARTICLE_INITIAL_RANGE_R - PARTICLE_INITIAL_RANGE_L) * randomDouble() + PARTICLE_INITIAL_RANGE_L;
					}
				}
			}
		}
		
		if (swarm_centroids == null) {
			// baricentros
			swarm_centroids = new double[k][DIMENSION];
			int[] count = new int[k];
			for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
				for (int dimension = 0; dimension < DIMENSION; dimension++) {
					swarm_centroids[swarm_clusters[particle]][dimension] += particle_position[particle][dimension];
				}
				count[swarm_clusters[particle]]++;
			}
			
			for (int cluster = 0; cluster < k; cluster++) {
				for (int dimension = 0; dimension < DIMENSION; dimension++) {
					swarm_centroids[cluster][dimension] /= count[cluster];
				}	
			}
		}
		
		for (int p = 0; p < swarm_clusters.length; p++) {
			swarm_clusters[p] = -1;
		}
		

		
		boolean change = false;
		do {
			change = false;
			for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
				double distance = Double.MAX_VALUE;
				int own_cluster = -1;
				for (int cluster = 0; cluster < k; cluster++) {
					double current_distance = euclidianDistance(swarm_centroids[cluster], particle_position[particle]);
					 if (distance > current_distance) {
						 distance = current_distance;
						 own_cluster = cluster;
					 }
				}
				if (swarm_clusters[particle] != own_cluster) {
					change = true;
				}
				swarm_clusters[particle] = own_cluster;
			}
			// recalcula centroids
			double[][] new_centroids = new double[k][DIMENSION];
			int[] count = new int[k];
			
			for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
				for (int d = 0; d < DIMENSION; d++) {
					new_centroids[swarm_clusters[particle]][d] += particle_position[particle][d];
				}
				count[swarm_clusters[particle]]++;
			}
			for (int cluster = 0; cluster < k; cluster++) {
				for (int d = 0; d < DIMENSION; d++) {
					swarm_centroids[cluster][d] = new_centroids[cluster][d]/((double)count[cluster]);
				}	
			}
		} while (change == false);
		
		return swarm_clusters;
	}
	
	private double euclidianDistance(double[] ds, double[] ds2) {
		double sum = 0.0f;
		for (int i = 0; i < ds2.length; i++) {
			sum += Math.pow(ds[i] - ds2[i], 2);
		}
		sum = Math.sqrt(sum);
		return sum;
	}

	
	public static void main(String[] args) throws FileNotFoundException {
		/*double a[] = { 1, 1};
		double b[] = { -1, -1};
		PSO pso = new PSO(args);
		double scalar = pso.scalarProduct(a, b);
		System.out.println("cos 0 = "+pso.cosBetweenVectors(a, b));*/
	
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
		
//		String filename = "psoResults" + swarm_initial_topology.toString();
//		if (swarm_topology_mechanism != null){
//			filename += "_"+ swarm_topology_mechanism.toString();
//		}
//		filename += ".txt";
//                File psoFile = new File(filename);
//                System.out.println("File being created at " + psoFile.getAbsolutePath() +"...");
//            
//                PrintWriter printWriter = new PrintWriter(psoFile)) {
                Diversity diversity = new Diversity(this, printWriter);
                DAPSO dapso = new DAPSO(this, printWriter);
                PsoDD psoDD = new PsoDD(this, printWriter);
                
                
                double run_final_values[] = new double[RUNS];
                for (int run = 0; run < RUNS; run++) {
                    this.initializeRNG();
                    this.initializePSO();
                    
                    //this.kohonen(swarm_number_of_clusters);

                    /// plotters:
                    //swarm_fitness_plotter.setLabel("Fitness");
                    //swarm_fitness_plotter.init();
                    //swarm_fev_plotter.setLabel("evolutionary factor");
                    //swarm_fev_plotter.init();
                    
                    double[] eigenvalues = null;
                    double[] eigenvaluesADJ = null;
                    double[] eigenvaluesADJTOP = null;
                    //this.k_means(swarm_number_of_clusters);
                    double centroids_sum = -1;
                    do {
                        if(current_iteration%100==0)
                            System.out.println(current_iteration/100 + " de 30");
                        this.evaluatePositionAndUpdatePersonal();
                        this.swarm_gbest = this.getGBest();
                        this.findBestNeighbours();
                        //this.updateParticlesNeighboursDistances();
                        //this.updateParticlesDistances();
                        //this.updateParticlesEvolutionaryState_2(true);
                        //this.updateParticlesEvolutionaryState_2(false);
                        //this.updateEvolutionaryState_1();
                        this.updateNeighbourhood();
                        this.updateVelocity();
                        this.updatePosition();
                        this.updateInfluenceGraph();
                        
                        String it = "it:#" + this.current_iteration + " " + particle_best_value[this.swarm_gbest];
//				System.out.println("it:#" + this.current_iteration+ " " + particle_best_value[this.swarm_gbest]);
                        printWriter.println(it);
                        
                        String igLine = "ig:#" + this.current_iteration + " ";
//				System.out.print("ig:#" + this.current_iteration + " ");
                        
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
                            //System.out.print(this.particle_best_value[i] + " ");
                            
                            printWriter.print(this.particle_best_value[i] + " ");
                            
                        }
                        
                        printWriter.println();
                        
                        diversity.iterate();
                        dapso.iterate();
                        psoDD.iterate();
                        
//				PSODD METRIC PRINT
//				
//				System.out.println("Iteration: " + this.current_iteration);
//				System.out.println("Best particle: " + particle_best_value[this.swarm_gbest]);
//				System.out.println("Previous best particle: " + psoDD.getPreviousFbest());
//				System.out.println("Current average velocity: " + psoDD.getAvg_particles_velocity());
//				System.out.println("Previous average velocity: " + psoDD.getPrevious_avg_particles_velocity());
//				System.out.println("Stagnation Ratio: " + psoDD.getR() +"\n");
				
//				DAPSO METRIC PRINT
                        
                        
                        
                        
                        //System.out.print("pbest:#" + this.current_iteration + " ");
                        
                        
                        
                        
                        //System.out.println();
                        //              
//				tesst = this.getInfluenceGraph(this.swarm_influence_history_maximum_length);
//				for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
//					for (int j = 0; j < NUMBER_OF_PARTICLES; j++) {
//						System.out.print(test[i][j]+" ");
//					}
//					//System.out.println();
//				}
                        //System.out.println();
                        
                        /// particle position
//				System.out.println();
                        printWriter.println();
                        for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
                            //System.out.print("position:#" + this.current_iteration + "#");
                            
                            //                      printWriter.print("position:#" + this.current_iteration + "#");
                            
                            //System.out.print(i + " ");
                            //                  printWriter.print(i + " ");
                            
                            for (int d = 0; d < DIMENSION; d++) {
                                //System.out.print(this.particle_position[i][d] + " ");
                                
                                //                          printWriter.print(this.particle_position[i][d] + " ");
                                
                            }
                            //System.out.println();
                            //                  printWriter.println();
                        }
                        //System.out.println();
                        //              printWriter.println();
                        
                        {
                            int ranking[] = new int[NUMBER_OF_PARTICLES];
                            int ranking_of_particle[] = new int[NUMBER_OF_PARTICLES];
                            int particle_index_shuffled[] = new int[NUMBER_OF_PARTICLES];
                            for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
                                ranking[particle] = particle;
                                particle_index_shuffled[particle] = particle;
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
                            
//					if (this.current_iteration == 10) {
//						InfluenceGraph_graphview.toPNG(swarm_neighborhood_graph, swarm_influence_graph, ranking_of_particle);
//						//System.exit(0);
//					}
                            
                        }
                        
                        //this.updateInfluenceGraph(iteration, 1); // just a 'k' window
                        //this.kohonen(swarm_number_of_clusters);
                        /*
                        try {
                        System.out.println("natural conectivity: "+NaturalConnectivity.getInstance().calculate(getAdjacency(this.swarm_neighborhood_graph)));
                        } catch (Exception e) {
                        System.out.println("natural conectivity: 0.0");
                        }
                        
                        try {
                        System.out.println("Algebraic conectivity: "+AlgebraicConnectivity.getInstance().calculate(getAdjacency(this.swarm_neighborhood_graph)));
                        
                        } catch (Exception e) {
                        System.out.println("algebraic conectivity: 0.0");
                        }
                        */
                        
                        //System.out.println("it:#"+run+" "+ this.current_iteration + " " + decimal_format.format(particle_best_value[swarm_gbest]));
                        //System.out.println("it:#"+run+" "+ this.current_iteration + " fev1: " + decimal_format.format(swarm_evolutionary_state_1));
                        //System.out.println("it:#"+run+" "+ this.current_iteration + " fev2: " + decimal_format.format(swarm_evolutionary_state_2));
                        //System.out.println("it:#"+run+" "+ iteration + " fev2r: " + decimal_format.format(swarm_evolutionary_state_2_reduced));
                        //topologyViewerSwapGBest(gbest);
                        //swarm_fitness_plotter.addPoint(iteration, particle_best_value[gbest]);
                        if (((this.current_iteration+1) % 100) == 0) {
                            //k_means(10);//	k_meanNeighbourhoodCreate(10);
                            //initializeTopology();
                            swarm_clusters = null;
                            swarm_clusters_kohonen = null;
                            swarm_centroids = null;
                            //this.kohonen(swarm_number_of_clusters);
                            
                        }
                        /*
                        System.out.print("run:#"+run+"it: "+iteration+" deg: ");
                        for (int vertex = 0; vertex < NUMBER_OF_PARTICLES; vertex++) {
                        System.out.print(numberOfNeighbours(vertex)+" ");
                        }
                        System.out.println();
                        
                        System.out.print("run:#"+run+"it: "+iteration+" fev_p1: ");
                        for (int vertex = 0; vertex < NUMBER_OF_PARTICLES; vertex++) {
                        System.out.print(particle_fev_2[vertex]+" ");
                        }
                        System.out.println();
                        
                        System.out.print("run:#"+run+"it: "+iteration+" fev_p2: ");
                        for (int vertex = 0; vertex < NUMBER_OF_PARTICLES; vertex++) {
                        System.out.print(particle_fev_2_reduced[vertex]+" ");
                        }
                        System.out.println();
                        */
                        /*
                        System.out.println("edge count: "+this.swarm_neigbourhood_graph_jung.getEdgeCount());
                        Map<Integer, Double> average_distances = getAverageDistances(this.swarm_neigbourhood_graph_jung);
                        double average_distance = 0.0;
                        for (int vertex = 0; vertex < NUMBER_OF_PARTICLES; vertex++) {
                        double vertex_average_distance = average_distances.get((Integer) vertex);
                        System.out.print(vertex_average_distance+" ");
                        average_distance += vertex_average_distance;
                        }
                        average_distance /= (double) NUMBER_OF_PARTICLES;
                        System.out.println("\naverage_distance: "+average_distance);
                        
                        Map<Integer, Double> clustering = Metrics.clusteringCoefficients(this.swarm_neigbourhood_graph_jung);
                        double average_clustering = 0.0;
                        for (int vertex = 0; vertex < NUMBER_OF_PARTICLES; vertex++) {
                        double clustering_factor = clustering.get((Integer) vertex);
                        System.out.print(clustering_factor+" ");
                        average_clustering += clustering_factor;
                        }
                        average_clustering /= (double) NUMBER_OF_PARTICLES;
                        System.out.println("\naverage_clustering: "+average_clustering);
                        */
                        /*
                        try {
                        Thread.sleep(50);
                        } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                        }
                        */
                        //this.updateEvolutionaryState();
                        //swarm_fev_plotter.addPoint(iteration, swarm_evolutionary_state);
                        this.current_iteration = this.current_iteration + 1;
                        try {
                            Thread.sleep(1);
                        } catch (InterruptedException ex) {
                            Logger.getLogger(PSO.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    } while (this.current_iteration < MAXITER);
                    //System.out.println("MAXITER: " + MAXITER);
                    //System.out.println("it:#"+run+" " + iteration + " " + particle_best_value[this.swarm_gbest]);
                    run_final_values[run] = particle_best_value[this.swarm_gbest];
                    System.out.println("The file is ready.");
                    //System.out.println("random_number_generator_seed: "+random_number_generator_seed);
                }
                printWriter.close();
                //printFinals(run_final_values);
                //printAvgStdDev(run_final_values);
            
	}

	private void updateInfluenceGraphHistory(int iteration) {
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			swarm_influence_history[particle][iteration%swarm_influence_history_maximum_length] = particle_best_neighbour[particle];
		}
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
	/*
	private void updateInfluenceGraph(int iteration, int steps_back) {
		if (swarm_influence_graph_weighted) {
			if (swarm_influence_history != null) {
				for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
					for (int j = 0; j < NUMBER_OF_PARTICLES; j++) {
						swarm_influence_graph[i][j] = 0;
					}
				}
			}
		}
		
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			
			if (swarm_influence_graph_weighted) {
				if (swarm_influence_history == null) {
					swarm_influence_graph[particle_best_neighbour[particle]][particle]++;// [i][j] = i sends information to j
				} else {
					
					steps_back++; // 0 is the current iteration, just a trick to the loop above works right
					
					if (steps_back < 0) {
						steps_back = Math.abs(steps_back);
					} 
										
					if ((steps_back > iteration)) {
						steps_back = iteration;
					}
					
					for (int history = 0; history < steps_back; history++) {
						int history_step = (iteration - history)%swarm_influence_history_maximum_length;
						swarm_influence_graph[swarm_influence_history[particle][history_step]][particle]++;
					}
				}
			}	
		}
		
	}
	*/
	
	private int[][] toUndirect(boolean[][] influenceGraph) {
		int[][] undirect = new int[influenceGraph.length][influenceGraph.length];
		for (int i = 0; i < influenceGraph.length; i++) {
			for (int j = 0; j < influenceGraph.length; j++) {
				if (influenceGraph[i][j] || influenceGraph[j][i]) {
					undirect[i][j] = 1;
					undirect[j][i] = 1;
				} else {
					undirect[i][j] = 0;
					undirect[j][i] = 0;
				}
			}
		}
		return undirect;
	}
	
	private int[][] toUndirect(int[][] influenceGraph) {
		int[][] undirect = new int[influenceGraph.length][influenceGraph.length];
		for (int i = 0; i < influenceGraph.length; i++) {
			for (int j = 0; j < influenceGraph.length; j++) {
				int value = influenceGraph[i][j] >= influenceGraph[j][i] ? influenceGraph[i][j] : influenceGraph[j][i];
				undirect[i][j] = value;
				undirect[j][i] = value;
			}
		}
		return undirect;
	}
	
	private Map<Integer, Double> getAverageDistances(Graph<Integer, String> graph) {
		DijkstraDistance<Integer,String> dijkstra_distances = new DijkstraDistance<Integer, String>(graph);
		
		Map<Integer, Double> map = new HashMap<Integer, Double>();
		for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
			double average_size = 0.0;
			int neighbourhood_size = 0;
			Map<Integer, Number> map_distance_from_i = dijkstra_distances.getDistanceMap((Integer) i);
			for (int j = 0; j < NUMBER_OF_PARTICLES; j++) {
				if (map_distance_from_i.get((Integer) j) != null){
					average_size += map_distance_from_i.get((Integer) j).doubleValue();
					neighbourhood_size++;
				}
			}
			map.put((Integer) i, (Double) (average_size/(neighbourhood_size-1)));
		}
		return map;
	}
	private int[][] getAdjacency(boolean[][] swarm_neighborhood_graph2) {
		int[][] ret = new int[swarm_neighborhood_graph2.length][swarm_neighborhood_graph2.length];
		for (int i = 0; i < ret.length; i++) {
			for (int j = 0; j < ret.length; j++) {
				if (swarm_neighborhood_graph2[i][j]) {
					ret[i][j] = 1; 
				} else {
					ret[i][j] = 0;
				}
			}
		}
		return ret;
	}
	private void updateParticlesNeighboursDistances() {
		for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
			for (int j = i; j < NUMBER_OF_PARTICLES; j++) {
				if (swarm_neighborhood_graph[j][i]) {
					swarm_distance_particles[i][j] = euclidianDistance(particle_position[i], particle_position[j]);
					swarm_distance_particles[j][i] = swarm_distance_particles[i][j];
				}
			}
			
		}
	}
	
	private double matrixNorma(double [] matrix) {
		double sum = 0.0;
		for (int element = 0; element < matrix.length; element++) {
			sum += Math.pow(matrix[element], 2.0);
		}
		sum = Math.sqrt(sum);
		return sum;
	}
	private double cosBetweenVectors(double[] vectorA, double[] vectorB) {
		double cos = this.scalarProduct(vectorA, vectorB); 
		cos = cos/(matrixNorma(vectorA)*matrixNorma(vectorB));
		return cos;
	}
	private double scalarProduct(double[] matrixA, double[] matrixB) {
		double sum = 0.0;
		
		if (matrixA.length == matrixB.length) {
			for (int element = 0; element < matrixB.length; element++) {
				sum += matrixA[element]*matrixB[element];
			}
		}
		// (-1, 2) (1, -2)
		return sum;		
	}
	/*
	 private void updateParticlesEvolutionaryState() {
		DecimalFormat df = new DecimalFormat("#.##");
		// we admit that updateParticlesDistances() was called here.
		int number_of_connections = 0;
		int number_of_neighbours[] = new int[NUMBER_OF_PARTICLES];
		this.swarm_evolutionary_state = 0d;
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			double d_g = 0d;
			double d_min = Double.MAX_VALUE;
			double d_max = Double.MIN_VALUE;
			this.particle_fev[particle] = 0.0d;
			number_of_neighbours[particle] = 0;
			for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
				if (particle == neighbour)
					continue;
				if (swarm_neighborhood_graph[particle][neighbour]) {
					if (this.swarm_distance_particles[particle][neighbour] > d_max) {
						d_max = this.swarm_distance_particles[particle][neighbour];
					} 
					if (this.swarm_distance_particles[particle][neighbour] < d_min) {
						d_min = this.swarm_distance_particles[particle][neighbour];
					}
					if (particle_best_neighbour[particle] == neighbour) {
						if (particle_best_neighbour[particle] == particle)
							System.out.println("QUE ZICA EH ESSA!?!?!?!?");
						d_g = this.swarm_distance_particles[particle][neighbour];
					}
					//this.particle_fev[particle] += this.swarm_distance_particles[particle][neighbour];
					number_of_neighbours[particle]++;
					number_of_connections++;
				}
			}
			//this.particle_fev[particle] /= (double) number_of_neighbours[particle];
			this.particle_fev[particle] = (d_g - d_min)/(d_max - d_min);
		}
		
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			this.swarm_evolutionary_state += (this.particle_fev[particle]*((double) number_of_neighbours[particle]))/((double) number_of_connections);
			//System.out.print(df.format(this.particle_fev[particle])+"*"+number_of_neighbours[particle]+" ");
			if (this.particle_fev[particle] < 0.01) {
				//removeNeighbourCarefully(particle, particle_best_neighbour[particle]);
				particle_failures[particle] += 25;
				
			}
		}
		//System.out.println();
		
	}
	 */
	
	private void updateParticlesEvolutionaryState_2(boolean reduced) {
		// we admit that updateParticlesDistances() was called here.
		int number_of_connections = 0;
		int number_of_neighbours[] = new int[NUMBER_OF_PARTICLES];
		if (reduced) {
			this.swarm_evolutionary_state_2_reduced = 0d;
		} else {
			this.swarm_evolutionary_state_2 = 0d;
		}
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			int best_neighbour = particle_best_neighbour[particle];
			double d_g = this.swarm_distance_particles[best_neighbour][particle];
			double d_min = Double.MAX_VALUE;
			double d_max = Double.MIN_VALUE;	
			double fev = 0.0;
			number_of_neighbours[particle] = 0;
			for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
				if (swarm_neighborhood_graph[particle][neighbour]) {
					double distance = 0.0;
					
					// distancia dessa particula para o melhor vizinho da particula:
					distance = this.swarm_distance_particles[best_neighbour][neighbour];
					
					if (!swarm_neighborhood_graph[best_neighbour][neighbour]) {
						if (reduced) {
							distance = d_g + this.swarm_distance_particles[particle][neighbour];
						}
					}
					
					if (neighbour == best_neighbour) {
						distance = d_g;
					}
					
					if (distance >= d_max) {
						d_max = distance;
					} 
					
					if (distance <= d_min) {
						d_min = distance;
					}
					
					number_of_neighbours[particle]++;
					number_of_connections++;
				}
			}
			if (d_max == d_min) {
				fev = 1.0d;
			} else {
				fev = (d_g - d_min)/(d_max - d_min);
			}		
			if (fev < 0.0) {
				System.out.println("d_g "+d_g+"d_min"+d_min+"d_max: "+d_max);
				System.out.println("the best neigh of "+particle+"is "+best_neighbour);;
				System.exit(0);
			}
			
			if (reduced) {
				this.particle_fev_2_reduced[particle] = fev;
			} else {
				this.particle_fev_2[particle] = fev;
			}
		}
		if (reduced) {
			for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
				this.swarm_evolutionary_state_2_reduced += (this.particle_fev_2_reduced[particle]*((double) number_of_neighbours[particle]))/((double) number_of_connections);
				//this.swarm_evolutionary_state += (this.particle_fev[particle])/((double) NUMBER_OF_PARTICLES);
				//System.out.print(df.format(this.particle_fev[particle])+"*"+number_of_neighbours[particle]+" ");
				//System.out.print(this.particle_fev[particle]+"*"+number_of_neigbours[particle]+" ");
				if (this.particle_fev_2_reduced[particle] < 0.05) {
					//removeNeighbourCarefully(particle, particle_best_neighbour[particle]);
					//System.out.println(df.format(this.particle_fev[particle])+"*"+number_of_neighbours[particle]+" ");
					//particle_failures[particle] += 10;
					
				}
			}
		} else {
			for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
				this.swarm_evolutionary_state_2 += (this.particle_fev_2[particle]*((double) number_of_neighbours[particle]))/((double) number_of_connections);
				//this.swarm_evolutionary_state += (this.particle_fev[particle])/((double) NUMBER_OF_PARTICLES);
				//System.out.print(df.format(this.particle_fev[particle])+"*"+number_of_neighbours[particle]+" ");
				//System.out.print(this.particle_fev[particle]+"*"+number_of_neigbours[particle]+" ");
				if (this.particle_fev_2[particle] < 0.05) {
					//removeNeighbourCarefully(particle, particle_best_neighbour[particle]);
					//System.out.println(df.format(this.particle_fev[particle])+"*"+number_of_neighbours[particle]+" ");
					//particle_failures[particle] += 10;
					
				}
			}
			
		}
	}

	
	private boolean removeNeighbourCarefully(int particle, int i) {
		boolean removed = false;
		if (numberOfNeighbours(i) > 2) {
			if (numberOfNeighbours(particle) > 2) {
				removeNeighbour(particle, i);
				removed = true;
			}
		}
		return removed;
		
	}
	private void updateParticlesDistances() {
		for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
			for (int j = i + 1; j < NUMBER_OF_PARTICLES; j++) {
				swarm_distance_particles[i][j] = euclidianDistance(particle_position[i], particle_position[j]);
				swarm_distance_particles[j][i] = swarm_distance_particles[i][j];
			}
			
		}
		
	}
	private void printAvgStdDev(double[] run_final_values) {
		double average = 0.0;
		double std_dev = 0.0;
		for (int run = 0; run < run_final_values.length; run++) {
			average += run_final_values[run];
		}
		average /= (double) run_final_values.length;
		
		for (int run = 0; run < run_final_values.length; run++) {
			std_dev += (run_final_values[run] - average)*(run_final_values[run] - average);
		}
		std_dev /= (double) (run_final_values.length - 1);
		std_dev = Math.sqrt(std_dev);
		System.out.println("avg: "+average);
		System.out.println("std_dev: "+std_dev);
	}
	private void printFinals(double[] run_final_values) {
		for (int run = 0; run < run_final_values.length; run++) {
			System.out.println("final#"+run+" "+run_final_values[run]);
		}
		
	}
	private int[] kohonen(int neuronNumber) {
		boolean changed;
		
		double learning_factor = 0.01;
		double initial_neighborhood_size = 1;
		double initial_learning_factor = 0.1;
		double time_constant_1 = 1000/Math.log(initial_neighborhood_size);
		double time_constant_2 = 1000;
		
		if (swarm_clusters == null) {
			swarm_clusters = new int[NUMBER_OF_PARTICLES];
		}
		
		int iteration = 1;
		if (swarm_clusters_kohonen == null) {
			swarm_clusters_kohonen = new double[neuronNumber][DIMENSION];
			// initialize them
			for (int i = 0; i < swarm_clusters_kohonen.length; i++) {
				for (int j = 0; j < swarm_clusters_kohonen[i].length; j++) {
					swarm_clusters_kohonen[i][j] = randomDouble();
				}
			}
		}
		
		do {
			changed = false;
			learning_factor = initial_learning_factor*Math.exp(-((double) iteration)/time_constant_2);
			for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
				// competition
				int neuron_winner = getTheNearest(particle_position[particle], swarm_clusters_kohonen);
				// cooperation
				
				// adaptation
				for (int neuron = 0; neuron < swarm_clusters_kohonen.length; neuron++) {
					// topological neighborhood
					double distance_to_winner = euclidianDistance(swarm_clusters_kohonen[neuron], swarm_clusters_kohonen[neuron_winner]);
					double neihborhood_size = initial_neighborhood_size*Math.exp(-((double) iteration)/time_constant_1);
					double topological_neighborhood = Math.exp(-(Math.pow(distance_to_winner, 2))/(2*Math.pow(neihborhood_size, 2)));
					for (int dimension = 0; dimension < swarm_clusters_kohonen[neuron].length; dimension++) {
						double delta_to_winner = swarm_clusters_kohonen[neuron_winner][dimension] - swarm_clusters_kohonen[neuron][dimension]; 
						swarm_clusters_kohonen[neuron][dimension] += learning_factor*topological_neighborhood*delta_to_winner;
					}
				}
				if (swarm_clusters[particle] != neuron_winner) {
					swarm_clusters[particle] = neuron_winner;
					changed = true;
				}
//				
			} 
			iteration++;
/*			System.out.println("another iteration... " + iteration);
			for (int clusters = 0; clusters < NUMBER_OF_PARTICLES; clusters++) {
				System.out.print(" " + swarm_clusters[clusters]);
			}
			System.out.println();
*/		} while (changed);
		
		return swarm_clusters;
	}

	private int getTheNearest(double[] ds, double[][] neuronsWeight) {
		int nearest = 0;
		double nearest_distance = euclidianDistance(ds, neuronsWeight[0]);
		for (int i = 1; i < neuronsWeight.length; i++) {
			double current_distance = euclidianDistance(ds, neuronsWeight[i]);
			if (current_distance < nearest_distance) {
				nearest = i;
				nearest_distance = current_distance;
			}
		}
		return nearest;
	}
	
	private void updateEvolutionaryState_1() {
		int particle;
		double d_min = Double.MAX_VALUE;
		double d_max = Double.MIN_VALUE;
		double d_g = 0;
		
		double[] mean_distance = new double[NUMBER_OF_PARTICLES];
		
		for (particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			int neighbour;	// fazer depois com apenas os vizinhos!
			// mean_distance[i] = (1/N-1) * sum ( euclidian distance from i to all other particles )
			mean_distance[particle] = 0.0;
			for (neighbour = 0; neighbour < NUMBER_OF_PARTICLES; ++neighbour) {
				if (particle == neighbour) {
					continue;
				}
				mean_distance[particle] += swarm_distance_particles[particle][neighbour];
			}
	
			mean_distance[particle] /= NUMBER_OF_PARTICLES - 1;
	
			if (mean_distance[particle] > d_max) {
				d_max = mean_distance[particle];
			}
	
			if (mean_distance[particle] < d_min) {
				d_min = mean_distance[particle];
			}
		}
	
		d_g = mean_distance[swarm_gbest];
		swarm_evolutionary_state_1 = (d_g - d_min)/(d_max - d_min);
	
	}

	private void updateNeighbourhood() {
		if (swarm_topology_mechanism != null) {
			switch (swarm_topology_mechanism) {
			case BARABASI_BASED_DYNAMIC_A:
				this.updateNeighbourhood_BA_DYNAMIC_A();
				break;
			case BARABASI_BASED_DYNAMIC_B:
				this.updateNeighbourhood_BA_DYNAMIC_B();
				break;
			case BARABASI_BASED_DYNAMIC_C:
				this.updateNeighbourhood_BA_DYNAMIC_C();
				break;			
			case BARABASI_BASED_DYNAMIC_D:
				this.updateNeighbourhood_BA_DYNAMIC_D();
				break;
			case BARABASI_BASED_DYNAMIC_DISSERTACAO:
				this.updateNeighbourhood_BA_DYNAMIC_DISSERTACAO();
				break;
			case BARABASI_BASED_DYNAMIC_LOG:
				this.updateNeighbourhood_BA_DYNAMIC_LOG();
				break;			
			default:
				break;
			}
		}
		
	}
	
	
private void updateNeighbourhood_BA_DYNAMIC_LOG() {
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
	// ranking[0] eh o index da PIOR particula
	// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
	// ranking_of_particle[i] eh a posicao da particula i no ranking, valor = 0 indica PIOR POSICAO 
	
	//random_number_generator = this.random_number_generator_independent;
	// in the same order of particle updating, the swarm behaviour is different
	// when different seeds are used. 
	
	//// -> esse laco pega uma particula aleatoriamente no enxame!
	for (int particle_ordered = 0; particle_ordered < ranking_of_particle.length; particle_ordered++) {
		
		int particle = particle_index_shuffled[particle_ordered];
		double r;// = random_number_generator.nextDouble(); 
		int who = 0;
		float sum = 0;
		double power = 1;
		
		/*
		if (particle_failures[particle] < particles_failures_threshold_particular[particle]) {
			continue;
		}
		*/
		///////// LOG
		double log_k = 0.01;
		double p_failures = (double) particle_failures[particle];
		double p_reconnect  = Math.max(0, Math.log10(log_k*p_failures)); /// min seria sempre 0 ou negativo
		r = random_number_generator.nextDouble();
		if (r > p_reconnect) {
			continue; // nao executa o mecanismo de reconexao
		}
 
					
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
	
private void updateNeighbourhood_BA_DYNAMIC_DISSERTACAO() {
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
	// ranking[0] eh o index da PIOR particula
	// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
	// ranking_of_particle[i] eh a posicao da particula i no ranking, valor = 0 indica PIOR POSICAO 
	
	//random_number_generator = this.random_number_generator_independent;
	// in the same order of particle updating, the swarm behaviour is different
	// when different seeds are used. 
	
	//// -> esse laco pega uma particula aleatoriamente no enxame!
	for (int particle_ordered = 0; particle_ordered < ranking_of_particle.length; particle_ordered++) {
		
		int particle = particle_index_shuffled[particle_ordered];
		
		if (particle_failures[particle] < particles_failures_threshold_particular[particle]) {
			continue;
		}

		double r;// = random_number_generator.nextDouble(); 
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

	// the order of the particles being updated changes the behaviour of
	// the swarm. (all the rng have the same seed here)
	/*
	it:#0 1500 1561306.5549688472   <- when the same rng is used to shuffle the particle updating order
	random_number_generator_seed: 1339596216600
	 */
	/*
	it:#0 1500 1374333.8417697118
	random_number_generator_seed: 1339596216600
	 */
	/*
	it:#0 1500 429.4389155150009
	random_number_generator_seed: 1339596216600
	 */
	private void updateNeighbourhood_BA_DYNAMIC_D() {
		
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
		// ranking[0] eh o index da PIOR particula
		// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
		// ranking_of_particle[i] eh a posicao da particula i no ranking
		
		//random_number_generator = this.random_number_generator_independent;
		// in the same order of particle updating, the swarm behaviour is different
		// when different seeds are used. 
		/*
		it:#0 1500 1561306.5549688472  
		random_number_generator_seed: 1339596216600
		 */
		/*
		it:#0 1500 77.22749784073
		random_number_generator_seed: 1339596216600
		*/
		/*
		it:#0 1500 1646246.3367618397
		random_number_generator_seed: 1339596216600 
		 */
		for (int particle_ordered = 0; particle_ordered < ranking_of_particle.length; particle_ordered++) {
			
			int particle = particle_index_shuffled[particle_ordered];
			
			if (particle_failures[particle] < particles_failures_threshold_particular[particle]) {
				continue;
			}
			//particles_failures_threshold_particular[particle] = particles_failures_threshold_particular[particle] - 1;
			if (particles_failures_threshold_particular[particle] == 0) {
				//particles_failures_threshold_particular[particle] = 50;
			}
			double r;// = random_number_generator.nextDouble(); 
			int who = 0;
			float sum = 0;
			double power = 1; //particle_failures_power[particle];
			//particle_failures_power[particle] += 1.0;
			
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
			}
			//who = chooseByRouletteWheelRankingBased_exponencial(); //chooseByRouletteWheelFitnessBased();
			// se vai conectar a um melhor, tudo bem.
			//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] == 0)) {
			//if ((ranking_of_particle[who] < ranking_of_particle[particle]) && (particle_failures[who] < particle_failures[particle])) {
			//if ((ranking_of_particle[who] < ranking_of_particle[particle])) {
			if (particle_failures[who] == 0) {
				if (!swarm_neighborhood_graph[who][particle]) {
					//if (particle_failures[who] == 0) {
						//	addNeighbour(particle, who);
							particle_failures[particle] = 0; 
					//}
				}
			} else {
				if (swarm_neighborhood_graph[who][particle]) {
					removeNeighbour(particle, who);
					//System.out.println("#################################");
				}
			}
			
			boolean ex_neighbours[];
			/*
			if (particle_failures_power[particle] > 10) {
				int particles_mutated = 0;
				double mutation_rate = random_number_generator.nextDouble();
				particle_failures_power[particle] = 1;
				System.out.println("tryin!");
				if (mutation_rate < 0.1) {
					System.out.println("BUM!");
					particle_failures_power[particle] = 1;
					particle_failures[particle] = 1;
					initializeParticle(particle);
					particles_mutated++;
					for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
						//if (swarm_neighborhood_graph[particle][neighbour] || swarm_neighborhood_graph[particle_best_neighbour[particle]][neighbour]) {
						if (swarm_neighborhood_graph[particle][neighbour]) {
							
//							for (int n_neighbour = 0; n_neighbour < NUMBER_OF_PARTICLES; n_neighbour++) {
//								particle_failures_power[particle] = 1;
//								initializeParticle(neighbour);
//								particles_mutated++;
//							}
							particle_failures_power[particle] = 1;
							penalizeNeigboursFailures(neighbour, 0);
							initializeParticle(neighbour);
							particles_mutated++;
						}
					}
					System.out.println("particles mutated: "+particles_mutated);
				}
				
			}
			*/
					
			if (particle_failures[particle] == 0) {
				boolean g_best_removed = removeNeighbourCarefully(particle, particle_best_neighbour[particle]);
				//removeNeighbourCarefully(particle, particle_best_neighbour[particle]);
				double mutation_rate = random_number_generator.nextDouble(); 
				
				if (mutation_rate < 0.30) {
				/*
					if (mutation_rate < 0.15) {
						removeAllNeigboursOfCarefullyBut(who, particle);
					} else {
						removeAllNeigboursOfCarefullyBut(particle, who);
					}
					penalizeNeigboursFailures(who, 0.5);
					particle_failures[who] = 0;
					*/
				}
								
				addNeighbour(who, particle);
				
				// remove todos os vizinhos da particula
				//ex_neighbours = removeAllNeigboursOf(particle);

				// para cada ex_vizinho
				/*
				int count = 0;
				for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
					boolean ex_n_neighbours[];
					if (count > 3) {
						break;
					}
					// ajusta para zero a quantidade de erros
					if (ex_neighbours[neighbour]) {
						particle_failures[neighbour] = 0;
						count ++;
					}
					// remove todos os vizinhos
					ex_n_neighbours = removeAllNeigboursOf(neighbour);
					// ajusta para zero a quantidade de erros dos vizinhos desse vizinho
					for (int n_neighbour = 0; n_neighbour < ex_n_neighbours.length; n_neighbour++) {
						if (ex_n_neighbours[n_neighbour]) {
							particle_failures[n_neighbour] = 0;
						} 
					}
				}
				*/
				setNeigboursFailuresTo(particle, 0);
				//setNeigboursFaliuresTo(who, 0); // TODO
				// e reconectar a nova vizinhanca
				for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
					if (swarm_neighborhood_graph[particle][neighbour]) {
						penalizeNeigboursFailures(neighbour, 0.5);
						for (int n_neighbour = 0; n_neighbour < NUMBER_OF_PARTICLES; n_neighbour++) {
							if (swarm_neighborhood_graph[n_neighbour][neighbour]) {
							//	penalizeNeigboursFailures(neighbour, 0.75);
							}
						}
						
						//setNeigboursFaliuresTo(neighbour, 0);
					}
					//if (ex_neighbours[neighbour]) {
//						addNeighbour(neighbour, particle);
	//				}
				}
			if (g_best_removed) {
				particle_failures[particle_best_neighbour[particle]] = 50;
			}

				
			}	
			/*
			/// LET's CONNECT TO TWO PARTICLES
			r = randomDouble();
			sum = 0;
			who = -1;
			for (wheel = (NUMBER_OF_PARTICLES - 1); wheel >= 0 ; wheel--) {
				// ESCOLHE PELA ROULETTE WHEEL UM PARA SE CONECTAR
				sum += ((float) wheel + 1) / (pa_sum);
				if (ranking[wheel] == particle) {
					continue;
				}
				if (sum >= r) {
					who = ranking[wheel];
					break;
				}
			}

			if ((ranking_of_particle[who] < ranking_of_particle[particle])) {
				if (!swarm_neighborhood_graph[who][particle]) {
					if (particle_failures[who] == 0) {
							particle_failures[particle] = 0;
							addNeighbour(who, particle);
					}
				}boolean[][] 
			} else {
				removeNeighbour(particle, who);
			}	
			*/
		}

	}
	
	
	private double exponencialGame(double x, double k) {
		return Math.exp(-k*x);
	}
	
	//exp(-2*$1/100
	
	private int chooseByRouletteWheelRankingBased_exponencial() {
		double r = randomDouble(); 
		int ranking[] = new int[NUMBER_OF_PARTICLES];
		int ranking_of_particle[] = new int[NUMBER_OF_PARTICLES];
		int who = 0;
		double sum = 0;
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			ranking[particle] = particle;
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
		// ranking[0] eh o index da PIOR particula
		// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
		// ranking_of_particle[i] eh a posicao da particula i no ranking
		r = randomDouble();
		sum = 0;
		who = NUMBER_OF_PARTICLES - 1;
		int wheel;
		for (wheel = (NUMBER_OF_PARTICLES); wheel > 0; wheel--) {
			// ESCOLHE PELA ROULETTE WHEEL UM PARA SE CONECTARexp(-2*$1/100)}
			//sum = ((float) exponencialGame(ranking[wheel]/(NUMBER_OF_PARTICLES/2),2));
			sum = Math.exp(-2*((double)wheel + 1)/100);
			
			if (sum >= r) {
				who = ranking[NUMBER_OF_PARTICLES - wheel];
				break;
			}			
		}

		return who;
	}

	
	private int chooseByRouletteWheelFitnessBased() {
		double r = randomDouble(); 
		int ranking[] = new int[NUMBER_OF_PARTICLES];
		int ranking_of_particle[] = new int[NUMBER_OF_PARTICLES];
		int who = 0;
		float sum = 0;
		float pa_sum = 0;
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			ranking[particle] = particle;
		}
		for (int particle = 0; particle < particle_best_value.length; particle++) {
			pa_sum += particle_best_value[particle];
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
		// ranking[0] eh o index da PIOR particula
		// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
		// ranking_of_particle[i] eh a posicao da particula i no ranking
		r = randomDouble();
		sum = 0;
		who = -1;
		int wheel;
		for (wheel = (NUMBER_OF_PARTICLES - 1); wheel >= 0 ; wheel--) {
			// ESCOLHE PELA ROULETTE WHEEL UM PARA SE CONECTAR
			sum += ((float) particle_best_value[ranking[wheel]]) / (pa_sum);
			
			if (sum >= r) {
				who = ranking[wheel];
				break;
			}
			
		}
		return who;
	}
	
	private boolean[] removeAllNeigboursOfCarefullyBut(int particle, int who) {
		boolean[] neighbours = new boolean[NUMBER_OF_PARTICLES];
		for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
			neighbours[neighbour] = swarm_neighborhood_graph[particle][neighbour];
			if (who == neighbour) 
				continue;
			if (swarm_neighborhood_graph[particle][neighbour]) {
				if (numberOfNeighbours(neighbour) > 1){
					removeNeighbour(neighbour, particle);
				}
			}
		}
		return neighbours;
	}
	
	private boolean[] removeAllNeigboursOf(int particle) {
		boolean[] neighbours = new boolean[NUMBER_OF_PARTICLES];
		for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
			neighbours[neighbour] = swarm_neighborhood_graph[particle][neighbour];
			if (swarm_neighborhood_graph[particle][neighbour]) {
				removeNeighbour(neighbour, particle);
			}
		}
		return neighbours;
	}
	
	
	
	private void penalizeNeigboursFailures(int particle, double value) {
		for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
			if (swarm_neighborhood_graph[particle][neighbour]) {
				particle_failures[neighbour] = (int) (((double) particle_failures[neighbour])*value);
			}
		}
	}
	
	private void setNeigboursFailuresTo(int particle, int value) {
		for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
			if (swarm_neighborhood_graph[particle][neighbour]) {
				particle_failures[neighbour] = 0;
			}
		}
	}
	
	private void updateNeighbourhood_BA_DYNAMIC_C() {
		int ranking[] = new int[NUMBER_OF_PARTICLES];
		int ranking_of_particle[] = new int[NUMBER_OF_PARTICLES];
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			ranking[particle] = particle;
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
		// ranking[0] eh o index da PIOR particula
		// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
		// ranking_of_particle[i] eh a posicao da particula i no ranking

		for (int particle = 0; particle < ranking_of_particle.length; particle++) {
			
			if (particle_failures[particle] < particles_failures_threshold) {
				continue;
			}
			
			double r = randomDoubleIndependent(); 
			int who = 0;
			float sum = 0;
			float pa_sum = (float) pa(NUMBER_OF_PARTICLES, 1);
			r = randomDoubleIndependent();
			sum = 0;
			who = -1;
			int n;
			//do {
				for (n = (NUMBER_OF_PARTICLES - 1); n >= 0 ; n--) {
					// ESCOLHE PELA ROULETTE WHEEL UM PARA SE CONECTAR
					sum += ((float) n + 1) / (pa_sum);
					if (ranking[n] == particle) {
						continue;
					}
					if (sum >= r) {
						if (numberOfNeighbours(ranking[n]) >= 20) {
						//	printf("ops, nao rolan neighbours: %d\n", numberOfNeighbours(ranking[n]));
							continue;
						}
						who = ranking[n];
						break;
					}
				}
			//} while (who == -1);


			// se vai conectar a um melhor, tudo bem.
			//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] == 0)) {
			//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] < failures[part])) {
			if (n > -1) {
				if ((ranking_of_particle[who] < ranking_of_particle[particle])) {
					if (particle_failures[who] == 0) {
						if (numberOfNeighbours(who) <= 20) {
							addNeighbour(particle, who);
							particle_failures[particle] = 0; 
						}
					}
					//	printf("connecting %d to: %d = ranking, who: %d, n = %d, pbest: %f\n", part, particles_ranking[who], who, n, pbest[who]);
				} else {
						removeNeighbour(particle, who);
				}
			}
			if (particle_failures[particle] == 0) {
				for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
					if (particle == neighbour)
						continue;
					if (swarm_neighborhood_graph[particle][neighbour]) {
						particle_failures[neighbour] /= 2;
					}
				}
			}
		}
	
	}
	
	private void updateNeighbourhood_BA_DYNAMIC_B() {
		int ranking[] = new int[NUMBER_OF_PARTICLES];
		int ranking_of_particle[] = new int[NUMBER_OF_PARTICLES];
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			ranking[particle] = particle;
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
		// ranking[0] eh o index da PIOR particula
		// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
		// ranking_of_particle[i] eh a posicao da particula i no ranking

		for (int particle = 0; particle < ranking_of_particle.length; particle++) {
			
			if (particle_failures[particle] < particles_failures_threshold) {
				continue;
			}
			
			double r = randomDoubleIndependent(); 
			int who = 0;
			float sum = 0;
			float pa_sum = (float) pa(NUMBER_OF_PARTICLES, 1);
			r = randomDoubleIndependent();
			sum = 0;
			who = -1;
			int n;
			//do {
				for (n = (NUMBER_OF_PARTICLES - 1); n >= 0 ; n--) {
					// ESCOLHE PELA ROULETTE WHEEL UM PARA SE CONECTAR
					sum += ((float) n + 1) / (pa_sum);
					if (ranking[n] == particle) {
						continue;
					}
					if (sum >= r) {
						if (numberOfNeighbours(ranking[n]) >= 20) {
						//	printf("ops, nao rolan neighbours: %d\n", numberOfNeighbours(ranking[n]));
							continue;
						}
						who = ranking[n];
						break;
					}
				}
			//} while (who == -1);


			// se vai conectar a um melhor, tudo bem.
			//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] == 0)) {
			//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] < failures[part])) {
			if (n > -1) {
				if ((ranking_of_particle[who] < ranking_of_particle[particle])) {
					if (particle_failures[who] == 0) {
						if (numberOfNeighbours(who) <= 20) {
							addNeighbour(particle, who);
							particle_failures[particle] = 0; 
						}
					}
					//	printf("connecting %d to: %d = ranking, who: %d, n = %d, pbest: %f\n", part, particles_ranking[who], who, n, pbest[who]);
				} else {
						removeNeighbour(particle, who);
				}
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

	private void updateNeighbourhood_BA_DYNAMIC_A() {
		/// CRIA O RANKING DAS PARTICULAS - INICIO
		int ranking[] = new int[NUMBER_OF_PARTICLES];
		int ranking_of_particle[] = new int[NUMBER_OF_PARTICLES];
		for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
			ranking[particle] = particle;
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
		// ranking[0] eh o index da PIOR particula
		// ranking[NUMBER_OF_AGENTS - 1] eh o index da MELHOR particula
		// ranking_of_particle[i] eh a posicao da particula i no ranking

		for (int a = 0; a < NUMBER_OF_PARTICLES; a++) {
			if (particle_failures[a] < particles_failures_threshold)
				continue;
	
			for (int particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
				for (int neighbour = 0; neighbour < NUMBER_OF_PARTICLES; neighbour++) {
					if (neighbour == particle)
						continue;
					double r = randomDoubleIndependent();
					int who = 0;
					double sum = 0;
					double pa_sum = (double) pa(NUMBER_OF_PARTICLES, 1);
					for (int n = (NUMBER_OF_PARTICLES - 1); n >= 0 ; n--) {
						// ESCOLHE PELA ROULETTE WHEEL UM PARA SE CONECTAR
						sum += ((float) n + 1) / (pa_sum);
						if (ranking[n] == particle) {
							continue;
						}
						if (sum >= r) {
							who = ranking[n];
							break;
						}
					}
					// se vai conectar a um melhor, tudo bem.
					//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] == 0)) {
					//if ((particles_ranking[who] > particles_ranking[part]) && (failures[who] < failures[part])) {
					if ((ranking_of_particle[who] < ranking_of_particle[particle]) && (particle_failures[who] < particle_failures[particle]) ) {
						if ((who == neighbour)) {
							addNeighbour(particle, who);
						}
					} else {
						removeNeighbour(particle, who);
					}
				}
				removeNeighbour(particle, particle_best_neighbour[particle]); // modified
				particle_failures[a] = 0;
			}
		}
		
	}

	public int getGBest() {
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
			double[] cos_test = new double[DIMENSION];
			if (particle == 0) {
				cos_test = particle_velocity[0].clone();
			}
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
			if (particle == 0) {
				//System.out.println("cos 0 = " + this.cosBetweenVectors(cos_test, particle_velocity[0]));
			}
		}
	}
	
	
	
	
	

	private double randomDouble() {
		return this.random_number_generator.nextDouble();
	}
	
	private double randomDoubleIndependent() {
		return this.random_number_generator_independent.nextDouble();
	}
	
	private void findBestNeighbours() {
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
		double zero = 0.000000000000000000000000000001;
		// evaluate function and update pbest
		for (int particle = 0; particle < particle_position.length; particle++) {
			
			double current_fitness = FUNCTION.compute(particle_position[particle]);
			if (current_fitness < particle_best_value[particle]) {
				double delta = Math.abs(current_fitness -  particle_best_value[particle]);
				if (delta < 0.00001D) {
					///System.out.println("DELTA EH PEQUENO!");
					particle_failures[particle]++;
				} else {
					particle_failures[particle] = 0;
					particle_best_value[particle] = current_fitness;
					for (int dimension = 0; dimension < particle_best_value_position[particle].length; dimension++) {
						particle_best_value_position[particle][dimension] = particle_position[particle][dimension];
					}
				} 
			} else {
				particle_failures[particle]++;
			}
			//}
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
			// TODO: need tests here. Checked: I think it is OK.
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
		case K_MEANS_BASED:
			k_meanNeighbourhoodCreate(swarm_number_of_clusters);
			break;
		case KOHONEN_BASED:
			kohonenNeighbourhoodCreate(swarm_number_of_clusters);
			break;
		default:
			break;
		}	
	}

	private void kohonenNeighbourhoodCreate(int k) {
		swarm_clusters = null;
		swarm_centroids = null;
		swarm_clusters_kohonen = null;
		int[] clusters = this.kohonen(k);
		int[] counter = new int[clusters.length];
		int[] last_connected = new int[clusters.length];
		int[] first_connected = new int[clusters.length];
		/*
		for (int i = 0; i < k; i++) {
			counter[i] = 0;
		}
		for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
			counter[clusters[i]]++;
		}
		for (int i = 0; i < k; i++) {
			System.out.println(counter[i]);
		}*/
		for (int i = 0; i < swarm_neighborhood_graph.length; i++) 	{
			for (int j = i+1; j < swarm_neighborhood_graph[i].length; j++) {
				removeNeighbour(i, j);
				if (i == j) {
					continue;
				}
				
				if (clusters[i] == clusters[j]) {
					// GLOBAL
					/* 
					addNeighbour(i, j);
					*/
					/*
					// LOCAL
					if (first_connected[clusters[i]] == 0) {
						first_connected[clusters[i]] = i;
					}
					addNeighbour(i, j);
					last_connected[clusters[j]] = j;
					break;
					*/
					
					
					// global but limited 
					if (numberOfNeighbours(j) < swarm_initial_maximum_neighbors) {
						addNeighbour(i, j);
					}
				}
			}
		}
		/*
		// LOCAL
		for (int cluster = 0; cluster < clusters.length; cluster++) {
			addNeighbour(first_connected[cluster], last_connected[cluster]);
		}
		*/
		//swarm_centroids = null;
		//swarm_clusters = null;
		
	}
	private void k_meanNeighbourhoodCreate(int k) {
		int[] clusters = this.k_means(k);
		int[] last_connected = new int[clusters.length];
		int[] first_connected = new int[clusters.length];
		for (int i = 0; i < swarm_neighborhood_graph.length; i++) 	{
			for (int j = i+1; j < swarm_neighborhood_graph[i].length; j++) {
				removeNeighbour(i, j);
				if (i == j) {
					continue;
				}
				if (clusters[i] == clusters[j]) {
					/*
					if (numberOfNeighbours(j) < swarm_initial_maximum_neighbors) {
						addNeighbour(i, j);
					}
					*/
				}
			}
		}
		swarm_centroids = null;
		swarm_clusters = null;
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


