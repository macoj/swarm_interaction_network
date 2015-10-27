from giant_component_analysis_plotter import GiantComponentDeathPlotter

__author__ = 'marcos'
from swarm_parser import SwarmParser
from giant_component_analysis import GiantComponentDeath
from pre_callbacks import PreCallback


class SwarmAnalyzer:
    def __init__(self):
        pass

    @staticmethod
    def read_files_and_plot(filenames, windows_size=None, calculate_on=None):
        # windows_size = [10, 40, 50, 60, 70, 80, 90, 100]
        # windows_size = [1, 5, 10, 15, 20] #, 40, 50, 60, 70, 80, 90, 100]
        # windows_size = [300, 400, 500]
        # calculates_on = [1300, 1400, 1500]
        # calculate_on = 1500
        fitness_grep = 'it\:#[0-9]*'
        influence_graph_grep = 'ig\:#[0-9]*'
        pre_callback = PreCallback.to_symmetric
        all_graph_matrices = {}
        for filename in filenames:
            title, filename = filename
            graph_matrices = []
            #for calculate_on in calculates_on:
            if type(windows_size) == int:
                windows_size = [windows_size]
            for window_size in windows_size:
                graph_index = window_size
                pos_callback = lambda x, y: graph_matrices.append((graph_index, x[1]))
                SwarmParser.read_file_and_measure(filename,
                                                  calculate=None,
                                                  influence_graph_grep=influence_graph_grep,
                                                  fitness_grep=fitness_grep,
                                                  window_size=window_size,
                                                  pre_callback=pre_callback,
                                                  pos_callback=pos_callback,
                                                  calculate_on=calculate_on)
            all_graph_matrices[title] = graph_matrices
            ### create the GiantComponentDeath analysis
        pd_datas = GiantComponentDeath.create_giant_component_curves(all_graph_matrices)
        GiantComponentDeathPlotter.create_giant_component_death_curve(calculate_on, pd_datas, windows_size)
        #create_strength_distribution_curves_windows_comparison(all_graph_matrices, calculate_on, windows_size)
        #create_heatmap_plot(all_graph_matrices, calculate_on)
        #create_strength_distribution_curves(all_graph_matrices, calculate_on)

        # pd_data = (title, pd_data)
            # pd_datas_2.append(pd_data)

            # # but 'graphs' is actually igraph.graph, but we need
            # # networkx graphs, dammit! (because just nx.graph can be plot with matplotlib :( -- it seems)
            # nx_graphs = []
            # graph = None
            # for graph in graph_matrix:
            #     graph_component_histogram = graph[2].components().sizes()
            #     nx_graph = from_igraph_to_nxgraph(graph[2], only_connected_nodes=True)
            #     title = str(graph[1]) + " ("+str(graph[0])+") [" + str(nx.number_of_nodes(nx_graph)) \
            #                           + "/" + str(graph[2].vcount()) + "]"
            #     nx_graphs.append((title, nx_graph, graph_component_histogram))
            # if not nx_graphs:
            #     nx_graphs = None
            #
            # ### here is the fitness data
            # pd_data_1 = None
            # if fitness is not None:
            #     pd_data_1 = pd.DataFrame({'x': range(len(fitness)), 'y': fitness})
            #     pd_data_1 = ('Fitness', pd_data_1)

            ### create the histograms data
            # gets the last graph in 'graphs' and plot the degree distribution of it
      #return graph_matrix

    @staticmethod
    def do_it():
        windows_size = [10, 100, 500, 1000]
        # return SwarmAnalyzer.read_files_and_plot([
        #     ('FSS1',
        #      '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original'),
        #     ('FSS2',
        #      '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original'),
        #     ('FSS3',
        #      '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original'),
        #     ('FSS4',
        #      '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original')],
        #     windows_size=-1,
        #     calculate_on=1000)/mnt/50_particles_simulations
        return SwarmAnalyzer.read_files_and_plot([
            ('Dynamic',
             '/mnt/50_particles_simulations/pso_dynamic_initial_ring_F6_16'),
            ('Ring',
             '/mnt/50_particles_simulations/pso_ring_F6_13'),
            ('Global',
             '/mnt/50_particles_simulations/pso_global_F6_16'),
            ('von Neumann',
             '/mnt/50_particles_simulations/pso_dynamic_initial_ring_F6_16')],
            windows_size=[10, 500, 1000],
            calculate_on=1000)
        # return

"""
execfile("swarm_analyzer.py")
SwarmAnalyzer.do_it()
"""