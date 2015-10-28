import sys
from giant_component_analysis_plotter import GiantComponentDeathPlotter
__author__ = 'marcos'
from swarm_parser import SwarmParser
from giant_component_analysis import GiantComponentDeath
from pre_callbacks import PreCallback
from plotter import Plotter
import pandas as pd
from scipy import interpolate
import numpy as np


class SwarmAnalyzer:
    def __init__(self):
        pass

    @staticmethod
    def read_file_and_plot(filename, windows_size=None, calculate_on=None):
        fitness_grep = 'it\:#[0-9]*'
        influence_graph_grep = 'ig\:#[0-9]*'
        pre_callback = PreCallback.to_symmetric
        if type(windows_size) == int:
            windows_size = [windows_size]
        all_graph_matrices = SwarmParser.read_files_and_measure(calculate_on, [filename], fitness_grep,
                                                                influence_graph_grep, pre_callback, windows_size)
        title, filename = filename
        graph_matrices = all_graph_matrices[title]
        curves = GiantComponentDeath.get_giant_component_curves_areas(graph_matrices)
        return curves

    @staticmethod
    def read_files_and_export_hdf(windows_size, basename):
        filenames = [basename+"%02d" % i for i in range(1, 30)]
        for filename in filenames:
            graphs = SwarmAnalyzer.read_file_and_plot(('None', filename), windows_size=windows_size, calculate_on=-1)
            graphs = graphs[:10]
            areas = []
            delta = 0.001
            tx = np.arange(0, 1 + delta, delta)
            for graph in graphs:
                x, y = list(graph.x), list(graph.y)
                x.append(1.0)
                y.append(y[len(y)-1])
                f = interpolate.interp1d(x, y, kind='nearest')
                ty = map(float, map(f, tx))
                areas.append(sum(ty * tx)/len(tx))
                del x
                del y
                del f
                del ty
            df = pd.DataFrame({'x': range(windows_size, windows_size+len(areas)), 'y': areas})
            df.to_hdf(filename+"_"+str(windows_size)+".hdf", 'df')
            del df
    """
execfile("swarm_analyzer.py")
basename = "/mnt/50_particles_simulations/pso_global_F6_"
basename = "/mnt/50_particles_simulations/pso_ring_F6_"
basename = "/mnt/50_particles_simulations/pso_dynamic_initial_ring_F6_"
SwarmAnalyzer.read_files_and_export_hdf(100, basename)
    """
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
        if type(windows_size) == int:
            windows_size = [windows_size]
        all_graph_matrices = SwarmParser.read_files_and_measure(calculate_on, filenames, fitness_grep,
                                                                influence_graph_grep, pre_callback, windows_size)
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

if __name__ == "__main__":
    SwarmAnalyzer.read_files_and_export_hdf(int(sys.argv[1]), sys.argv[2])

"""
execfile("swarm_analyzer.py")
SwarmAnalyzer.do_it()
"""