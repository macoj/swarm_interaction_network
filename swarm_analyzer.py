__author__ = 'marcos'
import sys
import time
import igraph
import copy
import pandas as pd
import numpy as np
from plotter import Plotter
from scipy import interpolate
from callbacks import Callback
from swarm_parser import SwarmParser
from giant_component_analysis import GiantComponentDeath
from giant_component_analysis_plotter import GiantComponentDeathPlotter
from scipy.stats import pearsonr
import scipy.stats as st
from statsmodels.distributions.empirical_distribution import ECDF
from functools import reduce
from operator import add


class SwarmAnalyzer:
    def __init__(self):
        pass

    @staticmethod
    def get_giant_component_destruction_curves(filename, window_size, until=-1):
        influence_graph_grep = 'ig\:#'
        pos_callback = Callback.to_symmetric
        all_graph_matrices, _ = SwarmParser.read_files_and_measures([('', filename)],
                                                                    influence_graph_grep=influence_graph_grep,
                                                                    pos_callback=pos_callback,
                                                                    windows_size=[window_size], until=until)
        graph_matrices = all_graph_matrices[''][window_size]
        # the maximum is 2*tw, where tw is the window size, but if the graph is from an iteration t that is less than
        # tw, then the maximum is 2*t, therefore:
        weight_normalize = [2.0 * i if i < window_size else 2.0 * window_size for i in range(len(graph_matrices))]
        curves = GiantComponentDeath.create_giant_component_curves(graph_matrices, adjusted=True, include_zero=False,
                                                                   weight_normalize=weight_normalize)
        return curves

    @staticmethod
    def get_areas_under_curves(curves, delta=0.001, normalize=False, normalize_max_y=None, normalize_constant=None):
        areas = []
        tx = np.arange(0, 1 + delta, delta)
        normalization = None
        if normalize:
            assert normalize_max_y is not None or normalize_constant is not None, "[ERROR] Missing normalizing factor!"
            normalization = float(normalize_max_y) * len(tx)
        for graph in curves:
            x, y = list(graph.x), list(graph.y)
            x.append(1.0)
            y.append(y[len(y) - 1])
            f = interpolate.interp1d(x, y, kind='nearest')
            ty = map(float, map(f, tx))
            # areas.append(sum(ty * tx)/len(tx))
            areas.append(sum(ty))
            del x
            del y
            del f
            del ty
        if normalization:
            areas = [a / normalization for a in areas]
        return areas

    @staticmethod
    def get_giant_component_destruction_area(filename, window_size, number_of_individuals=100, until=-1):
        graphs = SwarmAnalyzer.get_giant_component_destruction_curves(filename, window_size=window_size, until=until)
        areas = SwarmAnalyzer.get_areas_under_curves(graphs, normalize=True, normalize_max_y=number_of_individuals)
        # df = pd.DataFrame({'x': range(window_size, window_size + len(areas)), 'y': areas})
        df = pd.DataFrame({'x': range(len(areas)), 'y': areas})
        return df

    """
    execfile("swarm_analyzer.py")
    filename = "/mnt/pso_100_particles/global_F06_00"
    df = SwarmAnalyzer.get_giant_component_destruction_area(filename, 100)
    """

    @staticmethod
    def get_number_of_components_of_graph(graph, min_weight=None, pre_callback=None):
        if pre_callback:
            graph = pre_callback(graph)
        igraph_graph = igraph.Graph.Weighted_Adjacency(graph.tolist(), mode=igraph.ADJ_MAX)
        if min_weight is not None:
            GiantComponentDeath.remove_weighted_edges(igraph_graph, min_weight)
        components = len(igraph_graph.components())
        return components

    @staticmethod
    def get_number_of_components(filename, window_size, min_weight, **kargs):
        influence_graph_grep = 'ig\:#'
        pos_callback = lambda x: SwarmAnalyzer.get_number_of_components_of_graph(x,
                                                                                 min_weight=min_weight * 2 * window_size,
                                                                                 pre_callback=Callback.to_symmetric)
        all_graph_matrices, _ = SwarmParser.read_file_and_measures(filename,
                                                                   influence_graph_grep=influence_graph_grep,
                                                                   pos_callback=pos_callback,
                                                                   window_size=window_size, **kargs)
        return all_graph_matrices

    @staticmethod
    def create_influence_graph_graphml(filename, output_file_name, window_size=1000, calculate_on=1000):
        influence_graph_grep = 'ig\:#'
        pre_callback = Callback.to_symmetric
        # for calculate_on in calculates_on:
        graph, _ = SwarmParser.read_file_and_measures(filename, influence_graph_grep=influence_graph_grep,
                                                      window_size=window_size, pre_callback=pre_callback,
                                                      calculate_on=calculate_on)
        igraph_graph = igraph.Graph.Weighted_Adjacency(graph[0][1].tolist(), mode=igraph.ADJ_MAX)
        igraph.Graph.write_graphml(igraph_graph, output_file_name)
    """
    filename = './data/100_particles/vonneumann_F06_15.teste'
    output_file_name = './vonneumann.graphml'
    """

    @staticmethod
    def get_swarm_informations_from_file(filename, informations_grep, information_map=float, **kargs):
        _, informations = SwarmParser.read_files_and_measures([('', filename)], informations_grep=informations_grep,
                                                              information_map=information_map, **kargs)
        informations = informations[''][-1]  # there is no window here!
        # let's get the longest information sequence
        max_information_length = float("-inf")
        max_information_key = None
        for information_grep in informations:
            if len(informations[information_grep]) > max_information_length:
                max_information_length = len(informations[information_grep])
                max_information_key = information_grep
        iterations = [iteration for (iteration, _) in informations[max_information_key]]
        iterations.sort()
        df = pd.DataFrame({'x': iterations})
        for information_grep in informations:
            dict_information = dict(informations[information_grep])
            df[information_grep] = [dict_information[i] if i in dict_information else float("nan") for i in iterations]
        return df

    @staticmethod
    def read_hdfs_and_plot(basename):
        filenames = [basename + "%02d" % i for i in range(1, 30)]
        windows_size = 1000
        for filename in filenames:
            df = pd.read_hdf(filename + "_" + str(windows_size) + ".hdf", 'df')
            Plotter.plot_curve(df, figsize=(18, 6))

    @staticmethod
    def read_files_and_plot(filenames, windows_size, calculate_on):
        influence_graph_grep = 'ig\:#'
        pre_callback = Callback.to_symmetric
        graph_matrices, _ = SwarmParser.read_files_and_measures(
            filenames, influence_graph_grep=influence_graph_grep, pos_callback=pre_callback, windows_size=windows_size,
            calculate_on=calculate_on)
        normalize = [2 * i for i in windows_size]
        pd_datas = []
        for title, _ in filenames:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(graphs, weight_normalize=normalize)
            pd_datas.append((title, dict(zip(windows_size, curves_areas))))
        GiantComponentDeathPlotter.giant_component_death_curve(
            calculate_on, pd_datas, windows_size, xlim=(0, 1.0), figsize=(4.5, 4))
    """
    execfile("swarm_analyzer.py")
    execfile("giant_component_analysis_plotter.py")
    filenames = [('Global', "./data/global_F06_15"), ('Ring', "./data/ring_F06_15"), ('Von Neumann', "./data/vonneumann_F06_15"), ('Dynamic', "./data/dynamicring_F06_15")]
    df = SwarmAnalyzer.read_files_and_plot(filenames, windows_size=[100, 1000], calculate_on=1000)
    """
