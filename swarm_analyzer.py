__author__ = 'marcos'
import igraph
import pandas as pd
from opt.callbacks import Callback
from swarm_parser import SwarmParser
from influence_graph import InfluenceGraph
from giant_component_analysis import GiantComponentDeath


class SwarmAnalyzer:
    def __init__(self):
        pass

    @staticmethod
    def create_influence_graph(filename, influence_graph_grep='ig\:#', window_size=1000, calculate_on=1000):
        pre_callback = Callback.to_symmetric
        graph, _ = SwarmParser.read_file_and_measures(
            filename, influence_graph_grep=influence_graph_grep, window_size=window_size, pre_callback=pre_callback,
            calculate_on=calculate_on)
        igraph_graph = InfluenceGraph.create_graph_from_matrix(graph[0][1])
        return igraph_graph

    @staticmethod
    def create_influence_graph_graphml(filename, output_file_name, window_size=1000, calculate_on=1000):
        igraph_graph = SwarmAnalyzer.create_influence_graph(
            filename, window_size=window_size, calculate_on=calculate_on)
        InfluenceGraph.to_graphml(igraph_graph, output_file_name)

    @staticmethod
    def get_giant_component_destruction_curves(filename, window_size, until=-1, calculate_on=-1, count='components'):
        filenames = [('', filename)]
        all_graph_matrices = SwarmAnalyzer.get_graph_matrices_from_files(
            filenames, windows_size=[window_size], until=until, calculate_on=calculate_on)
        graph_matrices = all_graph_matrices[''][window_size]
        # the maximum is 2*tw, where tw is the window size, but if the graph is from an iteration t that is less than
        # tw, then the maximum is 2*t, therefore:
        weight_normalize = [2.0 * i if i < window_size else 2.0 * window_size for i in range(len(graph_matrices))]
        curves = GiantComponentDeath.create_giant_component_curves(
            graph_matrices, adjusted=True, include_zero=False, weight_normalize=weight_normalize, count=count)
        return curves
    """
    execfile("swarm_analyzer.py")
    filename = './data/vonneumann_F06_15'
    window_size = 10
    until = 100
    curves = SwarmAnalyzer.get_giant_component_destruction_curves(filename, window_size, calculate_on=100, count='count')
    import matplotlib.pyplot as plt
    plt.plot(curves[0]['x'], curves[0]['y'])
    plt.show()
    """

    @staticmethod
    def get_graph_matrices_from_files(filenames, influence_graph_grep='ig\:#', **kargs):
        pos_callback = Callback.to_symmetric
        graph_matrices, _ = SwarmParser.read_files_and_measures(
            filenames, influence_graph_grep=influence_graph_grep, pos_callback=pos_callback, **kargs)
        return graph_matrices
    """
    execfile("swarm_analyzer.py")
    filenames = [('Global', "./data/global_F06_15"), ('Ring', "./data/ring_F06_15")]
    df = SwarmAnalyzer.get_graph_matrices_from_files(filenames, windows_size=[100, 1000], calculate_on=1000)
    df['Global'][100][0]
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
    def get_number_of_components(filename, window_size, min_weight, influence_graph_grep='ig\:#', **kargs):
        pos_callback = lambda x: SwarmAnalyzer.get_number_of_components_of_graph(
            x, min_weight=min_weight * 2 * window_size, pre_callback=Callback.to_symmetric)
        all_graph_matrices, _ = SwarmParser.read_file_and_measures(
            filename, influence_graph_grep=influence_graph_grep, pos_callback=pos_callback,
            window_size=window_size, **kargs)
        return all_graph_matrices

    @staticmethod
    def get_swarm_informations_from_file(filename, informations_grep, information_map=float, **kargs):
        _, informations = SwarmParser.read_files_and_measures(
            [('', filename)], informations_grep=informations_grep, information_map=information_map, **kargs)
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
