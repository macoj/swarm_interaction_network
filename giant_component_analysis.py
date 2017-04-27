__author__ = 'marcos'
import pandas as pd
from influence_graph import InfluenceGraph


class GiantComponentDeath:
    def __init__(self):
        pass

    @staticmethod
    def low_edges_weight_removal(
            igraph_graph, graphs_keeper_giant_size_threshold=None, include_zero=True, count='components'):
        """ Removes the edges progressively by the weight starting with the ones with the lowest weight and returns
        the size of the giant component after each removal.

        :param igraph_graph:
        :param graphs_keeper_giant_size_threshold:
        :return:
        """
        assert count == "components" or count == "size", \
            "I can only count the number of components or giant component size!"
        total_nodes = float(igraph_graph.vcount())
        # it can keep graph structures for certain thresholds, thus:
        graphs = []
        graphs_kept = []
        if graphs_keeper_giant_size_threshold:
            graphs_kept = list(set(graphs_keeper_giant_size_threshold))
            graphs_kept.sort(reverse=True)
            graphs_kept = map(lambda x: (False, x), graphs_kept)
        death_evolution_weight = []
        death_evolution_size = []
        number_of_components = []
        weight_values = list(set(igraph_graph.es['weight']))
        if include_zero:
            weight_values.append(0)
        weight_values.sort()
        graph_copy = igraph_graph.copy()
        for weight in weight_values:
            GiantComponentDeath.remove_weighted_edges(graph_copy, weight)
            death_evolution_weight.append(weight)
            size_perc = graph_copy.components().giant().vcount()/total_nodes
            death_evolution_size.append(size_perc)
            number_of_components.append(len(graph_copy.components()))
            ### this will keep the graphs
            if graphs_keeper_giant_size_threshold:
                index_graph = 0
                for graph in graphs_kept:
                    if not graph[0]:
                        if graph[1] >= size_perc:
                            graphs.append((graph[1], size_perc, graph_copy.copy()))
                            graphs_kept[index_graph] = (True, graph[1])
                    index_graph += 1
        # the result dataframe is sorted by the weight, low -> high
        if count == 'components':
            pd_df = pd.DataFrame({'x': death_evolution_weight, 'y': number_of_components})
        else:
            pd_df = pd.DataFrame({'x': death_evolution_weight, 'y': death_evolution_size})
        if graphs is not []:
            result = (pd_df, graphs)
        else:
            result = pd_df
        return result

    @staticmethod
    def remove_weighted_edges(igraph_graph, threshold=0.0):
        remove_them = []
        for e in igraph_graph.es:
            if e['weight'] <= threshold:
                remove_them.append(e.index)
        igraph_graph.delete_edges(remove_them)

    @staticmethod
    def create_giant_component_curve(
            graph_matrix, return_graphs_with_giant_sizes=None, normalize=None, adjusted=False, **kargs):
        igraph_graph = InfluenceGraph.create_graph_from_matrix(graph_matrix)
        # create the graph objects as well as the death analysis
        pd_data, graphs = GiantComponentDeath.low_edges_weight_removal(
            igraph_graph, return_graphs_with_giant_sizes, **kargs)
        # pd_data, graphs = GiantComponentDeath.nodes_degree_removal(igraph_graph, return_graphs_with_giant_sizes)
        # the weights leading to the destruction of the graph can be normalized.
        # when there is a time window tw, the maximum weight of an edge is equal
        # to 2*tw, this is the case of two particles sharing information all the
        # time.
        if normalize:
            pd_data['x'] /= normalize
        if adjusted:
            pd_data['x'] -= min(pd_data.x)  # to have zero (we admit that the min value is positive)
        return pd_data

    @staticmethod
    def create_giant_component_curves(graph_matrices, adjusted=False, weight_normalize=None, **kargs):
        pd_datas = []
        normalize_index = 0
        if type(graph_matrices) == dict:
            graph_matrices = graph_matrices.values()
        for graph_matrix in graph_matrices:
            if type(graph_matrix) == tuple:
                _, graph_matrix = graph_matrix
            if type(weight_normalize) == list:
                normalize_c = weight_normalize[normalize_index]
                normalize_index += 1
            else:
                normalize_c = weight_normalize
            pd_data = GiantComponentDeath.create_giant_component_curve(
                graph_matrix, normalize=normalize_c, adjusted=adjusted, **kargs)
            pd_datas.append(pd_data)
        return pd_datas

