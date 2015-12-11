__author__ = 'marcos'
import pandas as pd
import igraph
from giant_component_analysis_plotter import GiantComponentDeathPlotter


class GiantComponentDeath:
    def __init__(self):
        pass

    @staticmethod
    def nodes_degree_removal(igraph_graph, graphs_keeper_giant_size_threshold=None, start_with_highest=True):
        #total_nodes = float(igraph_graph.vcount())
        graphs = []
        graphs_kept = []
        if graphs_keeper_giant_size_threshold:
            graphs_kept = list(set(graphs_keeper_giant_size_threshold))
            graphs_kept.sort(reverse=True)
            graphs_kept = map(lambda x: (False, x), graphs_kept)

        death_evolution_degree = []
        death_evolution_size = []
        graph_copy = igraph_graph.copy()
        removals = 0
        while graph_copy.vcount() > 0:
            removals += 1
            # highest_degree = GiantComponentDeath.remove_nodes_with_highest_degree(graph_copy)
            total_nodes = float(graph_copy.vcount())
            print str(total_nodes)
            if total_nodes > 0:
                size_perc = graph_copy.components().giant().vcount()/total_nodes
            else:
                size_perc = 0
            death_evolution_size.append(size_perc)
            death_evolution_degree.append(removals)
            if graphs_keeper_giant_size_threshold:
                index_graph = 0
                for graph in graphs_kept:
                    if not graph[0]:
                        if graph[1] >= size_perc:
                            graphs.append((graph[1], size_perc, graph_copy.copy()))
                            graphs_kept[index_graph] = (True, graph[1])
                    index_graph += 1
        pd_df = pd.DataFrame({'x': death_evolution_degree, 'y': death_evolution_size})
        if graphs is not []:
            result = (pd_df, graphs)
        else:
            result = pd_df
        return result

    @staticmethod
    def low_edges_weight_removal(igraph_graph, graphs_keeper_giant_size_threshold=None, include_zero=True):
        """ Removes the edges progressively by the weight starting with the ones with the lowest weight and returns
        the size of the giant component after each removal.

        :param igraph_graph:
        :param graphs_keeper_giant_size_threshold:
        :return:
        """
        total_nodes = float(igraph_graph.vcount())
        # it can keep graph structures for certain threhsolds, thus:
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
        #pd_df = pd.DataFrame({'x': death_evolution_weight, 'y': death_evolution_size})
        pd_df = pd.DataFrame({'x': death_evolution_weight, 'y': number_of_components})
        #print str(number_of_components)
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
        #return (remove_them)

    @staticmethod
    def remove_nodes_with_highest_degree(igraph_graph):
        highest_degree = max(igraph_graph.strength(weights='weight'))
        GiantComponentDeath.remove_nodes_with_degree(igraph_graph, highest_degree)
        return highest_degree

    @staticmethod
    def remove_nodes_with_degree(igraph_graph, degree):
        remove_them = []
        if 'weight' in igraph_graph.edge_attributes():
            degrees = igraph_graph.strength(weights='weight')
        else:
            degrees = igraph_graph.vs.degree()
        for v in igraph_graph.vs.indices:
            if degrees[v] == degree:
                remove_them.append(v)
        a = igraph_graph.vcount()
        igraph_graph.delete_vertices(remove_them)
        a -= igraph_graph.vcount()
        print 'removing all with degree ' + str(degree) + ' - removed ' + str(a)

    @staticmethod
    def create_giant_component_curve(graph_matrix, return_graphs_with_giant_sizes=None,
                                     normalize=None, adjusted=False, include_zero=True):
        # igraph_graph = igraph.Graph.Weighted_Adjacency(graph_matrix.tolist(), mode=igraph.ADJ_MAX)
        igraph_graph = igraph.Graph.Weighted_Adjacency(graph_matrix.tolist(), mode=igraph.ADJ_PLUS)
        # create the graph objects as well as the death analysis
        pd_data, graphs = GiantComponentDeath.low_edges_weight_removal(igraph_graph,
                                                                       return_graphs_with_giant_sizes,
                                                                       include_zero=include_zero)
        # pd_data, graphs = GiantComponentDeath.nodes_degree_removal(igraph_graph, return_graphs_with_giant_sizes)
        # the weights leading to the destruction of the graph can be normalized.
        # when there is a time window tw, the maximum weight of an edge is equal
        # to 2*tw, this is the case of two particles sharing information all the
        # time.
        if normalize:
            pd_data['x'] /= normalize
        if adjusted:
            pd_data['x'] -= min(pd_data.x)
        return pd_data

    @staticmethod
    def create_giant_component_curves(graph_matrices, adjusted=False, include_zero=True, weight_normalize=None):
        pd_datas = []
        normalize_index = 0
        if type(graph_matrices) == dict:
            graph_matrices = graph_matrices.values()
        for graph_matrix in graph_matrices:
            print graph_matrix
            if type(graph_matrix) == tuple:
                _, graph_matrix = graph_matrix
                print graph_matrix
            if type(weight_normalize) == list:
                normalize_c = weight_normalize[normalize_index]
                normalize_index += 1
            else:
                normalize_c = weight_normalize
            pd_data = GiantComponentDeath.create_giant_component_curve(graph_matrix,
                                                                       normalize=normalize_c,
                                                                       adjusted=adjusted,
                                                                       include_zero=include_zero)
            pd_datas.append(pd_data)
        return pd_datas

