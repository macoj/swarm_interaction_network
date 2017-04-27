__author__ = 'marcos'
import igraph


class InfluenceGraph():
    def __init__(self):
        pass

    @staticmethod
    def create_graph_from_matrix(graph_matrix):
        igraph_graph = igraph.Graph.Weighted_Adjacency(graph_matrix.tolist(), mode=igraph.ADJ_MAX)
        return igraph_graph