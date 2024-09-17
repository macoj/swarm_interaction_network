#! /usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'marcos'
import igraph


class InfluenceGraph():
    def __init__(self):
        pass

    @staticmethod
    def create_graph_from_matrix(graph_matrix):
        igraph_graph = igraph.Graph.Weighted_Adjacency(graph_matrix.tolist(), mode=igraph.ADJ_MAX)
        return igraph_graph

    @staticmethod
    def to_graphml(igraph_graph, output_file_name):
        igraph.Graph.write_graphml(igraph_graph, output_file_name)