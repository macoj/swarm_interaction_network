import re
import igraph
import matplotlib
import sys
import numpy
import scipy
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import chain, cycle
from scipy.spatial.distance import pdist, squareform, euclidean
from scipy.cluster.hierarchy import dendrogram
from abc import abstractmethod, ABCMeta
from fastcluster import *
from matplotlib.lines import Line2D
import statsmodels.api as sm
import powerlaw


class Measure:
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def calculate(self, matrix):
        pass


class Spectra(Measure):
    @staticmethod
    def calculate(self, matrix):
        #ev, = numpy.linalg.eig(matrix)
        return numpy.linalg.eig(matrix)


class SpectraUndirect(Measure):
    @staticmethod
    def calculate(self, matrix):
        #ev, = numpy.linalg.eig(matrix)
        size = len(matrix)
        for i in range(size):
            for j in range(size):
                if matrix[i][j] == matrix[j][i]:
                    continue
                if matrix[i][j] == 0 or matrix[j][i] == 0:
                    weight = max(matrix[i][j], matrix[j][i])
                    matrix[i][j] = weight
                    matrix[j][i] = weight
                else:
                    # the weights are different and nonzero, 
                    # many approaches can be used here: min,
                    # avg, max, etc. Let's do it 'max'. 
                    weight = max(matrix[i][j], matrix[j][i])
                    matrix[i][j] = weight
                    matrix[j][i] = weight
        teste = 0
        for i in range(size):
            for j in range(size):
                if matrix[i][j] != matrix[j][i]:
                    print "############################################################################## ERRR"
                    teste = 1
        if teste == 0:
            print "############################################################################## OK"
        print (str(matrix))
        return numpy.linalg.eig(matrix)






def create_igraph_from_matrix(matrix, mode=igraph.ADJ_UNDIRECTED):
    if mode == igraph.ADJ_UNDIRECTED:
        # we are going to sum i,j to j,i:
        matrix = to_symmetric(matrix)
        # so, in the swarm, it means how many times a particle i and 
        # a particle j shared information (i->j and j->i are two information exchange) 
    g = igraph.Graph.Weighted_Adjacency(matrix.tolist())
    return g
    
        



def create_distance_matrix_from_positions(particle_positions, distance=euclidean):
    distance_matrix = None
    if particle_positions is not None:
        number_of_particles = len(particle_positions)
        distance_matrix = numpy.zeros((number_of_particles, number_of_particles))
        for i in particle_positions:
            for j in particle_positions:
                #todo: please, you should improve this...
                particle_distance = distance(particle_positions[i], particle_positions[j])
                distance_matrix[int(i)][int(j)] = particle_distance
    return distance_matrix


def remove_diagonal(array):
    array_result = array
    if array.shape[0] == array.shape[1]:
        array_reshaped = array.reshape((1, array.size))[0]
        diagonal_elements = range(0, array.size, array.shape[0]+1)
        array_removed = numpy.delete(array_reshaped, diagonal_elements)
        array_result = array_removed  # .reshape((array.shape[0]-1, array.shape[0]))
    return array_result


def calculate_correlation_influence_position(influence_graph, particles_positions):
    #todo: we may need to add the idea of transitivity, but we will need
    #todo: to explain why we do not take into account in the analyses
    distance_matrix = create_distance_matrix_from_positions(particles_positions)
    distance_matrix /= distance_matrix.max()
    distance_matrix = 1 - distance_matrix
    influence_graph /= influence_graph.max()
    #distance_reshaped = pd.Series(distance_matrix.reshape((1, distance_matrix.size))[0])
    #influence_reshaped = pd.Series(influence_graph.reshape((1, distance_matrix.size))[0])
    distance_reshaped = pd.Series(remove_diagonal(distance_matrix))
    influence_reshaped = pd.Series(remove_diagonal(influence_graph))
    correlation = distance_reshaped.corr(influence_reshaped)
    return correlation


def calculate_correlation_evolution(filename, iterations, windows_size, absolute_value=True, output_filename=None):
    pd_datas = []
    for iteration in iterations:
        correlation_evolution = []
        absolute_windows = []
        for window_size in windows_size:
            if type(window_size) is float:
                window_size = int(iteration * window_size)
            absolute_windows.append(window_size)
            ig, pos = EasyGraph.get_influence_graph_and_particles_position(filename, position_grep="position:#",
                                                                           influence_graph_grep="ig\:#[0-9]*",
                                                                           window_size=window_size,
                                                                           calculate_on=iteration)
            correlation = calculate_correlation_influence_position(ig, pos)
            correlation_evolution.append(correlation)
            print "iteration: " + str(iteration) + " with window: " + \
                  str(window_size) + " correlation: " + str(correlation)
        if absolute_value:
            pd_datas.append(pd.DataFrame({'x': absolute_windows, 'y': correlation_evolution}))
        else:
            pd_datas.append(pd.DataFrame({'x': windows_size, 'y': correlation_evolution}))
    Plotter.plot_curve(pd_datas,
                       output_filename=output_filename,
                       x_label="window size",
                       y_label="correlation with particle position",
                       legends=map(str, iterations))


def do_it():
    windows_size = [10, 100, 500, 1000]
    return read_files_and_plot([
        ('FSS1',
         '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original'),
        ('FSS2',
         '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original'),
        ('FSS3',
         '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original'),
        ('FSS4',
         '/home/marcos/PhD/research/pso_influence_graph_communities/fss_F6_original')],
        windows_size=-1,
        calculate_on=1000)

        # ('Ring',
        #  '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_ring_F6_13'),
        # ('Global',
        #  '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_global_F6_16'),
        # ('von Neumann',
        #  '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_neumann_F6_18')],
        # windows_size=-1,
        # calculate_on=1000)
    # return
    for calculate_on in range(100, 2501, 100):
        print str(calculate_on)
        read_files_and_plot(
            [
                ('Dynamic',
                 '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_dynamic_initial_ring_F6_16'),
                ('Ring',
                 '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_ring_F6_13'),
                ('Global',
                 '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_global_F6_16'),
                ('von Neumann',
                 '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_neumann_F6_18')],
            windows_size=-1,
            calculate_on=calculate_on)
    return
    windows_size = [10, 40, 50, 60, 70, 80, 90, 100]
    for iteration in range(100, 1501, 100):
        print str(iteration)
        read_files_and_plot([
            ('Dynamic',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_dynamic_initial_ring_F6_16'),
            ('Ring',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_ring_F6_13'),
            ('Global',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_global_F6_16'),
            ('von Neumann',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_neumann_F6_18')],
            windows_size,
            iteration)

    for iteration in range(500, 1501, 100):
        print str(iteration)
        windows_size = range(iteration - 4*100, iteration + 1, 100)
        read_files_and_plot([
            ('Dynamic',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_dynamic_initial_ring_F6_16'),
            ('Ring',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_ring_F6_13'),
            ('Global',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_global_F6_16'),
            ('von Neumann',
             '/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_neumann_F6_18')],
            windows_size,
            iteration)







def create_strength_distribution_curves_windows_comparison(all_graph_matrices, calculate_on, windows_size):
    print str(all_graph_matrices)
    data_hists = {}
    for title in all_graph_matrices:
        for graph_matrix in all_graph_matrices[title]:
            title_legend, graph_matrix = graph_matrix
            # mode=igraph.ADJ_MAX, otherwise igraph sums up!
            igraph_graph = igraph.Graph.Weighted_Adjacency(graph_matrix.tolist(), mode=igraph.ADJ_MAX)
            # data_hists[title] = igraph_graph.degree()
            data_hist = igraph_graph.strength(weights='weight')
            if title_legend not in data_hists:
                data_hists[title_legend] = {}
            data_hists[title_legend][title] = data_hist
            #data_hists[title] = [w/float(calculate_on) for w in data_hists[title]]
            # and the weight distribution of it
            #all_edges_weights = list(chain(*igraph_graph.get_adjacency(attribute='weight')))
            #data_hists[title] = [w for w in all_edges_weights if w != 0.0]
    font = {'family': 'normal',
            'weight': 'normal',
            'size': 8}
    matplotlib.rc('font', **font)
    # plot the results
    fig = plt.figure(figsize=(9, 6))
    plot_gridspec = gridspec.GridSpec(3, 3, width_ratios=[1, 0.001, 1], height_ratios=[1, 0.001, 1])
    graphs_index = 0
    ylabel = r"$p(X\geq x)$"
    #ylabel = "Frequency"
    xlabel = "Node strength"
    ax3 = None
    lines = ["-", "--", ":", "-."]  # matplotlib.markers.MarkerStyle.markers.keys() #
    linecycler = cycle(lines)
    titles_sort = data_hists.keys()
    titles_sort.sort()
    titles_sort.reverse()
    for title in titles_sort:
        x = graphs_index / (len(data_hists) / 2)
        y = graphs_index % (len(data_hists) / 2)
        y = 2 if y == 1 else y
        x = 2 if x == 1 else x
        ax3 = fig.add_subplot(plot_gridspec[x, y])
        ax3.set_axisbelow(True)
        for title_legend in data_hists[title]:
            # plot with power law fit:
            fit = powerlaw.Fit(data_hists[title][title_legend], discrete=True)
            fit_fig = fit.plot_ccdf(label=title_legend, linestyle=next(linecycler),
                                    linewidth=1.6) #, marker='x', markersize=2)
            fit_fig.set_xscale("linear")
            fit_fig.set_yscale("linear")
        plt.tick_params(axis='both', which='major', labelsize=9)
        plt.tick_params(axis='both', which='minor', labelsize=9)
        plt.title("Window size = " + str(title))
        plt.grid(True)
        graphs_index += 1
    plt.legend(loc=1, fontsize=8)
    ax3 = fig.add_subplot(plot_gridspec[1, 0])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_frame_on(False)
    plt.ylabel(ylabel, fontsize=13, labelpad=20)
    ax3 = fig.add_subplot(plot_gridspec[2, 1])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_frame_on(False)
    plt.xlabel(xlabel, fontsize=13, labelpad=20)
    plt.suptitle("Snapshot of the " + str(calculate_on) +
                 "th iteration", fontsize=14)
    plt.savefig('/home/marcos/PhD/research/pso_influence_graph_communities/ccdf_comp' +
                str(calculate_on) + '-' +
                '_'.join(map(str, windows_size)) + '.pdf', bbox_inches='tight')
    plt.close()
    # plt.show()


def create_strength_distribution_curves(all_graph_matrices, calculate_on):
    data_hists = {}
    for title in all_graph_matrices:
        for graph_matrix in all_graph_matrices[title]:
            title_legend, graph_matrix = graph_matrix
            # mode=igraph.ADJ_MAX, otherwise igraph sums up!
            igraph_graph = igraph.Graph.Weighted_Adjacency(graph_matrix.tolist(), mode=igraph.ADJ_MAX)
            ####
            igraph.Graph.write_graphml(igraph_graph,
                                       '/home/marcos/PhD/research/pso_influence_graph_communities/' + title + "_"
                                       + str(title_legend) + "_" + str(calculate_on) + '.graphml')
            # data_hists[title] = igraph_graph.degree()
            data_hists[title] = igraph_graph.strength(weights='weight')
            #data_hists[title] = [w/float(calculate_on) for w in data_hists[title]]
            # and the weight distribution of it
            #all_edges_weights = list(chain(*igraph_graph.get_adjacency(attribute='weight')))
            #data_hists[title] = [w for w in all_edges_weights if w != 0.0]
    #return
    font = {'family': 'normal',
            'weight': 'normal',
            'size': 8}
    matplotlib.rc('font', **font)
    # plot the results
    fig = plt.figure(figsize=(9, 6))
    plot_gridspec = gridspec.GridSpec(3, 3, width_ratios=[1, 0.001, 1], height_ratios=[1, 0.001, 1])
    graphs_index = 0
    #ylabel = r"$p(X\geq x)$"
    ylabel = "Frequency"
    xlabel = "Node strength"
    xlim = [numpy.inf, 0]
    lines = ["-", "--", "-.", ":"]  # matplotlib.markers.MarkerStyle.markers.keys() #
    linecycler = cycle(lines)
    for title in data_hists:
        xlim[0] = min(xlim[0], min(data_hists[title]))
        xlim[1] = max(xlim[1], max(data_hists[title]))
    for title in data_hists:
        x = graphs_index / (len(data_hists) / 2)
        y = graphs_index % (len(data_hists) / 2)
        y = 2 if y == 1 else y
        x = 2 if x == 1 else x
        ax3 = fig.add_subplot(plot_gridspec[x, y])
        ax3.set_axisbelow(True)
        # plot with power law fit:
        fit = powerlaw.Fit(data_hists[title], discrete=True)
        fit_fig = fit.plot_ccdf(label='Empirical Data', marker='.') #linestyle=next(linecycler))
        fit_fig.set_xscale("linear")
        fit_fig.set_yscale("linear")

        # xvalues = numpy.linspace(min(data_hists[title]), max(data_hists[title]))
        # fit.power_law.plot_ccdf(ax=fit_fig, data=xvalues, color='r', linestyle='--', label='Power law fit')
        # print title
        # R, p = fit.distribution_compare('power_law', 'exponential', normalized_ratio=True)
        # if R < 0:
        #     print 'it is a power law with p = ' + str(p)
        # elif R > 0:
        #     print 'it is not a power law with p = ' + str(p)
        # else:
        #     print 'R = 0, p = ' + str (p)


        # ecdf = sm.distributions.ECDF(data_hists[title])
        # xvalues = numpy.linspace(min(data_hists[title]), max(data_hists[title]))
        # yvalues = 1 - ecdf(xvalues)
        # plt.scatter(xvalues, yvalues, c='blue', marker='.')
        # #plt.plot(xvalues, yvalues, c='blue', marker='.')
        # #plt.plot(ecdf.x, ecdf.y)
        # ax3.set_xscale("log")
        # ax3.set_yscale("log")


        #plt.xlim((xlim[0], xlim[1]))
        # plt.scatter(data_hists[title])
        # ax3.hist(data_hists[title],
        #          bins=15,
        #          facecolor='#758AB8',
        #          #edgecolor="none"
        #          #alpha=0.45
        #          )
        #plt.ylim(0, 12)
        plt.xlim(xlim[0], xlim[1])
        plt.tick_params(axis='both', which='major', labelsize=9)
        plt.tick_params(axis='both', which='minor', labelsize=9)
        plt.title(title)
        plt.grid(True)
        graphs_index += 1
    #plt.legend(loc=4)
    ax3 = fig.add_subplot(plot_gridspec[1, 0])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_frame_on(False)
    plt.ylabel(ylabel, fontsize=13, labelpad=20)
    ax3 = fig.add_subplot(plot_gridspec[2, 1])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_frame_on(False)
    plt.xlabel(xlabel, fontsize=13, labelpad=20)
    plt.suptitle("Snapshot of the " + str(calculate_on) +
                 "th iteration", fontsize=14)
    plt.savefig('/home/marcos/PhD/research/pso_influence_graph_communities/ccdf' +
                str(calculate_on) + '.pdf', bbox_inches='tight')
    plt.close()
    #plt.show()


def create_heatmap_plot(all_graph_matrices, calculate_on):
    heatmap_dfs = {}
    for title in all_graph_matrices:
        for graph_matrix in all_graph_matrices[title]:
            title_legend, graph_matrix = graph_matrix
            heatmap_dfs[title] = pd.DataFrame(graph_matrix)
    font = {'family': 'normal',
            'weight': 'normal',
            'size': 8}
    matplotlib.rc('font', **font)
    # plot the results
    fig = plt.figure(figsize=(9, 7))
    plot_gridspec = gridspec.GridSpec(3, 3, width_ratios=[1, 0.001, 1], height_ratios=[1, 0.001, 1])
    graphs_index = 0
    ylabel = "Particles"
    xlabel = "Particles"
    xlim = [numpy.inf, 0]
    for title in heatmap_dfs:
        x = graphs_index / (len(heatmap_dfs) / 2)
        y = graphs_index % (len(heatmap_dfs) / 2)
        y = 2 if y == 1 else y
        x = 2 if x == 1 else x
        ax3 = fig.add_subplot(plot_gridspec[x, y])
        ordered = True
        matrixdf = heatmap_dfs[title]
        ordered = False
        if ordered:
            row_pairwise_dists = squareform(pdist(matrixdf))
            row_clusters = linkage(row_pairwise_dists, method='complete')
            row_dendogram = dendrogram(row_clusters, no_plot=True, count_sort='ascending')

        # calculate pairwise distances for columns
        if ordered:
            col_pairwise_dists = squareform(pdist(matrixdf.T))
            col_clusters = linkage(col_pairwise_dists, method='complete')
            col_dendogram = dendrogram(col_clusters, no_plot=True, count_sort='ascending')
        axi = ax3.imshow(matrixdf, interpolation='nearest', aspect='auto', origin='lower')
        # axi = ax3.imshow(matrixdf.ix[row_dendogram['leaves'], col_dendogram['leaves']],
        # interpolation='nearest', aspect='auto', origin='lower')

        ax3.get_xaxis().set_ticks([])
        ax3.get_yaxis().set_ticks([])
        # fit.lognormal.plot_ccdf(ax=fit_fig, color='g', linestyle='--', label='Lognormal fit')
        # fit.exponential.plot_cdf(ax=fit_fig, color='y', linestyle='--', label='Exponencial fit')
        # fit.truncated_power_law.plot_cdf(ax=fit_fig, color='b', linestyle='--', label='Truncated fit')

        # return data_hists[title]
        # ecdf = sm.distributions.ECDF(data_hists[title])
        # xvalues = numpy.linspace(min(data_hists[title]), max(data_hists[title]))
        # yvalues = 1 - ecdf(xvalues)
        # plt.scatter(xvalues, yvalues, c='blue')
        #plt.xlim((10 ** 3, 10 ** 5))
        # plt.scatter(data_hists[title])
        # ax3.hist(data_hists[title],
        #          bins=100,
        #          facecolor='blue',
        #          alpha=0.45)
        # plt.tick_params(axis='both', which='major', labelsize=9)
        # plt.tick_params(axis='both', which='minor', labelsize=9)
        plt.title(title)
        graphs_index += 1
    plt.legend(loc=4)
    ax3 = fig.add_subplot(plot_gridspec[1, 0])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_frame_on(False)
    plt.ylabel(ylabel, fontsize=13, labelpad=20)
    ax3 = fig.add_subplot(plot_gridspec[2, 1])
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_frame_on(False)
    plt.xlabel(xlabel, fontsize=13, labelpad=20)
    plt.suptitle("Snapshot of the " + str(calculate_on) +
                 "th iteration - Edges weight heatmap", fontsize=14)
    plt.savefig('/home/marcos/PhD/research/pso_influence_graph_communities/heatmap-' +
                str(calculate_on) + '.pdf', bbox_inches='tight')
    plt.close()


def read_files_and_plot(filenames, windows_size=None, calculate_on=None):
    # windows_size = [10, 40, 50, 60, 70, 80, 90, 100]
    # windows_size = [1, 5, 10, 15, 20] #, 40, 50, 60, 70, 80, 90, 100]
    # windows_size = [300, 400, 500]
    # calculates_on = [1300, 1400, 1500]
    # calculate_on = 1500
    fitness_grep = 'it\:#[0-9]*'
    influence_graph_grep = 'ig\:#[0-9]*'
    pre_callback = to_symmetric
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
            EasyGraph.read_file_and_measure(filename,
                                            calculate=None,
                                            influence_graph_grep=influence_graph_grep,
                                            fitness_grep=fitness_grep,
                                            window_size=window_size,
                                            pre_callback=pre_callback,
                                            pos_callback=pos_callback,
                                            calculate_on=calculate_on)
        all_graph_matrices[title] = graph_matrices
        ### create the GiantComponentDeath analysis
    create_giant_component_curves(all_graph_matrices, calculate_on, windows_size)
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



def str_and_print(obj):
    print (str(obj))
            

def from_igraph_to_nxgraph(igraph_graph, only_connected_nodes=False):
    nxgraph = nx.DiGraph()
    adjacency = list(igraph_graph.get_adjacency(attribute='weight'))
    if not only_connected_nodes:
        for i in range(len(adjacency)):
            nxgraph.add_node(i)
    for i in range(len(adjacency)):
        for j in range(len(adjacency)):
            if adjacency[i][j] != 0:
                nxgraph.add_edge(i, j, weight=adjacency[i][j])
    return nxgraph


def create_and_save_plot(matrix_out, fitness=None):
    ### FILENAME
    base_dir = '/home/marcos/PhD/research/pso_influence_graph_communities/50_particles/heatmaps/'
    if not ('append_filename' in globals()):
        append = ''
    else:
        # this is too dirty, but ok...
        global append_filename
        append = append_filename + '_'
    output_filename = base_dir+append+str(matrix_out[0])+'.png'
    print output_filename

    ### create the GiantComponentDeath analysis
    return_graphs_with_giant_sizes = [1.0, 0.9, 0.7, 0.5, 0.3]
    matrix = to_symmetric(matrix_out[1])
    igraph_graph = igraph.Graph.Weighted_Adjacency(matrix.tolist())

    # create the graph objects as well as the death analysis
    pd_data_2, graphs = GiantComponentDeath.low_edges_weight_removal(igraph_graph, return_graphs_with_giant_sizes)
    #pd_data_2, graphs = GiantComponentDeath.nodes_degree_removal(igraph_graph, return_graphs_with_giant_sizes)

    pd_data_2 = ('Giant Comp. Death', pd_data_2)

    # but 'graphs' is actually igraph.graph, but we need
    # networkx graphs, dammit! (because just nx.graph can be plot with matplotlib :( -- it seems)
    nx_graphs = []
    graph = None
    for graph in graphs:
        graph_component_histogram = graph[2].components().sizes()
        nx_graph = from_igraph_to_nxgraph(graph[2], only_connected_nodes=True)
        title = str(graph[1]) + " ("+str(graph[0])+") [" + str(nx.number_of_nodes(nx_graph)) \
                              + "/" + str(graph[2].vcount()) + "]"
        nx_graphs.append((title, nx_graph, graph_component_histogram))
    if not nx_graphs:
        nx_graphs = None

    ### here is the fitness data
    pd_data_1 = None
    if fitness is not None:
        pd_data_1 = pd.DataFrame({'x': range(len(fitness)), 'y': fitness})
        pd_data_1 = ('Fitness', pd_data_1)

    ### create the histograms data
    # gets the last graph in 'graphs' and plot the degree distribution of it
    data_hist_1 = ('Degree Distribution', graph[2].degree())
    all_edges_weights = list(chain(*graph[2].get_adjacency(attribute='weight')))
    # and the weight distribution of it
    data_hist_2 = ('Weight Distribution', [w for w in all_edges_weights if w != 0.0])
    Plotter.create_heatmap(matrix_out[1], main_title=str(matrix_out[0]),
                           output_filename=output_filename,
                           pd_data_1=pd_data_1,
                           pd_data_2=pd_data_2,
                           data_hist_1=data_hist_1,
                           data_hist_2=data_hist_2,
                           graphs=nx_graphs)





'''
execfile('/home/marcos/pyworkspace/SwarmParser/EasyGraphpy')
append_filename = '
   without window:
   BADpso_dynamic_initial_ring_F6_16:it:#2999     2140695.499368374   OK sym
   GOODpso_dynamic_initial_ring_F6_02:it:#2999     21.318622649096547  OK sym and non-sym
   BADpso_ring_F6_13:it:#2999                     2372248.830220812   OK sym
   GOODpso_ring_F6_21:it:#2999                     1155170.4496039404  OK sym
   BADpso_global_F6_16:it:#2999                   4905393.756248348   OK sym
   GOODpso_global_F6_07:it:#2999                   1561303.7169722272  OK sym
   BADpso_static_neumann_04
   GOODpso_static_neumann_18
   BADpso_random_F6_07
   GOODpso_random_F6_12



execfile('/home/marcos/pyworkspace/SwarmParser/easygraph.py)
append_filename = 'pso_dynamic_ring_16_noW'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_dynamic_initial_ring_F6_16',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/EasyGraphpy')
append_filename = 'pso_dynamic_ring_02_noW'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_dynamic_initial_ring_F6_02',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/easygraph.py)
append_filename = 'pso_static_ring_13_noW'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_ring_F6_13',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/easy_graph/easygraph.py)
append_filename = 'pso_static_ring_21_noW'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_ring_F6_21',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/easygraph.py)
append_filename = 'pso_static_global_16_noW'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_global_F6_16', calculate = None,
fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/EasyGraphpy')
append_filename = 'pso_static_global_07_noW'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_global_F6_07',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/easygraph.py)
append_filename = 'pso_static_neumann_18_w10'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_neumann_F6_18',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/EasyGraphpy')
append_filename = 'pso_static_neumann_04_w10'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_neumann_F6_04',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/EasyGraphpy')
append_filename = 'pso_random_F6_12_noW'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_random_F6_12',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)

execfile('/home/marcos/pyworkspace/SwarmParser/EasyGraphpy')
append_filename = 'pso_random_F6_07_noWs'
SwarmParser.read_file_and_measure(
'/home/marcos/PhD/research/pso_influence_graph_communities/100_particles/pso_random_F6_07',
calculate = None, fitness_grep = 'it\:#[0-9]*', influence_graph_grep = 'ig\:#[0-9]*',
pre_callback = to_symmetric, pos_callback = create_and_save_plot)
'''
