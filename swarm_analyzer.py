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
    def calculate_velocities_correlation(
            filename, dimensions=100, particles=1000, informations_grep="velocities\:#", kind='all', kind_average='all',
            absolute=False, save_hdf=None, **kargs):
        information_map = lambda x: np.array(map(float, x.split())).reshape(particles, dimensions)
        velocities = SwarmParser.read_file_and_measures(
            filename, influence_graph_grep=None,
            informations_grep=informations_grep, information_map=information_map, **kargs)
        velocities = [v[1] for v in velocities[1][informations_grep]]
        # normalize all velocities at each iteration.
        velocities = [v / np.linalg.norm(v) for v in velocities]
        correlation_t = []
        v_sum = []
        v_average = []

        if kind_average in ['all']:
            v_sum = reduce(add, (velocities[it] for it in range(len(velocities))))
            v_average = v_sum / float(len(velocities))
        elif kind_average in ['iteration']:
            v_sum = copy.deepcopy(velocities[0])

        for iteration in range(1, len(velocities)):
            if kind in ['all']:
                correlations = pd.DataFrame(np.rot90(velocities[iteration])).corr()
                # correlation = np.array(correlations).reshape(1, correlations.shape[0] * correlations.shape[1])[0]
                correlation = np.ravel(correlations)
                correlation_t.append(correlation)
            elif kind in ['average', 'fluctuations']:
                if kind_average in ['iteration']:
                    v_sum += velocities[iteration]
                    v_average = v_sum / float(iteration)
                correlations = None
                if kind == 'fluctuations':
                    # if iteration < len(velocities) - 1:
                    fluctuations = [velocities[iteration][p] - v_average[p] for p in range(particles)]
                    correlations = pd.DataFrame(np.rot90(fluctuations)).corr()
                else:
                    # average
                    correlations = pd.DataFrame(np.rot90(v_average)).corr()
                if correlations is not None:
                    # correlation = np.array(correlations).reshape(1, correlations.shape[0] * correlations.shape[1])[0]
                    correlation = np.ravel(correlations)
                    correlation_t.append(correlation)
        if save_hdf is not None:
            df = pd.DataFrame(correlation_t)
            df.to_hdf(save_hdf, 'df')
        return correlation_t

    """


DIEGO
import matplotlib.pyplot as plt

execfile("swarm_analyzer.py")
filename = 'watts_strogatz_1.00000_F23_29_min'
filename_hdf = 'watts_strogatz_1.00000_F23_29.hdf'
correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind='fluctuations', save_hdf=filename_hdf)
for particle in range(0):
    Plotter.plot_curve({'x': range(len(correlation_t[particle])), 'y': correlation_t[particle]}, dpi=72, figsize=(20, 5), tight_layout=[], x_label="Iteration", y_label="pearson correlation", title="Particle #%d" % particle, output_filename="plot.png")

filename_hdf = 'watts_strogatz_1.00000_F23_29.hdf'
df = pd.read_hdf(filename_hdf, 'df')

fig, a = plt.subplots()

for it in range(0,2):
        print it
        values = df.irow(it)
        individuals = int(np.sqrt(len(values)))
        unique = [values[j*100 + i] for i in range(individuals) for j in range(i+1, individuals)]
        unique = np.abs(unique)
        plt.hist(unique)
        #plt.xlim((-1, 1))
        plt.show()
        ax = fit.plot_ccdf(ls='', marker='.')
        ax = fit.power_law.plot_ccdf(ax=ax)
        ax = fit.lognormal.plot_ccdf(ax=ax)
        ax = fit.exponential.plot_ccdf(ax=ax)
        ax = fit.truncated_power_law.plot_ccdf(ax=ax)
        plt.show()

import matplotlib
matplotlib.use('Agg')
execfile("swarm_analyzer.py")
topologies = [("kregular%d" % i) for i in [3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90]]
topologies += ['noc2']
topologies += ['global']
topologies += ['ring']
topologies = [("kregular%d" % i) for i in range(12, 20, 2) + range(22, 30, 2)]
for topology in topologies:
    for function in [23]:
        for run in [0]:
            filename = "./%s_F%02d_%02d.with_positions" % (topology, function, run)
            kind = "average"
            print filename
            filename_hdf = filename + "_%s_correlation.hdf" % kind
            correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind=kind, save_hdf=filename_hdf)
            del correlation_t

    filename = "./data/100_particles/regular30_F21_00.with_positions_head_30"
    kind = "average"
    correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind=kind, until=until)
    for particle in range(100):
        Plotter.plot_curve({'x': range(len(correlation_t[particle])), 'y': correlation_t[particle]}, dpi=72, figsize=(20, 5), tight_layout=[], x_label="Iteration", y_label="%s pearson correlation" % kind, title="Particle #%d" % particle, output_filename="%s_particle_%d_%s_F%02d_%02d.png" % (kind, particle, topology, function, run))

    import powerlaw
    filename = "regular30_F21_00.with_positions_fluctuations_correlation"
    df = pd.read_hdf(filename + ".hdf", 'df')
    results = []
    for it in range(100, 1000, 100):
        print it
        values = df.irow(it)
        individuals = int(np.sqrt(len(values)))
        unique = [values[j*100 + i] for i in range(individuals) for j in range(i+1, individuals)]
        unique = np.abs(unique)
        plt.hist(unique)
        #plt.xlim((-1, 1))
        plt.show()
        ax = fit.plot_ccdf(ls='', marker='.')
        ax = fit.power_law.plot_ccdf(ax=ax)
        ax = fit.lognormal.plot_ccdf(ax=ax)
        ax = fit.exponential.plot_ccdf(ax=ax)
        ax = fit.truncated_power_law.plot_ccdf(ax=ax)
        plt.show()

topology, function = "kregular30", 23
topology, function = "regular30", 21
topology = "noc2"
topology, function = "global", 23
function = 21
function = 23
run = 0
kind = "fluctuations"
filename = "./%s_F%02d_%02d.with_positions" % (topology, function, run)
filename_hdf = filename + "_%s_correlation.hdf" % kind
print filename_hdf
df = pd.read_hdf(filename_hdf, 'df')

from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib
matplotlib.use('Agg')
execfile("swarm_analyzer.py")


x = np.arange(-1, 1, 0.01)
curves.append({'x': x, 'y': ecdf(x)})

all_values[0]
df.tail(0)
d
values = matrix[0]
iteration = 1520

iteration = 1000
iteration = 3988 # flat
iteration = 3970 # v
iteration = 3971 # v
iteration = 3973 # v
iteration = 3710 # v
iteration = 4190 # 0 with tail
iteration = 4500 # 0 not much tail
iterations = [3988, 3971, 4190, 4500]
Plotter.plos_style()
for iteration in iterations:
    values = np.array(df.irow(iteration))
    values = values.reshape(100, 100)
    values = values[np.triu_indices(n=100, k=1)]
    print iteration
    plt.figure(figsize=(1.65, 1.3))
    plt.hist(values, bins=15, facecolor='#a6bddb')
    plt.xlim(-1, 1)
    plt.ylabel("Frequency")
    plt.xlabel("Correlation")
    plt.tight_layout(rect=[-0.08, -0.11, 1.09, 1.09])
    # plt.savefig("%s_it%d_correlation_%s.pdf" % (filename, iteration, kind))
    plt.savefig("histogram_%d.pdf" % (iteration))
    # plt.clf()

    plt.show()
    """

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

    """
execfile("swarm_analyzer.py")
topology = "kregular"
filename = "./data/100_particles/"+topology+"_F06_06.teste"
until = 2000
ws = list(np.arange(0, 1.01, 0.01))
tws = range(1, 1001, 10) + [1000]
components = dict([(i, np.zeros((len(ws), len(tws)))) for i in range(until + 1)])
for i in range(len(ws)):
    for j in range(len(tws)):
        w = ws[i]
        tw = tws[j]
        for it, cs in SwarmAnalyzer.get_number_of_components(filename, tw, w, calculate_on=1000):
            components[it][i][j] = cs
        print i, j

print "Exporting: "
for it in components:
    pd.DataFrame(components[0]).to_hdf("./data/100_particles/number_of_components/heatmap_1000_%s_%05d.hdf" % (topology, it), 'df')


execfile("swarm_analyzer.py")
execfile("plotter.py")

filename = "./data/100_particles/vonneumann_F06_06.teste"
window_size = 100
w = 0.2
SwarmAnalyzer.get_number_of_components(filename, window_size, w, until=1000)

    execfile("swarm_analyzer.py")
    execfile("plotter.py")

    window_size = 100
    until = 3000
    ws = [0.7, 0.5, 0.2]
    ws = [0.9, 0.7, 0.5, 0.3]

    filename = "./data/100_particles/vonneumann_F06_06.teste"
    curves_1 = []
    for w in ws:
        components = SwarmAnalyzer.get_number_of_components(filename, window_size, w, until=until)
        curve = {}
        curve['x'] = []
        curve['y'] = []
        for x, y in components:
            curve['x'].append(x)
            curve['y'].append(y)
        curves_1.append(curve)

    filename = "./data/100_particles/ring_F06_06.teste"
    curves_2 = []
    for w in ws:
        components = SwarmAnalyzer.get_number_of_components(filename, window_size, w, until=until)
        curve = {}
        curve['x'] = []
        curve['y'] = []
        for x, y in components:
            curve['x'].append(x)
            curve['y'].append(y)
        curves_2.append(curve)
    legends = [str(w) for w in ws]

    # smoothing the curves: getting points every 15 iterations
    curvess1 = []
    for c in curves_1:
        curve = {}
        curve['x'] = [c['x'][i] for i in range(0, len(c['x']), 15)]
        curve['y'] = [c['y'][i] for i in range(0, len(c['y']), 15)]
        curvess1.append(curve)

    curvess2 = []
    for c in curves_2:
        curve = {}
        curve['x'] = [c['x'][i] for i in range(0, len(c['x']), 15)]
        curve['y'] = [c['y'][i] for i in range(0, len(c['y']), 15)]
        curvess2.append(curve)
    execfile("plotter.py")
plt.clf()
    #Plotter.plot_curve(curves, legends=legends, marker=['.', '+', 'x'], markersize=4, font_size=10, markevery=17, alpha=0.9, tight_layout=[], figsize=(5,4), grid=True, linewidth=0.7, xlim=(window_size, 1400), y_label="Number of\ncomponents", x_label="Iterations")
    Plotter.plot_subplot_curve([curvess2, curvess1], titles=['Von Neumann', 'Ring'], legends=legends, marker=['.', '+', 'x'], markersize=4, font_size=13, markevery=17, alpha=0.9, tight_layout=[], figsize=(13, 4), grid=True, linewidth=0.7, xlim=(window_size, 3000), y_label="Number of\ncomponents", x_label="Iterations")
    Plotter.plot_subplot_curve([curvess2, curvess1], colors=['#009e73', '#0072b2', '#d55e00', '#e69f00'], titles=['Ring', 'Von Neumann'], legends=legends, marker=['v', 's', 'o', 'D'], markevery=10, markersize=3, mew=1.2, font_size=13, tight_layout=[-0.03, 0, 1.03, 1.06], figsize=(4.5, 1.9), grid=True, linewidth=1.2, xlim=(window_size, 2000), y_label="Number of components", x_label="Iterations", ylim=(0, 102), output_filename="number_of_components.pdf")
    Plotter.plot_subplot_curve([curves, curves_2], titles=['Von Neumann', 'Ring'], legends=legends, marker=['.', '+', 'x'], markersize=4, font_size=13, markevery=17, alpha=0.9, tight_layout=[], figsize=(13, 4), grid=True, linewidth=0.7, xlim=(window_size, 3000), ylim=(0, 17), y_label="Number of\ncomponents", x_label="Iterations", output_filename="number_of_components.pdf")
    """

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
    """
execfile("swarm_analyzer.py")
# plotting measures vs. fitness:
filename = "./data/100_particles/ring_F06_06.teste"
measures = ["radius:#", "aggregation_factor:#", "it:#", "average_of_average_around_all_particles:#", "normalized_average_around_center:#", "average_of_average_around_all_particles:#"]

df_info = SwarmAnalyzer.get_swarm_informations_from_file(filename, measures)
dfs = [{'y': df_info[k]/max(df_info[k]), 'x': df_info['x']} for k in measures]
Plotter.plot_curve(dfs, legends=measures, markersize=0, markevery=10, figsize=(20,6), grid=True, linewidth=1, ylim=(0, 1.0))


# plotting measures vs. fitness vs. AUC:
filename = "./data/100_particles/ring_F06_06.teste"
measures = ["radius:#", "aggregation_factor:#", "it:#", "average_of_average_around_all_particles:#", "normalized_average_around_center:#", "average_of_average_around_all_particles:#"]

df_info = SwarmAnalyzer.get_swarm_informations_from_file(filename, measures)
dfs = [{'y': df_info[k]/max(df_info[k]), 'x': df_info['x']} for k in measures]

auc_100 = SwarmAnalyzer.get_giant_component_destruction_area(filename, window_size=100, until=10000)
print time.localtime()
auc_200 = SwarmAnalyzer.get_giant_component_destruction_area(filename, window_size=200, until=10000)
print time.localtime()
measures = []
dfs = []
dfs += [{'y': auc_100['y'], 'x': auc_100['x']}]
dfs.append(SwarmAnalyzer.difference_n(auc_100, 2))
dfs.append(SwarmAnalyzer.difference_n(auc_200, 2))
# dfs += [{'y': auc_200['y'], 'x': auc_200['x']}]
# legends = measures + ['auc100', 'auc200']
legends = measures + ['auc100', 'auc200']
Plotter.plot_curve(dfs, legends=legends, markersize=0, markevery=10, figsize=(20,6), grid=True, linewidth=1)


execfile("swarm_analyzer.py")
execfile("swarm_parser.py")
execfile("plotter.py")

measures = ["distance:#"]
filename = "./global_F06_00.position"
filename = "./ring_F06_00.position"
filename = "./vonneumann_F06_00.position"
pre_callback = Callback.to_symmetric
iter = 999
# for iter in [999]:
for iter in [999]:
    graphs, _ = SwarmParser.read_file_and_measures(filename, influence_graph_grep="ig\:#", window_size=iter, calculate_on=iter, pos_callback=pre_callback)

    graph_1_ = graphs[0][1]
    df_info = SwarmAnalyzer.get_swarm_informations_from_file(filename, measures, information_map=lambda x: x, calculate_on=iter)  # position

    graph_2_ = SwarmParser.read_matrix_from_line(list(df_info.irow(0))[1])  # position
    graph_1 = [graph_1_[i][j] for i in range(len(graph_1_)) for j in range(i, len(graph_1_)) if i != j] # np.reshape(graph_1_, (1, 10000))  # graph
    graph_2 = [graph_2_[i][j] for i in range(len(graph_2_)) for j in range(i, len(graph_2_)) if i != j] # np.reshape(graph_2_, (1, 10000))  # position

    graph_1 = np.array(graph_1)  # graph
    graph_2 = np.array(graph_2)  # position

    graphs, _ = SwarmParser.read_file_and_measures(filename, influence_graph_grep="ig\:#", window_size=iter, calculate_on=iter)
    igraph_graph = igraph.Graph.Weighted_Adjacency(graphs[0][1].tolist(), mode=igraph.ADJ_PLUS)
    max_weight = max(igraph_graph.es['weight'])
    for w in igraph_graph.es:
        w['weight'] = max_weight - w['weight']
    distances = {}
    for i in igraph_graph.vs:
        if i not in distances:
            distances[i.index] = {}
        for j in igraph_graph.vs:
            distances[i.index][j.index] = igraph_graph.shortest_paths(source=i, target=j, weights='weight')[0][0]
    indices = distances.keys()
    indices.sort()
    final_distances = [[distances[i][j] if distances[i][j] != float('inf') else 100.0 for j in indices] for i in indices]
    final_distances = np.array(final_distances)
    graph_3 = [final_distances[i][j] for i in range(len(final_distances)) for j in range(i, len(final_distances)) if i != j]  # np.reshape(final_distances, (1, 10000))
    graph_3 = np.array(graph_3)

    import pandas as pd
    # df = pd.DataFrame({'x': graph_1[0]/max(graph_1[0]), 'y': graph_2[0]/max(graph_2[0])})
    #df = pd.DataFrame({'x': graph_1[0]/max(graph_1[0]), 'y': graph_3[0]/float(max(graph_3[0]))})
    df = pd.DataFrame({'y': graph_2/float(max(graph_2)), 'x': graph_3/float(max(graph_3))})
    xlabel = "Geodesic distance"
    ylabel = "Euclidean distance"
    title = "Ring"
    output_filename="ring_%04d.pdf" % iter
    output_filename = None
    Plotter.plot_curve(df, linewidth=0, figsize=(5, 3.8), tight_layout=[], output_filename=output_filename, annotate=(0.65, 0.35, "$r=%0.2f$" % (df.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors='#de2d26', grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel, xlim=(0,1), ylim=(0,1))
    # print "%03d:  %f" % (iter, df.corr()['y']['x'])


Plotter.plot_curve(df, linewidth=0, figsize=(5, 4), tight_layout=[], alpha=0.6, marker='.', annotate=(0.15, 0.84, "$R=%0.2f$" % (df.corr()['y']['x'])), font_size=11, colors='#de2d26', grid=True, markersize=6, title="Global", x_label=xlabel, y_label=ylabel, xlim=(0,1), ylim=(0,1))

graph_1 = graph_1_
import igraph
igraph_graph = igraph.Graph.Weighted_Adjacency(graph_1_.tolist(), mode=igraph.ADJ_PLUS)
import igraph
igraph_graph = igraph.Graph.Weighted_Adjacency(graph_1_.tolist(), mode=igraph.ADJ_PLUS)
distances = {}
for i in igraph_graph.vs:
    if i not in distances:
        distances[i.index] = {}
    for j in igraph_graph.vs:
        distances[i.index][j.index] = igraph_graph.shortest_paths(source=i, target=j, weights='weight')
        break

max_weight = max(igraph_graph.es['weight'])
for w in igraph_graph.es:
    w['weight'] = max_weight - w['weight']


indices = distances.keys()
indices.sort()
final_distances = [[distances[i][j] if distances[i][j] != float('inf') else 100.0 for j in indices] for i in indices]
final_distances = np.array(final_distances)
graph_3 = np.reshape(final_distances, (1, 10000))

    """

    @staticmethod
    def difference_n(df, n):
        array = np.array(df['y'])
        array = array[:len(array) - n] - array[n:]
        data = {'y': list(array), 'x': list(df['x'][:len(df['x']) - n])}
        return data

    @staticmethod
    def get_giant_component_destruction_area_from_files(basename, window_size, runs=30):
        filenames = [basename + "%02d" % i for i in range(1, runs + 1)]
        df = None
        run = 1
        for filename in filenames[:1]:
            df = SwarmAnalyzer.get_giant_component_destruction_area(filename, window_size=window_size)
        df.columns = ['x', "%02d" % run]
        run += 1
        for filename in filenames[1:]:
            df_run = SwarmAnalyzer.get_giant_component_destruction_area(filename, window_size=window_size)
            df["%02d" % run] = df_run['y']
            run += 1
            del df_run
        return df

    """
    execfile("swarm_analyzer.py")
    basename = "/mnt/pso_100_particles/global_F06_"
    df = SwarmAnalyzer.get_giant_component_destruction_area_from_files(basename, 100, runs=3)
    """

    @staticmethod
    def export_giant_component_destruction_areas(basename, window_size, runs=30):
        df = SwarmAnalyzer.get_giant_component_destruction_area_from_files(basename, runs=runs, window_size=window_size)
        df.to_hdf(basename + str(window_size) + ".hdf", 'df')

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
        graph_matrices, _ = SwarmParser.read_files_and_measures(filenames, influence_graph_grep=influence_graph_grep,
                                                                pos_callback=pre_callback, windows_size=windows_size,
                                                                calculate_on=calculate_on)
        normalize = [2 * i for i in windows_size]
        pd_datas = []
        for title, _ in filenames:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(graphs,
                                                                             weight_normalize=normalize)
            pd_datas.append((title, dict(zip(windows_size, curves_areas))))
        GiantComponentDeathPlotter.giant_component_death_curve(calculate_on, pd_datas, windows_size, xlim=(0, 1.0),
                                                               figsize=(4.5, 4))

    """
    execfile("swarm_analyzer.py")
    execfile("giant_component_analysis_plotter.py")
    filenames = [('Global', "./data/100_particles/global_F06_15.teste"), ('Ring', "./data/100_particles/ring_F06_15.teste"), ('Von Neumann', "./data/100_particles/vonneumann_F06_15.teste"), ('Dynamic', "./data/100_particles/dynamicring_F06_15.teste")]
    df = SwarmAnalyzer.read_files_and_plot(filename, windows_size=[10, 50, 100, 1000], calculate_on=1000, tight_layout=[])
    """

    @staticmethod
    def read_files_and_plot_with_area(filenames, windows_size, calculate_on, **kargs):
        influence_graph_grep = 'ig\:#'
        pre_callback = Callback.to_symmetric
        graph_matrices, _ = SwarmParser.read_files_and_measures(filenames, influence_graph_grep=influence_graph_grep,
                                                                pos_callback=pre_callback, windows_size=windows_size,
                                                                calculate_on=calculate_on)
        normalize = [2 * i for i in windows_size]
        pd_datas = []
        for title, _ in filenames:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(graphs,
                                                                             weight_normalize=normalize)
            pd_datas.append((title, dict(zip(windows_size, curves_areas))))
        GiantComponentDeathPlotter.giant_component_death_curve_with_area(pd_datas, xlim=(0, 1.0), figsize=(4.5, 1.9),
                                                                         mew=1.2,
                                                                         tight_layout=[-0.02, 0.01, 1.02, 1.06],
                                                                         output_filename="graph_destruction_area.pdf")  # , **kargs)

    """
    plt.clf()
    execfile("swarm_analyzer.py")
    execfile("giant_component_analysis_plotter.py")
    filenames = [('Global', "./data/100_particles/global_F06_06.teste"), ('Ring', "./data/100_particles/ring_F06_06.teste")]
    SwarmAnalyzer.read_files_and_plot_with_area(filenames, windows_size=[100, 1000], calculate_on=1000, output_filename="graph_destruction_area_10002.pdf")
    SwarmAnalyzer.read_files_and_plot_with_area(filenames, windows_size=[100, 1000], calculate_on=1000)
    """
    # create_strength_distribution_curves_windows_comparison(all_graph_matrices, calculate_on, windows_size)
    # create_heatmap_plot(all_graph_matrices, calculate_on)
    # create_strength_distribution_curves(all_graph_matrices, calculate_on)

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
    # return graph_matrix


# if __name__ == "__main__":
#    SwarmAnalyzer.read_files_and_export_hdf(int(sys.argv[1]), sys.argv[2])

"""
execfile("swarm_analyzer.py")
SwarmAnalyzer.do_it()
python2.7 ./swarm_analyzer.py 1 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ && python2.7 ./swarm_analyzer.py 5 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 10 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ && python2.7 ./swarm_analyzer.py 25 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 50 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ && python2.7 ./swarm_analyzer.py 100 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 250 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ && python2.7 ./swarm_analyzer.py 500 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 750 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ && python2.7 ./swarm_analyzer.py 1000 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 2000 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ && python2.7 ./swarm_analyzer.py 3000 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &


python2.7 ./swarm_analyzer.py 1 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ && python2.7 ./swarm_analyzer.py 5 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 10 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ && python2.7 ./swarm_analyzer.py 25 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 50 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ && python2.7 ./swarm_analyzer.py 100 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 250 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ && python2.7 ./swarm_analyzer.py 500 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 750 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ && python2.7 ./swarm_analyzer.py 1000 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 2000 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ && python2.7 ./swarm_analyzer.py 3000 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &



python2.7 ./swarm_analyzer.py 1 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ && python2.7 ./swarm_analyzer.py 5 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &
python2.7 ./swarm_analyzer.py 10 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ && python2.7 ./swarm_analyzer.py 25 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &
python2.7 ./swarm_analyzer.py 50 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ && python2.7 ./swarm_analyzer.py 100 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &
python2.7 ./swarm_analyzer.py 250 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ && python2.7 ./swarm_analyzer.py 500 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &
python2.7 ./swarm_analyzer.py 750 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ && python2.7 ./swarm_analyzer.py 1000 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &
python2.7 ./swarm_analyzer.py 2000 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ && python2.7 ./swarm_analyzer.py 3000 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &
p
python2.7 ./swarm_analyzer.py 200 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 200 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 200 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &

python2.7 ./swarm_analyzer.py 150 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 150 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 150 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &

python2.7 ./swarm_analyzer.py 125 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 125 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 125 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &

python2.7 ./swarm_analyzer.py 175 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 175 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 175 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &

python2.7 ./swarm_analyzer.py 225 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ &
python2.7 ./swarm_analyzer.py 225 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ &
python2.7 ./swarm_analyzer.py 225 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ &


python2.7 ./swarm_analyzer.py 1 /home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_ >> output &

python2.7 ./swarm_analyzer.py 1 /home/moliveira/pso/50_particles_simulations/pso_global_F6_ >> output &

python2.7 ./swarm_analyzer.py 1 /home/moliveira/pso/50_particles_simulations/pso_ring_F6_ >> output &

we_have = [5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 500, 750]

all = range(5, 1000, 5)
map(lambda x: all.remove(x), we_have)

iterations = []
line = 35
for i in range(0, (len(all)/line)*line, line):
    iterations.append([all[i+k] for k in range(line)])
iterations.append([all[i] for i in range((len(all)/line)*line, len(all))])

c = "/home/moliveira/pso/50_particles_simulations/pso_dynamic_initial_ring_F6_"
#c = "/home/moliveira/pso/50_particles_simulations/pso_global_F6_"
#c = "/home/moliveira/pso/50_particles_simulations/pso_ring_F6_"
for iters in iterations:
    command = "python2.7 ./swarm_analyzer.py %d %s" % (iters[0], c)
    for i in iters[1:]:
        command += " && python2.7 ./swarm_analyzer.py %d %s " % (i, c)
    print command + " & "

"""

"""
import numpy as np
import matplotlib.pyplot as plt
t = np.arange(-50, 50, 0.012)
sp = np.fft.fft(np.sin(t))
freq = np.fft.fftfreq(t.shape[-1])
plt.plot(freq, sp.real, freq, sp.imag)

func = np.sin(t)
func = np.cos(2*np.pi*(3*t))*np.exp(-np.pi*t*t)
sp = np.fft.fft(func)

import pandas as pd
plt.plot(t, func)
plt.plot(t, sp.)
plt.show()
plt.legend(loc=1)


# #, 'y2': func})
plt.plot(t, np.sin(t))
"""
"""
execfile("swarm_analyzer.py")
df3 = SwarmAnalyzer.get_giant_component_destruction_area_from_files("/mnt/50_particles_simulations/pso_dynamic_initial_ring_F6_", 1000, runs=2)
"""
"""
execfile("swarm_analyzer.py")
basename = "/mnt/50_particles_simulations/hdfs/pso_global_F6_"
basename = "/mnt/50_particles_simulations/hdfs/pso_ring_F6_"
basename = "/mnt/50_particles_simulations/hdfs/pso_dynamic_initial_ring_F6_"
SwarmAnalyzer.read_files_and_export_hdf(100, basename)
SwarmAnalyzer.read_hdfs_and_plot(basename)
window_size = 10
filename1 = "/mnt/50_particles_simulations/pso_dynamic_initial_ring_F6_16"
df1 = SwarmAnalyzer.get_giant_component_destruction_area(('None', filename1), window_size)
df1.to_hdf(filename1+"_"+str(window_size)+".hdf", 'df')
filename2 = "/mnt/50_particles_simulations/hdfs/pso_ring_F6_13"
df2 = SwarmAnalyzer.get_giant_component_destruction_area(('None', filename2), window_size)
df2.to_hdf(filename2+"_"+str(window_size)+".hdf", 'df')
filename3 = '/mnt/50_particles_simulations/hdfs/pso_global_F6_16'
df3 = SwarmAnalyzer.get_giant_component_destruction_area(('None', filename3), window_size)
df3.to_hdf(filename3+"_"+str(window_size)+".hdf", 'df')

Plotter.plot_curve([df1, df2, df3], figsize=(18, 6), legends=['dynamic', 'ring', 'global'])

filename = filename3
dfZ = pd.read_hdf(filename+"_"+str(10)+".hdf", 'df')
dfA = pd.read_hdf(filename+"_"+str(50)+".hdf", 'df')
dfB = pd.read_hdf(filename+"_"+str(100)+".hdf", 'df')
dfC = pd.read_hdf(filename+"_"+str(500)+".hdf", 'df')
dfD = pd.read_hdf(filename+"_"+str(1000)+".hdf", 'df')

Plotter.plot_curve([dfA, dfB, dfC, dfD], figsize=(18, 6), legends=['50', '100', '500', '1000'], linestyle="-", markersize=2, linewidth=.2)

dfff = pd.read_hdf("/mnt/50_particles_simulations/pso_ring_F6_13_1000.hdf", 'df')
Plotter.plot_curve([dfA, dfB], figsize=(18, 6), legends=['500', '1000'])

delta = 0.001
number_particles = 50
normalization = number_particles*len(np.arange(0, 1 + delta, delta))

execfile("swarm_analyzer.py")

basename = "/mnt/50_particles_simulations/hdfs/"
filename = "/mnt/50_particles_simulations/hdfs/pso_global_F6_10"
filename = "/mnt/50_particles_simulations/hdfs/pso_dynamic_initial_ring_F6_26"
filename = "/mnt/50_particles_simulations/hdfs/pso_ring_F6_10"
filenames = ["pso_global_F6_%02d" % i for i in range(1, 30)]
filenames = ["pso_dynamic_initial_ring_F6_%02d" % i for i in range(1, 30)]
filenames = ["pso_ring_F6_%02d" % i for i in range(1, 30)]
for filename in filenames[:6]:
    iterations = [5, 10, 25, 50, 75, 100]
    iterations = [500, 501, 502, 503]
    #iterations = [300, 350, 400, 450, 500, 550, 600, 1000, 2000]
    #iterations = [500, 1000, 2000, 3000] #, 3000, 4000, 5000]
    areass = []
    legends = []
    for iteration in iterations:
        areas = {'x': [], 'y': []}
        for tw in range(5, 1000, 5) + [2000, 3000, 4000, 5000]:
            if tw <= iteration:
                df = pd.read_hdf(basename+filename+"_"+str(tw)+".hdf", 'df')
                areas['x'].append(tw)
                areas['y'].append(float(df.y[df.x == iteration])/normalization)
        if areas['x']:
            areass.append(areas)
            legends.append(iteration)
    ylim = (0.5,0.7)
    ylim = (0.5,1.0)
    Plotter.plot_curve(areass, figsize=(23,7), output_filename="tws_"+filename+".png", legends=legends, x_label="Time Window", markersize=7, linewidth=2, alpha=0.7, grid=True, ylim=ylim)
xticks = range(5, 1000, 5)

Plotter.plot_curve(areass[0],  xticks_args=([100, 350], map(str, [100, 350])), legends=legends, figsize=(23,7), grid=True)

Plotter.plot_curve(grid=True, xticks_args=[xticks], legends=legends, linestyle=" ", markersize=4, linewidth=.2, marker=markers, colors=colors, ylim=(0,1), xlim=(0,500), alpha=0.7)

iterations = [5, 10, 25, 50, 75, 100]
iterations = [100, 200, 300, 400, 500, 600]
iterations = range(100, 1000, 50)
#iterations = [300, 350, 400, 450, 500, 550, 600, 1000, 2000]
#iterations = [500, 1000, 2000, 3000] #, 3000, 4000, 5000]
filename = "pso_global_F6_01"
filename = "pso_dynamic_initial_ring_F6_01"
filename = "pso_ring_F6_01"
areass = []
legends = []
for iteration in iterations:
    areas = {'x': [], 'y': []}
    for tw in range(5, 1000, 5) + [2000, 3000, 4000, 5000]:
        if tw <= iteration:
            df = pd.read_hdf(basename+filename+"_"+str(tw)+".hdf", 'df')
            areas['x'].append(tw)
            areas['y'].append(float(df.y[df.x == iteration])/normalization)
    if areas['x']:
        areass.append(areas)
        legends.append(iteration)
ylim = (0.5,0.7)
ylim = (0.5,1.0)
Plotter.plot_curve(areass, figsize=(23,7), legends=legends, x_label="Time Window", markersize=7, linewidth=2, alpha=0.7, grid=True, ylim=ylim)

execfile("swarm_analyzer.py")
tw = 50
for tw in [125, 150, 175, 200, 225, 250]:
for tw in [5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 500, 750, 1000, 2000, 3000, 4000, 5000]:
for tw in range(5, 1000, 5) + [2000, 3000, 4000, 5000]:
    print tw
    runs = 1
    basename1 = "/mnt/50_particles_simulations/hdfs/pso_dynamic_initial_ring_F6_"
    basename2 = "/mnt/50_particles_simulations/hdfs/pso_global_F6_"
    basename3 = "/mnt/50_particles_simulations/hdfs/pso_ring_F6_"
    files1 = [basename1+"%02d_%d.hdf" % (i, tw) for i in range(1, runs+1)]
    files2 = [basename2+"%02d_%d.hdf" % (i, tw) for i in range(1, runs+1)]
    files3 = [basename3+"%02d_%d.hdf" % (i, tw) for i in range(1, runs+1)]
    files = files1
    files = files2
    files = files3
    files = files1 + files2 + files3
    c = [pd.read_hdf(f, 'df') for f in files]
    delta = 0.001
    number_particles = 50
    normalization = number_particles*len(np.arange(0, 1 + delta, delta))
    for d in c:
        d.y = d.y/normalization

    colors = ["#ff0000"]*len(files1) + ["#237530"]*len(files1) + ["#0000FF"]*len(files3)
    legends = ["Dynamic"] + [None]*(len(files1)-1) + ["Global"] + [None]*(len(files1)-1) + ["Ring"] + [None]*(len(files1)-1)
    markers = ["x"]*len(files1) + ["s"]*len(files1) + ["."]*len(files3)
    #Plotter.plot_curve(c, figsize=(18, 6), output_filename="all_topologies_"+str(tw)+".png", title="Time window size = "+str(tw), x_label="Iterations", legends=legends, linestyle=" ", markersize=4, linewidth=.2, marker=markers, markevery=20, colors=colors, ylim=(0,1), xlim=(0,6000), alpha=0.7)
    xticks = range(500, 6001, 500)
    rect = [0, 0, 0.93, 0.99]
    #Plotter.plot_curve(c, figsize=(18, 6), tight_layout=rect, grid=True, xticks_args=[xticks], output_filename="all_topologies_"+str(tw)+".png", title="Time window size = "+str(tw), x_label="Iterations", legends=legends, linestyle=" ", markersize=4, linewidth=.2, marker=markers, markevery=20, colors=colors, ylim=(0,1), xlim=(0,6000), alpha=0.7)
    Plotter.plot_curve(c, figsize=(18, 6), tight_layout=rect, grid=True, xticks_args=[xticks], output_filename="one_rune_all_topologies_"+str(tw)+".png", title="Time window size = "+str(tw), x_label="Iterations", legends=legends, linestyle=" ", markersize=4, linewidth=.2, marker=markers, colors=colors, ylim=(0,1), xlim=(0,6000), alpha=0.7)

    #Plotter.plot_curve(c, figsize=(18, 6),tight_layout=rect, output_filename="all_topologies_"+str(tw)+".png", title="Time window size = "+str(tw), x_label="Iterations", legends=legends, linestyle=" ", markersize=6, linewidth=.2, marker=markers, markevery=20, colors=colors, ylim=(0,1), alpha=0.9)
    #(left, bottom, right, top)
    #rect = [0, 0, 0.93, 0.99]
    #Plotter.plot_curve(c, tight_layout=rect, figsize=(15, 6), grid=True, xticks_args=[xticks], title="Time window size = "+str(tw), x_label="Iterations", legends=legends, linestyle=" ", markersize=4, linewidth=.2, marker=markers, markevery=20, colors=colors, ylim=(0,1), xlim=(0,6000), alpha=0.7)


Plotter.plot_curve(c, figsize=(18, 6), legends=range(1, runs+1), linestyle="-", markersize=3, linewidth=.2, marker='.', markevery=20, ylim=(0,1), xlim=(0,6000))


execfile("swarm_analyzer.py")
tws = [5, 10, 25, 50, 75, 100, 500, 1000, 2000, 3000]
tws = [5, 10, 25, 50, 75, 100]
tws = [5, 10, 25, 50, 100, 500, 1000, 2000, 3000]
tws = [ 500, 1000, 2000, 3000]
tws = [1]
tws = [10]
tws = range(5, 1000, 5)
run = "/mnt/50_particles_simulations/hdfs/pso_dynamic_initial_ring_F6_01"
run = "/mnt/50_particles_simulations/hdfs/pso_global_F6_02"
run = "/mnt/50_particles_simulations/hdfs/pso_ring_F6_01"
files = [run+"_%d.hdf" % (tw) for tw in tws]
c = [pd.read_hdf(f, 'df') for f in files]
delta = 0.001
number_particles = 50
normalization = number_particles*len(np.arange(0, 1 + delta, delta))
for d in c:
    d.y = d.y/normalization

Plotter.plot_curve(c, figsize=(20, 6), legends=tws, linestyle="-", markersize=5, linewidth=.2, markevery=20, ylim=(0,1))

import matplotlib.pyplot as plt
cm = plt.get_cmap('seismic')
colors = [cm(1.*i/len(c)) for i in range(len(c))]
for i in range(0, len(tws)-1):
    Plotter.plot_curve(c[i:], figsize=(20, 6), legends=tws[i:], output_filename="spso_dynamic_"+str(i)+".png", xlim=(0,6000), linestyle=" ", markersize=5, marker='o', linewidth=.2, markevery=10, ylim=(0,1), colors=colors[i:], alpha=0.5, markeredgewidth=.5)


         markeredgecolor='red', markeredgewidth=5)
colors = sns.color_palette("BuGn_r", len(c))
colors.

execfile("swarm_analyzer.py")
df = SwarmAnalyzer.get_giant_component_destruction_area(('None', "/mnt/50_particles_simulations/pso_global_F6_02"), 1)
df[0]


execfile("swarm_analyzer.py")
import matplotlib.pyplot as plt
import scipy.stats
delta = 0.001
number_particles = 50
normalization = number_particles*len(np.arange(0, 1 + delta, delta))
curves = []
run = 20
basename1 = "/mnt/50_particles_simulations/hdfs/pso_dynamic_initial_ring_F6_"
basename2 = "/mnt/50_particles_simulations/hdfs/pso_global_F6_"
basename3 = "/mnt/50_particles_simulations/hdfs/pso_ring_F6_"
basenames = [basename1]
basenames = [basename1, basename2, basename3]
# basenames = [basename3]
start = 0
max_y = float('-inf')
min_y = float('inf')
before = 1
times_before = 1
for basename in basenames:
    for tw in range(5, 1000, 5):
        filename = basename+"%02d_%d.hdf" % (run, tw)
        df = pd.read_hdf(filename, 'df')
        df = df[df['x'] > start]
        df['y'] /= normalization
        for _ in range(times_before):
            df['y'] = abs(df['y'] - pd.DataFrame({'y': list(df[before:]['y'])+[None]})['y'])
            df = pd.DataFrame({'x': df['x'][before:], 'y': df[:len(df)-before]['y']})
        max_y = max(max_y, max(df['y']))
        min_y = min(min_y, min(df['y']))
for basename in basenames:
    num_bins = 100
    ys, xs = [], []
    for tw in range(5, 1000, 5):
        filename = basename+"%02d_%d.hdf" % (run, tw)
        df = pd.read_hdf(filename, 'df')
        df = df[df['x'] > start]
        df['y'] /= normalization
        for _ in range(times_before):
            #df['y'] = abs(df['y'] - pd.DataFrame({'y': list(df[before:]['y'])+[None]})['y'])
            df['y'] = df['y'] - pd.DataFrame({'y': list(df[before:]['y'])+[None]})['y']
            df = pd.DataFrame({'x': df['x'][before:], 'y': df[:len(df)-before]['y']})

        counts, bins = np.histogram(df['y'], bins=num_bins, range=(min_y, max_y))
        #counts, bins = np.histogram(df['y'], bins=num_bins) #, range=(0, 1.0))
        bins = bins[:-1] + (bins[1] - bins[0])/2
        probs = counts/float(counts.sum())
        # print probs.sum() # 1.0
        # plt.bar(bins, probs, 1.0/num_bins)
        # plt.show()
        ys.append(scipy.stats.entropy(probs))
        xs.append(tw)
    curves.append({'x': xs, 'y': ys})
Plotter.plot_curve(curves, figsize=(23,7), legends=['d', 'g', 'r'], x_label="Time Window", markersize=7, linewidth=2, alpha=0.7, grid=True, x_scale='log', y_scale='linear')
    """
