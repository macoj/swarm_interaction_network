import time
import igraph

__author__ = 'marcos'
import sys
import pandas as pd
import numpy as np
from plotter import Plotter
from scipy import interpolate
from callbacks import Callback
from swarm_parser import SwarmParser
from giant_component_analysis import GiantComponentDeath
from giant_component_analysis_plotter import GiantComponentDeathPlotter
from scipy.stats import pearsonr


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
        weight_normalize = [2.0*i if i < window_size else 2.0 * window_size for i in range(len(graph_matrices))]
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
            normalization = float(normalize_max_y)*len(tx)
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
            areas = [a/normalization for a in areas]
        return areas

    @staticmethod
    def get_giant_component_destruction_area(filename, window_size, number_of_individuals=100, until=-1):
        graphs = SwarmAnalyzer.get_giant_component_destruction_curves(filename, window_size=window_size, until=until)
        areas = SwarmAnalyzer.get_areas_under_curves(graphs, normalize=True, normalize_max_y=number_of_individuals)
        #df = pd.DataFrame({'x': range(window_size, window_size + len(areas)), 'y': areas})
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
            filename, dimensions=1000, particles=100, informations_grep="velocities\:#", kind='all', absolute=False,
            **kargs):
        information_map = lambda x: np.array(map(float, x.split())).reshape(particles, dimensions)
        velocities = SwarmParser.read_file_and_measures(
            filename, influence_graph_grep=None,
            informations_grep=informations_grep, information_map=information_map, **kargs)
        print 1
        velocities = [v[1] for v in velocities[1][informations_grep]]
        print 2
        velocities = [v/np.linalg.norm(v) for v in velocities]
        print 3
        if absolute:
            abs_cmp = lambda x, y: cmp(abs(x), abs(y))
        else:
            abs_cmp = cmp
        correlation_t = []
        v_sum = []
        for iteration in range(len(velocities)):
            if kind in ['all']:
                correlations = pd.DataFrame(np.rot90(velocities[iteration])).corr()
                # if kind == 'average_highest':
                #     correlation = np.average([sorted(correlations[p], cmp=abs_cmp, reverse=True)[1]
                #                               for p in range(particles)])
                # elif kind == 'average_average':
                #     correlation = np.average([np.average(sorted(correlations[p], cmp=abs_cmp, reverse=True)[1:])
                #                               for p in range(particles)])
                # else:
                correlation = np.array(correlations).reshape(1, correlations.shape[0]*correlations.shape[1])[0]
                correlation_t.append(correlation)
            elif kind in ['average', 'fluctuations']:
                if iteration == 0:
                    v_sum = velocities[0]
                else:
                    v_sum += velocities[iteration]
                v_average = v_sum / float(iteration+1)
                correlations = None
                if kind == 'fluctuations':
                    if iteration < len(velocities):
                        fluctuations = [velocities[iteration+1][p] - v_average[p] for p in range(particles)]
                        correlations = pd.DataFrame(np.rot90(fluctuations)).corr()
                else:
                    correlations = pd.DataFrame(np.rot90(v_average)).corr()
                if correlations is not None:
                    correlation = np.array(correlations).reshape(1, correlations.shape[0]*correlations.shape[1])[0]
                    correlation_t.append(correlation)
                print 11
        return correlation_t
# a = velocities[0][0]
# b = velocities[0][1]
# np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))
# pearsonr(a, b)[0]
#
# aa = a - np.average(a)
# bb = b - np.average(b)
# np.dot(aa, bb)/(np.linalg.norm(aa)*np.linalg.norm(bb))
    """
filename = "./data/100_particles/global_F21_00.with_positions"
filename = "./data/100_particles/regular30_F21_00.with_positions_head_30"
correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind='average')

execfile("swarm_analyzer.py")
info = "it\:#"

t = SwarmParser.read_file_and_measures(filename, informations_grep=[info])

curves = zip(*t[1][info])
curves = {'x': curves[0], 'y': curves[1]}
Plotter.plot_curve(curves)

import matplotlib
matplotlib.use('Agg')
execfile("swarm_analyzer.py")
for topology in ['global', 'regular30', 'ring']:
    for function in range(21, 28):
        for run in [0, 1]:
            filename = "./%s_F%02d_%02d.with_positions" % (topology, function, run)
            print filename
            kind = "fluctuations"
            print filename
            correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind=kind, until=1000)
            df = pd.DataFrame(correlation_t)
            df.to_hdf(filename + "_fluctuations_correlation.hdf", 'df')


kind = "highest"
kind = "average"
correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind=kind, until=until)
for particle in range(100):
    Plotter.plot_curve({'x': range(len(correlation_t[particle])), 'y': correlation_t[particle]}, dpi=72, figsize=(20, 5), tight_layout=[], x_label="Iteration", y_label="%s pearson correlation" % kind, title="Particle #%d" % particle, output_filename="%s_particle_%d_%s_F%02d_%02d.png" % (kind, particle, topology, function, run))

kind = "average_highest"
correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind=kind, until=20)
Plotter.plot_curve({'x': range(len(correlation_t)), 'y': correlation_t}, dpi=72, figsize=(20, 5), tight_layout=[], x_label="Iteration", y_label="%s pearson correlation" % kind, title="All particles", output_filename="%s_particle_%s_F%02d_%02d.png" % (kind, topology, function, run))


kind = "average_average"
correlation_t = SwarmAnalyzer.calculate_velocities_correlation(filename, kind=kind, until=20)
Plotter.plot_curve({'x': range(len(correlation_t)), 'y': correlation_t}, dpi=72, figsize=(20, 5), tight_layout=[], x_label="Iteration", y_label="%s pearson correlation" % kind, title="All particles", output_filename="%s_particle_%s_F%02d_%02d.png" % (kind, topology, function, run))

execfile("plotter.py")
import matplotlib
matplotlib.use('Agg')
execfile("swarm_analyzer.py")
names = {21: "Ackley",  22: "Griewank", 23: "Rastrigin", 24: "Rosenbrock", 25: "Schwefel", 26: "Sphere", 27: "Weierstrass"}
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
for topology in ['global', 'regular30', 'ring']:
    for function in range(21, 28):
        for run in [0, 1]:
            filename = "./%s_F%02d_%02d.with_positions_correlation" % (topology, function, run)
            df = pd.read_hdf(filename + ".hdf", 'df')
            bins = np.arange(-1, 1.01, 0.01)
            iterations = 6000
            counts = []
            for i in range(iterations):
                count, _ = np.histogram(df.irow(i), bins=bins)
                count = count/float(max(count))
                counts.append(count)
            title = "%s - f:'%s' run:%d" % (topology, names[function], run)
            Plotter.plot_heatmap(np.rot90(counts), main_title=title, output_filename=filename + ".png", vmax=None, vmin=None, set_yticks=[0, len(bins)/2.0, len(bins)-1], titles_y=["-1", "0", "1"], tight_layout=[0, 0, 1, 0.98], figsize=(30, 4), colorbar_on=False)



title = ""
xlabel = "Average communication diversity"
ylabel = "Stagnation iteration"
Plotter.plot_curve(dss, legend_ncol=2, ylim=(200, 1450), xlim=(0.05, 0.376), output_filename=output_filename, linewidth=0, legends=legends, mew=1.2, figsize=(4.8, 3.05), loc=2, tight_layout=[-0.03, -0.03, 1.03, 1.03], annotate=(0.66, 0.86, "$r=%0.2f$" % p[0]), font_size=11, alpha=1, marker=['+', 'x', '>', '^', 'o', 'v', 'D', 's', '3', '2', '<'], colors=colors, grid=True, markersize=4, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))



np.rotate90(counts)
np.array(counts)

Plotter.plot_
    """

    @staticmethod
    def get_number_of_components(filename, window_size, min_weight, **kargs):
        influence_graph_grep = 'ig\:#'
        pos_callback = lambda x: SwarmAnalyzer.get_number_of_components_of_graph(x, min_weight=min_weight*2*window_size,
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
    def generate_swarm_analysis_df_from_file(filename, informations_grep=None, tws=None, normalize=False):
        if tws is None:
            tws = [10, 25, 50, 75, 100, 200, 300, 400, 500, 1000]
            # tws = [100, 200, 300]
        if informations_grep is None:
            informations_grep = ["radius:#", "diameter:#", "average_around_center:#", "it:#", "coherence:#",
                                 "normalized_average_around_center:#", "aggregation_factor:#",
                                 "average_of_average_around_all_particles:#"]
        print "Generating swarm analysis dataframe from '%s'" % filename
        print " > information parsing: %s" % str(informations_grep)
        df_info = SwarmAnalyzer.get_swarm_informations_from_file(filename, informations_grep)
        if normalize:
            df_columns = {k: dict(zip(df_info['x'], df_info[k]/max(df_info[k]))) for k in informations_grep}
        else:
            df_columns = {k: dict(zip(df_info['x'], df_info[k])) for k in informations_grep}
        print "  > OK"
        print " > calculating AUC for tw: %s" % str(tws)
        for tw in tws:
            print "  > tw: %d (%s)" % (tw, time.strftime("%H:%M:%S-%d/%m/%y"))
            auc = SwarmAnalyzer.get_giant_component_destruction_area(filename, window_size=tw)
            print "   > OK"
            df_columns[tw] = dict(zip(auc['x'], auc['y']))
        print "  > OK"
        print " > creating data frame"
        keys = []
        for k in df_columns:
            keys += df_columns[k].keys()
        keys = list(set(keys))
        keys.sort()
        df = {'x': keys}
        for k in df_columns:
            df[k] = [df_columns[k][x] if x in df_columns[k] else None for x in keys]
        real_df = pd.DataFrame(df)
        real_df.index = real_df['x']
        del real_df['x']
        return real_df
    """
execfile("swarm_analyzer.py")
import time
filename = "./data/100_particles/ring_F06_06.teste"
df = SwarmAnalyzer.generate_swarm_analysis_df_from_file(filename)
    """

    @staticmethod
    def create_swarm_analysis_df_file_from_file(filename, output_hdf, informations_grep=None, tws=None,
                                                normalize=False):
        df = SwarmAnalyzer.generate_swarm_analysis_df_from_file(filename, informations_grep, tws, normalize)
        df.to_hdf(output_hdf, 'df')
    """
10 20 30 40 50 60 70 80 90 100
10 3
50 3
70 3

1 2 3
4 5 6 7 8
9 10 11 12 13
14 15 16 17 18
19 20
execfile("swarm_analyzer.py")
output_base_name = './data/'
topologies = ["ring", "vonneumann"] + ["kregular%d" % d for d in range(5, 10) + range(10, 100, 10)] + ["global"]
function = 2
0, 4
4, 9
9, 13
13, 17
for topology in topologies[13:17]:
    if topology in ["ring", "vonneumann", "global"]:
        suffix = ".teste"
    else:
        suffix = ""
    filenames = ['%s_F%02d_%02d%s' % (topology, function, run, suffix) for run in range(30)]
    for filename in filenames:
        print time.strftime(">> %H:%M:%S-%d/%m/%y")
        print filename, output_base_name+filename+'.hdf'
        df = SwarmAnalyzer.create_swarm_analysis_df_file_from_file(filename, output_base_name+filename+'.hdf')


    """

    @staticmethod
    def retrieve_stagnation_iteration_derivative(
            filename, smooth=50, diff_threshold=0.02, iter_threshold=500, head=None):
        df = pd.read_hdf(filename, 'df')
        if head:
            df = df.head(head)
        x = df.index[smooth:len(df)]
        fitness = np.array(df['it:#'][smooth:])
        # f_delta = [f(t) - f(t+1)]/f(t+1)
        # f_delta = f(t)/f(t+1) - 1
        fitness_diff = fitness[:len(fitness)-1] / fitness[1:len(fitness)]
        fitness_diff -= 1

        # noinspection PyTypeChecker
        fitness_diff = list(fitness_diff > diff_threshold)
        counts = []
        last = fitness_diff[0]
        last_c = 1
        for c in fitness_diff[1:]:
            if c == last:
                last_c += 1
            else:
                counts.append(last_c)
                last_c = 1
            last = c
        if last_c != 1:
            counts.append(last_c)
        threshold_i = 0
        for threshold_i in range(len(counts)):
            if counts[threshold_i] > iter_threshold:
                break
        stagnation_iteration = sum(counts[:threshold_i-1])
        return stagnation_iteration
    """
    execfile("swarm_analyzer.py")
    execfile("plotter.py")
    it = "it:#"
    topologies = ['ring', 'vonneumann', 'kregular5', 'kregular6', 'kregular7', 'kregular8', 'kregular9',
                  'kregular10', 'kregular20', 'kregular30', 'kregular40', 'kregular50',
                  'kregular60', 'kregular70', 'kregular80', 'kregular90', 'global']
    function = 14
    for t in topologies:
        for i in range(5):
            filename = "./data/100_particles/%s_F%02d_%02d.hdf" % (t, function, r)
            diff_threshold, iter_threshold = 0.00001, 500
            stagnation_iteration = SwarmAnalyzer.retrieve_stagnation_iteration_derivative(filename, 10, diff_threshold=diff_threshold, iter_threshold=iter_threshold)
            df = pd.read_hdf(filename, 'df')
            Plotter.plot_curve({'x': df.index, 'y': df[it]}, grid=True, vline_at=stagnation_iteration, figsize=(10, 4), xlim=(0, 2000))

    """

    # noinspection PyTypeChecker
    @staticmethod
    def calculate_cd_and_dfit(function=14, mean=True, tws=None, topologies=None, runs=30,
                              diff_threshold=10**-5, iter_threshold=500, fitness_only=False,
                              max_iteration_mode=None):
        it = "it:#"
        if not topologies:
            topologies = ['ring', 'vonneumann', 'kregular5', 'kregular6', 'kregular7', 'kregular8', 'kregular9',
                          'kregular10', 'kregular20', 'kregular30', 'kregular40', 'kregular50',
                          'kregular60', 'kregular70', 'kregular80', 'kregular90', 'global']
        if not tws:
            tws = [10, 25, 50, 75, 100, 200, 300, 400, 500, 1000]
        dfits = []
        fitnesses = []
        cds = []
        for t in topologies:
            fitness_t = []
            for r in range(runs):
                filename = "./data/100_particles/%s_F%02d_%02d.hdf" % (t, function, r)
                if max_iteration_mode is None:
                    stagnation_iteration = SwarmAnalyzer.retrieve_stagnation_iteration_derivative(
                        filename, diff_threshold=diff_threshold, iter_threshold=iter_threshold)
                else:
                    stagnation_iteration = max_iteration_mode
                df = pd.read_hdf(filename, 'df')
                fitness = np.array(df[it])
                if fitness_only:
                    fitness_t.append(fitness[-1])
                else:
                    # f_delta = [f(t) - f(t+1)]/f(t+1) = f(t)/f(t+1) - 1
                    fitness_diff = fitness[:len(fitness) - 1] / fitness[1:len(fitness)] - 1
                    cd = 1 - sum([df[c] for c in tws]) / float(len(tws))
                    if mean:
                        cds.append(np.mean(list(cd[1:stagnation_iteration + 1])))
                        if not fitness_only:
                            dfits.append(np.mean(list(fitness_diff[:stagnation_iteration])))
                    else:
                        cds += list(cd[1:stagnation_iteration + 1])
                        dfits += list(fitness_diff[:stagnation_iteration])
                del df
            if mean:
                fitnesses.append(np.mean(fitness_t))
            else:
                fitnesses.append(fitness_t)
        if fitness_only:
            return fitnesses
        else:
            return cds, dfits
    """
execfile("swarm_analyzer.py")
execfile("plotter.py")
function = 2
cds, dfits = SwarmAnalyzer.calculate_cd_and_dfit(function=function, runs=30, diff_threshold=10**-4, mean=True)
indices = range(0, len(cds), 30)
curves = []
bb = []
for i,j in zip(ds_index, ds_index[1:]):
    curves.append({'x': cds[i:j], 'y': dfits[i:j]})
    bb.append(dfits[i:j])

pearson_correlation = pearsonr(cds, dfits)
print pearson_correlation

topologies = ['Ring', 'Neumann'] + ['$k=%d$' % k for k in [5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90]] + ["Global"]
colorblind = [(0, 0, 0), (230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66), (0, 114, 178), (213, 94, 0), (204, 121, 167)]
colorblind = ['#%02x%02x%02x' % (c[0], c[1], c[2]) for c in colorblind]
colors = colorblind[3], colorblind[5], colorblind[6], colorblind[1], colorblind[2], colorblind[7]
colors = colorblind[1:]
Plotter.plos_style()
Plotter.plot_curve(curves, output_filename="dfits_cds_f_%02d.pdf" % function, legends=topologies, annotate=(0.78, 0.86, "$r=%0.2f$" % pearson_correlation[0]), tight_layout=[-0.025, -0.035, 1.02, 1.03], legend_ncol=3, colors=colors, linestyle="", figsize=(4.7, 2.7), loc=3, x_label="Average communication diversity", mew=1.2, y_label="Average fitness improvement", marker=['+', 'x', '>', '^', 'o', 'v', 'D', 's', '3', '2', '<'], grid=True, markersize=4)



execfile("swarm_analyzer.py")
execfile("plotter.py")
function = 2
bb = SwarmAnalyzer.calculate_cd_and_dfit(function=function, runs=30, diff_threshold=10**-4, fitness_only=True)
indices = range(0, len(fitness), 30)
bb = []
for i,j in zip(ds_index, ds_index[1:]):
    bb.append(dfits[i:j])


boxes_kargs = {'color': 'black', 'linewidth': 1.3, 'zorder': 3, 'fillstyle': 'full', 'facecolor': '#a6cee3'}
means_kargs = {'color': 'black', 'fillstyle': 'full', 'markerfacecolor': "black", 'marker': "s", 'markersize': 3, 'mew': 1, 'mec': 'black', 'zorder': 5}
fliers_kargs = {'color': 'black', 'marker': ".", 'markersize': 3, 'mew': 1.5, 'mec': 'black'}
whiskers_kargs = {'color': 'black', 'linewidth': 1.2, 'zorder': 2, 'linestyle': '-', 'alpha': 1.0, 'mec': 'black'}
medians_kargs = {'color': 'black', 'linewidth': 1.6, 'zorder': 5, 'alpha': 0.3}
caps_kargs = {'linewidth': 1.5, 'color': 'black',}

topologies = [2, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
topologies_legend = map(str, topologies)
topologies_legend.reverse()
# topologies_legend[0] = 'Ring'
# topologies_legend[len(topologies)-1] = 'Global'
legends = [range(1, len(cds)/30 + 1), topologies_legend, {"rotation": 0}]
Plotter.plot_boxplots(bb, output="fitness_f_%02d_r.pdf" % function, ylim=(6000, 11500), xticks_args=legends, grid_only='y', whis=1.5, xlim=(0.3, len(topologies)+0.6), grid=True, widths=0.7, tight_layout=[-0.03, -0.03, 1.025, 1.03], bootstrap=2000, boxes_kargs=boxes_kargs, showmeans=False, ylabel="Fitness", xlabel="$k$-regular topologies", size=(4.5, 2.8), showfliers=True, fliers_kargs=fliers_kargs, means_kargs=means_kargs, whiskers_kargs=whiskers_kargs, medians_kargs=medians_kargs, caps_kargs=caps_kargs)
plt.clf()
Plotter.plot_boxplots(bb)




for function in [2, 6, 14, 19]:
    print function
    fitness = SwarmAnalyzer.calculate_cd_and_dfit(function=function, runs=30, diff_threshold=10**-4, fitness_only=True)
    jump_every = 6
    jump_every_i = 0
    tex = ""
    for fitness_ in fitness:
        f_fitness = re.sub("e\+[0]", "\\\\times10^{", "{0:.4e}".format(fitness_))
        tex += "& $%s}$  " % f_fitness
        jump_every_i += 1
        if jump_every_i % jump_every == 0:
            print tex + "\\\\"
            tex = ""
    print tex + "\\\\"



#output_filename="f_%02d.png" % function
title = ""
xlabel = "Average communication diversity"
ylabel = "Stagnation iteration"
Plotter.plot_curve(dss, legend_ncol=2, ylim=(200, 1450), xlim=(0.05, 0.376), output_filename=output_filename, linewidth=0, legends=legends, mew=1.2, figsize=(4.8, 3.05), loc=2, tight_layout=[-0.03, -0.03, 1.03, 1.03], annotate=(0.66, 0.86, "$r=%0.2f$" % p[0]), font_size=11, alpha=1, marker=['+', 'x', '>', '^', 'o', 'v', 'D', 's', '3', '2', '<'], colors=colors, grid=True, markersize=4, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))



plt.clf()
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

weights = np.ones_like(dfits)/len(dfits)
fig = plt.figure(figsize=(4.5, 2.6))
plt.yscale('log')
plt.hist2d(cds, dfits, bins=100, cmap=plt.cm.gist_heat_r, weights=weights)
plt.show()

    """

    @staticmethod
    def calculate_correlation_cd_dfit(**kargs):
        cds, dfits = SwarmAnalyzer.calculate_cd_and_dfit(**kargs)
        return pearsonr(cds, dfits)
    """
execfile("swarm_analyzer.py")
SwarmAnalyzer.calculate_correlation_cd_dfit(function=6, runs=10)


>>> SwarmAnalyzer.calculate_correlation_cd_dfit(function=6)
(-0.92356894964013847, 1.2505047550893154e-213)
>>> SwarmAnalyzer.calculate_correlation_cd_dfit(function=14)
(-0.33756838240407916, 4.6573148666352012e-15)
>>> SwarmAnalyzer.calculate_correlation_cd_dfit(function=19)
(-0.34651911636639665, 7.8013584786842649e-16)
>>> SwarmAnalyzer.calculate_correlation_cd_dfit(function=2)
(-0.78771894359497063, 5.9069928507128958e-109)

    """

    @staticmethod
    def analyze_analysis_df_file_stagnation_derivative(filename, tw, smooth=50, diff_threshold=0.02,
                                                       iter_threshold=500, head=None, plot=False, tws=None):
        df = pd.read_hdf(filename, 'df')
        print filename
        if head:
            df = df.head(head)
        x = df.index[smooth:len(df)]
        if not tws:
            communication_diversity_smoothed = [sum(df[tw][n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
        else:
            sum_tws = sum([df[c] for c in tws])/float(len(tws))
            communication_diversity_smoothed = [sum(sum_tws[n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
        communication_diversity_smoothed = 1 - np.array(communication_diversity_smoothed)
        fitness = np.array(df['it:#'][smooth:])
        #
        # f_delta = [f(t) - f(t+1)]/f(t+1)
        # f_delta = f(t)/f(t+1) - 1
        #
        fitness_diff = fitness[:len(fitness)-1] / fitness[1:len(fitness)]
        # Plotter.plot_curve({'x': range(len(fitness_diff)), 'y': fitness_diff})
        # Plotter.plot_curve({'x': range(len(fitness)), 'y': fitness}, vline_at=stagnation_iteration)

        fitness_diff -= 1
        # noinspection PyTypeChecker
        fitness_diff = list(fitness_diff > diff_threshold)
        counts = []
        last = fitness_diff[0]
        last_c = 1
        for c in fitness_diff[1:]:
            if c == last:
                last_c += 1
            else:
                counts.append(last_c)
                last_c = 1
            last = c
        if last_c != 1:
            counts.append(last_c)
        threshold_i = 0
        for threshold_i in range(len(counts)):
            if counts[threshold_i] > iter_threshold:
                break
        stagnation_iteration = sum(counts[:threshold_i-1])
        communication_diversity_until_stagnation = communication_diversity_smoothed[:stagnation_iteration]
        if plot:
            Plotter.plot_curve({'x': x, 'y': fitness/max(fitness)}, vline_at=stagnation_iteration, figsize=(15, 4),
                               markersize=3, ylim=(0, 1.0), grid=True)

        return sum(communication_diversity_until_stagnation)/float(stagnation_iteration), float(stagnation_iteration)
    """
    execfile("swarm_analyzer.py")
    i = 0
    filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
    filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
    filename = "./data/100_particles/kregular6_F06_%02d.teste.hdf" % i
    sum_cd, stagnation_it = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative(filename, 100, head=3000, diff_threshold=0.01, plot=True)

    """

    @staticmethod
    def analyze_analysis_df_file_stagnation_derivative_and_measures(filename, tw, smooth=50, diff_threshold=0.0001,
                                                                    iter_threshold=500, head=None, tws=None, plot=False,
                                                                    **kargs):
        df = pd.read_hdf(filename, 'df')
        if head:
            df = df.head(head)
        x = df.index[smooth:len(df)]
        if not tws:
            communication_diversity_smoothed = [sum(df[tw][n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
        else:
            sum_tws = sum([df[c] for c in tws])
            communication_diversity_smoothed = [sum(sum_tws[n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
        curves = [{'x': x, 'y': communication_diversity_smoothed}]
        fitness = np.array(df['it:#'][smooth:])
        curves += [{'x': x, 'y': fitness/fitness.max()}]
        measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',
                    'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']

        fitness_diff = fitness[:len(fitness)-1] / fitness[1:len(fitness)]
        fitness_diff -= 1
        # noinspection PyTypeChecker
        fitness_diff = list(fitness_diff > diff_threshold)
        counts = []
        last = fitness_diff[0]
        last_c = 1
        for c in fitness_diff[1:]:
            if c == last:
                last_c += 1
            else:
                counts.append(last_c)
                last_c = 1
            last = c
        if last_c != 1:
            counts.append(last_c)
        threshold_i = 0
        for threshold_i in range(len(counts)):
            if counts[threshold_i] > iter_threshold:
                break
        stagnation_iteration = sum(counts[:threshold_i-1])
        for m in measures:
            ddx = df.index[smooth:len(df)]
            ddy = [sum(df[m][n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
            ddy /= max(ddy)
            curves += [{'x': ddx, 'y': ddy}]
        legends = ['cd', 'it']+measures
        if plot:
            Plotter.plot_curve(curves, legends=legends, vline_at=stagnation_iteration, figsize=(15, 4),
                               markersize=3, ylim=(0, 1.0), grid=True, loc=1, **kargs)
        else:
            result = [stagnation_iteration, np.mean(curves[0]['y'][:stagnation_iteration])]
            for c in curves[2:]:
                result.append(np.mean(c['y'][:stagnation_iteration]))
            return result

    """
    execfile("swarm_analyzer.py")
    results = []
    tw = 10
    topologies = ['global', 'ring', 'vonneumann']
    for t in topologies:
        for i in range(30):
            # filename = "./data/100_particles/dynamicring_F06_%02d.teste.hdf" % i
            # filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
            # filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
            # filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
            filename = "./data/100_particles/%s_F06_%02d.teste.hdf" % (t, i)
            # output_filename="auc_smoothed_ring_%03d_%02d.png" % (tw, i)
            # output_filename="auc_smoothed_vonneumann_%03d_%02d.png" % (tw, i)
            # output_filename="auc_smoothed_global_%03d_%02d.png" % (tw, i)
            # output_filename="auc_smoothed_dynamicring_%03d_%02d.png" % (tw, i)
            # SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_and_measures(filename, tw, head=3000, diff_threshold=0.02, output_filename=output_filename, smooth=1)
            # SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_and_measures(filename, tw, head=3000, diff_threshold=0.02, smooth=1, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
            results.append(SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_and_measures(filename, tw, head=3000, diff_threshold=0.02, smooth=1, plot=False, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])))
    # SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_and_measures(filename, tw, head=3000, diff_threshold=0.02, smooth=1, plot=True)
    df = pd.DataFrame(np.array(results))
    runs = 30
    ds_index = range(0, len(topologies)*runs + 1, runs)
    for measure in range(1, 7):
        ds = []
        # measure = 3
        away_std = 1.5
        for i,j in zip(ds_index, ds_index[1:]):
            d = pd.DataFrame({'x': df[measure][i:j], 'y': df[0][i:j]})
            dd = d[d.y < np.mean(d.y) + away_std*np.std(d.y)]
            dd = dd[dd.y > np.mean(dd.y) - away_std*np.std(dd.y)]
            dd = dd[dd.x < np.mean(dd.x) + away_std*np.std(dd.x)]
            dd = dd[dd.x > np.mean(dd.x) - away_std*np.std(dd.x)]
            ds.append(dd)
        ddd = ds[0]
        for d in ds[1:]:
            ddd = ddd.append(d)
        legends = topologies
        print pearsonr(ddd.x, ddd.y)
        Plotter.plot_curve(ds, linewidth=0, legends=legends, title=str(measure), figsize=(5, 3.8), loc=2, tight_layout=[], annotate=(0.49, 0.85, "$r=%0.2f$" % (ddd.corr()['x']['y'])), font_size=11, alpha=0.6, marker='.', colors=["#e41a1c", "#4daf4a", "#377eb8", "#fdae61"], grid=True, markersize=6)

        ddd = ds[1]
        Plotter.plot_curve(ddd, linewidth=0, legends=legends, title=str(measure), figsize=(5, 3.8), loc=2, tight_layout=[], annotate=(0.49, 0.85, "$r=%0.2f$" % (ddd.corr()['x']['y'])), font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6)

from scipy.stats import pearsonr

    """

    @staticmethod
    def analyze_analysis_df_file_measures_and_cd(filename, tw, iteration=None, **kargs):

        measures_and_cd = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_measures_and_cd(filename, tw,
                                                                                                       mean=False,
                                                                                                       **kargs)
        for k in measures_and_cd[0]:
            measures_and_cd[0][k] = measures_and_cd[0][k][:-1]
        df = pd.DataFrame(measures_and_cd[0])
        for m in measures_and_cd[1:]:
            ddy = m[m.keys()[0]]
            ddy = np.array(ddy[1:])/np.array(ddy[:-1])
            df[m.keys()[0]] = ddy
        if iteration:
            df = df.head(iteration + 1)
            df = df.tail(1)
        return df
    """
    execfile("swarm_analyzer.py")
    #for iteration in range(1, 50, 10):
    iteration = None
    results = []
    topologies = ['global', 'ring', 'vonneumann']
    for t in topologies:
        for i in range(30):
            filename = "./data/100_particles/%s_F06_%02d.teste.hdf" % (t, i)
            print filename
            result = SwarmAnalyzer.analyze_analysis_df_file_measures_and_cd(filename, tw=None, iteration=iteration, head=3000, diff_threshold=0.01, smooth=25, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
            results.append(result)

    df = results[0]
    for d in results[1:]:
        df = df.append(d)
    del df['x']

    print list(df.corr()['cd'])

    measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#', 'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']
    for m in measures:
        Plotter.plot_curve({'x': df.cd, 'y': df[m]}, linewidth=0, figsize=(5, 3.8), loc=2, tight_layout=[], font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6)

    by_topologies = []
    for r in range(0, 90, 30):
        res = results[r:r+30]
        dd = res[0]
        for d in res[1:]:
            dd = dd.append(d)
        by_topologies.append(dd)

    for ddd in by_topologies:
        for m in measures:
            Plotter.plot_curve({'x': ddd.cd, 'y': ddd[m]}, linewidth=0, figsize=(5, 3.8), loc=2, tight_layout=[], font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6)

    """

    @staticmethod
    def analyze_analysis_df_file_stagnation_derivative_measures_and_cd(filename, tw, smooth=50, diff_threshold=0.0001,
                                                                       iter_threshold=500, head=None, tws=None,
                                                                       mean=True, no_measure=False):
        df = pd.read_hdf(filename, 'df')
        if head:
            df = df.head(head)
        x = df.index[smooth:len(df)]
        if not tws:
            if smooth != 0:
                communication_diversity_smoothed = [sum(df[tw][n-smooth:n])/float(smooth)
                                                    for n in range(smooth, len(df))]
            else:
                communication_diversity_smoothed = df[tw]
        else:
            sum_tws = sum([df[c] for c in tws])/float(len(tws))
            if smooth != 0:
                communication_diversity_smoothed = [sum(sum_tws[n-smooth:n])/float(smooth)
                                                    for n in range(smooth, len(df))]
            else:
                communication_diversity_smoothed = sum_tws
        # changing the definition of communication diversity
        communication_diversity_smoothed = list(1.0 - np.array(communication_diversity_smoothed))
        fitness = np.array(df['it:#'][smooth:])
        measures = []
        if not no_measure:
            measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',
                        'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#',
                        'radius:#']
        # derivative
        fitness_diff = fitness[:len(fitness)-1] / fitness[1:len(fitness)]
        fitness_diff -= 1
        # noinspection PyTypeChecker
        fitness_diff = list(fitness_diff > diff_threshold)
        counts = []
        last = fitness_diff[0]
        last_c = 1
        for c in fitness_diff[1:]:
            if c == last:
                last_c += 1
            else:
                counts.append(last_c)
                last_c = 1
            last = c
        if last_c != 1:
            counts.append(last_c)
        threshold_i = 0
        for threshold_i in range(len(counts)):
            if counts[threshold_i] > iter_threshold:
                break
        stagnation_iteration = sum(counts[:threshold_i-1])
        if not mean:
            dfit = (fitness[:len(fitness)-1] - fitness[1:len(fitness)]) / fitness[1:len(fitness)]
            if smooth != 0:
                dfit = [sum(dfit[n-smooth:n])/float(smooth) for n in range(smooth, len(dfit))]
            curves = [{'x': x[:stagnation_iteration-smooth],
                       'cd': communication_diversity_smoothed[:stagnation_iteration-smooth],
                       'dfit': dfit[:stagnation_iteration-smooth]}]
        else:
            #curves = [{'cd': np.mean(communication_diversity_smoothed[:stagnation_iteration])}]
            curves = [{'cd': communication_diversity_smoothed[stagnation_iteration-smooth]}]
        for m in measures:
            if smooth != 0:
                ddy = [sum(df[m][n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
            else:
                ddy = df[m]
            # ddy /= max(ddy)
            #ddy = (np.array(ddy[1:]) - np.array(ddy[:-1])) / np.array(ddy[:-1])
            #ddy = np.array(ddy[1:]) / np.array(ddy[:-1])
            if not mean:
                curves += [{m: ddy[:stagnation_iteration-smooth]}]
            else:
                # curves += [{m: np.mean(ddy[:stagnation_iteration])}]
                curves += [{m: ddy[stagnation_iteration-smooth]}]
        return curves

    """

execfile("swarm_analyzer.py")
from scipy.stats import pearsonr

topologies = ['ring', 'vonneumann', 'kregular5', 'kregular6', 'kregular7', 'kregular8', 'kregular9',
              'kregular10', 'kregular20', 'kregular30', 'kregular40', 'kregular50',
              'kregular60', 'kregular70', 'kregular80', 'kregular90', 'global']
for f in [2, 6, 14, 19]:
    cds = []
    dfits = []
    tws = [10, 25, 50, 75, 100, 200, 300, 400, 500, 1000]
    for t in topologies:
        for i in range(20):
            filename = "./data/100_particles/%s_F%02d_%02d.hdf" % (t, f, i)
            df = pd.read_hdf(filename, 'df')
            cd = 1 - sum([df[c] for c in tws])/float(len(tws))
            fitness = np.array(df['it:#'])

            fitness_diff = fitness[:len(fitness)-1] / fitness[1:len(fitness)]
            fitness_diff -= 1
            max_it = 1000
            # cds += list(cd[1:max_it+1])
            # dfits += list(fitness_diff[:max_it])
            cds.append(np.mean(list(cd[1:max_it+1])))
            dfits.append(np.mean(list(fitness_diff[:max_it])))
    print f, max_it, pearsonr(cds, dfits)

Plotter.plot_curve({'x': cds, 'y': dfits}, linestyle='')

execfile("swarm_analyzer.py")
from scipy.stats import pearsonr

topologies = ['global']
topologies = ['ring']
topologies = ['vonneumann']
topologies = ['global', 'ring', 'vonneumann']
topologies = ['ring', 'kregular4', 'kregular5', 'kregular6', 'kregular7', 'kregular8', 'kregular9',
              'kregular10', 'kregular20', 'kregular30', 'kregular40', 'kregular50',
              'kregular60', 'kregular70', 'kregular80', 'kregular90', 'global']
topologies = ['ring', 'vonneumann', 'kregular5', 'kregular6', 'kregular7', 'kregular8', 'kregular9',
              'kregular10', 'kregular20', 'kregular30', 'kregular40', 'kregular50',
              'kregular60', 'kregular70', 'kregular80', 'kregular90', 'global']
f = 14
results = []
results_ = {}
for t in topologies:
    for i in range(5):
        filename = "./data/100_particles/%s_F%02d_%02d.hdf" % (t, f, i)
        print filename
        mean = False
        result = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_measures_and_cd(filename, mean=mean, tw=None, head=2000, smooth=50, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000], no_measure=True)

        if mean:
            for m in result:
                k = m.keys()[0]
                if k not in results_:
                    results_[k] = []
                results_[k].append(m.values()[0])
        else:
            results.append(pd.concat([pd.DataFrame(r) for r in result], axis=1))

if not mean:
    df = results[0]
    for d in results[1:]:
        df = df.append(d)
else:
    df = pd.DataFrame(results_)

i = 'dfit'
# for j in df.columns:
#     print "%s with %s: %s" % (i, j, pearsonr(df[i], df[j]))
#     Plotter.plot_curve({'x': df[j], 'y': df[i]}, output_filename="%s_%s.png" % (i, j), linewidth=0, y_label=i, title=str(pearsonr(df[i], df[j])), x_label=j, figsize=(7, 4.5), loc=2, tight_layout=[], font_size=11, alpha=0.4, marker='.', colors="#e41a1c", grid=True, markersize=6)


print [r for r in df.corr()['dfit']]
pearsonr(df.x, df.dfit)
Plotter.plot_curve({'x': df.cd, 'y': df.dfit}, y_scale='log', linewidth=0, figsize=(5, 3.8), loc=2, tight_layout=[], font_size=11, alpha=0.1, marker='.', colors="#e41a1c", grid=True, markersize=3)
Plotter.plot_curve({'x': df.dfit, 'y': df['radius:#']}, linewidth=0, figsize=(5, 3.8), loc=2, tight_layout=[], font_size=11, alpha=0.4, marker='.', colors="#e41a1c", grid=True, markersize=6)

# F14: [-0.30529389079926983, 1.0, -0.72614528289662417, nan, nan, nan, nan, nan, nan, nan]
# F19 [0.61180348812168583, 1.0, -0.72732930733369694, nan, nan, nan, nan, nan, nan, nan]

from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

weights = np.ones_like(df.dfit)/len(df.dfit)
fig = plt.figure(figsize=(4.5, 2.6))
# plt.hist2d(np.log(df.,cd), -np.log(df.dfit), bins=100, cmap=plt.cm.gist_heat_r, weights=weights)
#plt.xscale('log')
plt.yscale('log')
plt.hist2d(df.cd, df.dfit, bins=100, cmap=plt.cm.gist_heat_r, weights=weights)
plt.xlabel("Communication diversity")
plt.ylabel("Fitness improvement")
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return "%0.3f" % x
    # return r'${}\times10^{{{}}}$'.format(a, b)
plt.colorbar(format=ticker.FuncFormatter(fmt), pad=0.03, ticks=[0.0000, 0.0010, 0.0021, 0.0032])
plt.tight_layout(rect=[0, 0, 1.06, 1])
plt.grid()
fig.text(0.2, 0.85, "$r=%0.2f$" % (pearsonr(np.log(df['dfit']), df['cd'])[0]))
plt.ylim(0, 0.005)
print pearsonr(df.cd, -np.log(df.dfit))
print pearsonr(np.log(df.cd), -np.log(df.dfit))
print pearsonr(df.cd, df.dfit)
plt.show()
plt.savefig("hist2d_dfit_cd.pdf")
measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#', 'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']


by_topologies = []
for r in range(0, 90, 30):
    res = results[r:r+30]
    dd = res[0]
    for d in res[1:]:
        dd = dd.append(d)
    by_topologies.append(dd)

for ddd in by_topologies:
    for m in measures:
        Plotter.plot_curve({'x': ddd.cd, 'y': ddd[m]}, linewidth=0, figsize=(5, 3.8), loc=2, tight_layout=[], font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6)
for m in measures:
    Plotter.plot_curve({'x': df.cd, 'y': df[m]}, linewidth=0, figsize=(5, 3.8), loc=2, tight_layout=[], font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6)
    print pearsonr(df.cd, df[m])


runs = 30
ds_index = range(0, len(topologies)*runs + 1, runs)
for measure in range(1, 7):
    ds = []
    # measure = 3
    away_std = 1.5
    for i,j in zip(ds_index, ds_index[1:]):
        d = pd.DataFrame({'x': df[measure][i:j], 'y': df[0][i:j]})
        dd = d[d.y < np.mean(d.y) + away_std*np.std(d.y)]
        dd = dd[dd.y > np.mean(dd.y) - away_std*np.std(dd.y)]
        dd = dd[dd.x < np.mean(dd.x) + away_std*np.std(dd.x)]
        dd = dd[dd.x > np.mean(dd.x) - away_std*np.std(dd.x)]
        ds.append(dd)
    ddd = ds[0]
    for d in ds[1:]:
        ddd = ddd.append(d)
    legends = topologies
    print pearsonr(ddd.x, ddd.y)
    Plotter.plot_curve(ds, linewidth=0, legends=legends, title=str(measure), figsize=(5, 3.8), loc=2, tight_layout=[], annotate=(0.49, 0.85, "$r=%0.2f$" % (ddd.corr()['x']['y'])), font_size=11, alpha=0.6, marker='.', colors=["#e41a1c", "#4daf4a", "#377eb8", "#fdae61"], grid=True, markersize=6)

    ddd = ds[1]
    Plotter.plot_curve(ddd, linewidth=0, legends=legends, title=str(measure), figsize=(5, 3.8), loc=2, tight_layout=[], annotate=(0.49, 0.85, "$r=%0.2f$" % (ddd.corr()['x']['y'])), font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6)

# the average and std dev -- CORRELATION
execfile("swarm_analyzer.py")

from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.stats import pearsonr

results = []
topologies = ['global', 'ring', 'vonneumann']

topologies = ['global']
for t in topologies:
    for i in range(30):
        filename = "./data/100_particles/%s_F06_%02d.teste.hdf" % (t, i)
        mean = False
        result = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_measures_and_cd(filename, mean=mean, tw=None, head=2000, diff_threshold=0.01, smooth=50, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
        results.append(pd.concat([pd.DataFrame(r) for r in result], axis=1))


df = results[0]
for d in results[1:]:
    df = df.append(d)

measures = ['normalized_average_around_center:#',  'aggregation_factor:#', 'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']

titles = ['Normalized\naverage\naround center',  'Aggregation\nfactor', 'Average\naround\ncenter', 'Average of\naverage\nall particles', 'Diameter', 'Radius']

labels = []
corrs = np.zeros((len(titles) - 1,len(titles) - 1))
for i in range(len(avg_array)):
    labels_i = []
    for j in range(len(avg_array[i])):
        corrs[i][j], p_value = pearsonr(df[measures[i]], df[measures[j]])
        labels_i.append("%0.2f\n(%0.2f)" % (corrs[i][j], p_value))
    labels.append(labels_i)
df_labels = pd.DataFrame(labels)
df_labels.columns = df.columns
df_labels.index = df.index
titles = ["%d" % i for i in range(1, len(titles) + 2)]
output_filename = "correlation_measures_global.pdf"
output_filename = None
Plotter.plot_heatmap(corrs, output_filename=output_filename, colorbar_on=False, values_on_text=df_labels, titles=titles, ordered=False, tight_layout=[], figsize=(6,6), cmap=plt.cm.coolwarm, values_on=True, font_size=14)





execfile("swarm_analyzer.py")

from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.stats import pearsonr
# the average and std dev (with amortization)
all = None
each = []
for i in range(30):
    print i
    filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
    filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
    #filename = "./data/100_particles/dynamicring_F06_%02d.teste.hdf" % i
    #filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
    df = pd.read_hdf(filename, 'df')
    measures = ['aggregation_factor:#', 'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']
    titles = [ 'Aggregation\nfactor', 'Normalized\naverage\naround center', 'Average\naround\ncenter', 'Average of\naverage\nall particles', 'Diameter', 'Radius']
    curves = {}
    k = 20
    for m in measures:
        ddy = [sum(df[m][n-k:n])/k for n in range(k, len(df))]
        ddy /= max(ddy)
        curves[m] =  ddy
    df = pd.DataFrame(curves)
    each.append(df[measures].corr())


avg = sum(each)/30
std_dev = [np.power(i - avg, 2) for i in each]
avg_array = np.array(avg)
std_dev = np.array(np.power(sum(std_dev)/30.0, 0.5))
labels = []
for i in range(len(avg_array)):
    labels_i = []
    for j in range(len(avg_array[i])):
        labels_i.append("%0.2f\n(%0.2f)" % (avg_array[i][j], std_dev[i][j]))
    labels.append(labels_i)


df = pd.DataFrame(labels)
df.columns = avg.columns
df.index = avg.index
titles = ["%d" % i for i in range(1, len(titles) + 2)]
Plotter.plot_heatmap(avg, output_filename="correlation_measures_ring_.png", values_on_text=df, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=14)
Plotter.plot_heatmap(avg, output_filename="correlation_measures_ring_.pdf", values_on_text=df, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=14)
#Plotter.plot_heatmap(avg, values_on_text=df, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=14)




    """

    @staticmethod
    def analyze_stagnation_correlation_communication_diversity(
            tw, topologies=None, function=6, functions=None, runs=30, **kargs):
        if topologies is None:
            topologies = ['ring', 'vonneumann', 'kregular5', 'kregular6', 'kregular7', 'kregular8', 'kregular9',
                          'kregular10', 'kregular20', 'kregular30', 'kregular40', 'kregular50',
                          'kregular60', 'kregular70', 'kregular80', 'kregular90', 'global']
        if functions is None:
            functions = [function]
        sums = []
        stagnation_its = []
        for function in functions:
            for t in topologies:
                for r in range(runs):
                    filename = "./data/100_particles/%s_F%02d_%02d.hdf" % (t, function, r)
                    print filename
                    sum_cd, stagnation_it = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative(filename,
                                                                                                         tw, **kargs)
                    sums.append(sum_cd)
                    stagnation_its.append(stagnation_it)
        df = pd.DataFrame({'sum_cd': sums, 's_it': stagnation_its})
        return df
    """
tw = 10
df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw, head=2000, diff_threshold=0.01, smooth=1)
df[10] = df['sum_cd']
tws = [25, 50, 75, 100, 200, 300, 400, 500, 1000]
for tw in tws:
    df[tw] = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw, head=2000, diff_threshold=0.01, smooth=1)['sum_cd']


execfile("swarm_analyzer.py")
execfile("plotter.py")
df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw=None, function=2, head=2000, diff_threshold=0.01, smooth=1, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
#df.sum_cd /= max(df.sum_cd)
#df.s_it /= max(df.s_it)
title = ""
xlabel = "Average communication diversity"
ylabel = "Stagnation iteration"
d = pd.DataFrame({'x': df.sum_cd, 'y': df.s_it})
#Plotter.plot_curve(d, linewidth=0, figsize=(5, 3.8), tight_layout=[], output_filename=output_filename, annotate=(0.15, 0.84, "$R=%0.2f$" % (d.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors='#de2d26', grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel, xlim=(0,1), ylim=(0,1))
# topologies = ['Ring', 'Neumann',  '10-regular', '20-regular', '30-regular', '40-regular', '50-regular', '60-regular', '70-regular', '80-regular', '90-regular', 'Global']
# topologies = ['Ring', 'Neumann',  '5-regular', '6-regular', '7-regular', '8-regular', '10-regular', '20-regular', '30-regular', '40-regular', '50-regular', '60-regular', '70-regular', '80-regular', '90-regular', 'Global']
topologies = ['Ring', 'Neumann'] + ['$k=%d$' % k for k in [5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90]] + ["Global"]
# topologies = ['$k=%d$' % k for k in [2, 4, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]]
runs = 30
ds_index = range(0, len(topologies)*runs + 1, runs)
ds = []
for i,j in zip(ds_index, ds_index[1:]):
    ds.append(pd.DataFrame({'x': df.sum_cd[i:j], 'y': df.s_it[i:j]}))
legends = topologies
#Plotter.plot_curve(ds, linewidth=0, legends=legends, figsize=(5, 3.8), loc=1, tight_layout=[], annotate=(0.39, 0.84, "$r=%0.2f$" % (d.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors=["#e41a1c", "#4daf4a", "#377eb8", "#fdae61"], grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))

away_std = 1.5
ds_no_outliers = []
for d in ds:
    dd = d[d.y < np.mean(d.y) + away_std*np.std(d.y)]
    dd = dd[dd.y > np.mean(dd.y) - away_std*np.std(dd.y)]
    dd = dd[dd.x < np.mean(dd.x) + away_std*np.std(dd.x)]
    dd = dd[dd.x > np.mean(dd.x) - away_std*np.std(dd.x)]
    ds_no_outliers.append(dd)

ddd = ds_no_outliers[0]
for d in ds[1:]:
    ddd = ddd.append(d)

# ddd1 /= ddd.max()
# ddd2 /= ddd.max()
# ddd3 /= ddd.max()
# ddd4 /= ddd.max()
len(ddd)
output_filename="correlation_cd_si_no_outliers.pdf"
output_filename= None
from scipy.stats import pearsonr
p = pearsonr(ddd['x'], ddd['y'])
colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#e0f3f8", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"]
colors = ["#e41a1c", "#4daf4a", "#377eb8", "#e0f3f8", "#fdae61", "#66c2a5", "#f46d43", "#6a3d9a"]
legends = ['Ring', 'Neumann'] + ['$k=%d$' % k for k in [5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90]] + ["Global$\quad\quad$ "]
# dss = ds_no_outliers[:2]+ds_no_outliers[3:]
dss = ds_no_outliers
dss.reverse()
legends.reverse()
colorblind = [(0, 0, 0), (230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66), (0, 114, 178), (213, 94, 0), (204, 121, 167)]
colorblind = ['#%02x%02x%02x' % (c[0], c[1], c[2]) for c in colorblind]
colors = colorblind[3], colorblind[5], colorblind[6], colorblind[1], colorblind[2], colorblind[7]
colors = colorblind[1:]
title = ""
xlabel = "Average communication diversity"
ylabel = "Stagnation iteration"
Plotter.plot_curve(dss, legend_ncol=2, ylim=(200, 1450), xlim=(0.05, 0.376), output_filename=output_filename, linewidth=0, legends=legends, mew=1.2, figsize=(4.8, 3.05), loc=2, tight_layout=[-0.03, -0.03, 1.03, 1.03], annotate=(0.66, 0.86, "$r=%0.2f$" % p[0]), font_size=11, alpha=1, marker=['+', 'x', '>', '^', 'o', 'v', 'D', 's', '3', '2', '<'], colors=colors, grid=True, markersize=4, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))


for i in range(len(topologies)):
    p = pearsonr(ds_no_outliers[i]['x'], ds_no_outliers[i]['y'])
    print "%s: %0.2f (%0.3f)" % (topologies[i], p[0], p[1])

for i in range(len(topologies)):
    p = pearsonr(ds[i]['x'], ds[i]['y'])
    print "%s: %0.2f (%0.4f)" % (topologies[i], p[0], p[1])

















##########3

execfile("swarm_analyzer.py")
execfile("plotter.py")
from scipy.stats import pearsonr

df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw=None, function=19, head=2000, diff_threshold=0.02, smooth=1, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
topologies = ['Ring', 'Neumann'] + ['$k=%d$' % k for k in [5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90]] + ["Global"]
runs = 30
ds_index = range(0, len(topologies)*runs + 1, runs)
ds = []
for i,j in zip(ds_index, ds_index[1:]):
    ds.append(pd.DataFrame({'x': df.sum_cd[i:j], 'y': df.s_it[i:j]}))

away_std = 1.5
ds_no_outliers = []
for d in ds:
    dd = d[d.y < np.mean(d.y) + away_std*np.std(d.y)]
    dd = dd[dd.y > np.mean(dd.y) - away_std*np.std(dd.y)]
    dd = dd[dd.x < np.mean(dd.x) + away_std*np.std(dd.x)]
    dd = dd[dd.x > np.mean(dd.x) - away_std*np.std(dd.x)]
    ds_no_outliers.append(dd)

ddd = ds_no_outliers[0]
for d in ds[1:]:
    ddd = ddd.append(d)

p = pearsonr(ddd['x'], ddd['y'])


#F14 2k (0.48173211798028143, 1.3798970467458267e-30) 4k (0.47799341554669339, 4.4700103566216424e-30)

print p

output_filename= None
legends = ['Ring', 'Neumann'] + ['$k=%d$' % k for k in [5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90]] + ["Global$\quad\quad$ "]
dss = ds_no_outliers
dss.reverse()
legends.reverse()
colorblind = [(0, 0, 0), (230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66), (0, 114, 178), (213, 94, 0), (204, 121, 167)]
colorblind = ['#%02x%02x%02x' % (c[0], c[1], c[2]) for c in colorblind]
colors = colorblind[3], colorblind[5], colorblind[6], colorblind[1], colorblind[2], colorblind[7]
colors = colorblind[1:]
Plotter.plot_curve(dss, legend_ncol=2, output_filename=output_filename, linewidth=0, legends=legends, mew=1.2, figsize=(4.8, 3.05), loc=2, tight_layout=[-0.03, -0.03, 1.03, 1.03], annotate=(0.66, 0.86, "$r=%0.2f$" % p[0]), font_size=11, alpha=1, marker=['+', 'x', '>', '^', 'o', 'v', 'D', 's', '3', '2', '<'], colors=colors, grid=True, markersize=4, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))



pearsonr([1, 2, 3], [4, 3, 7])

Plotter.plot_curve(ds_no_outliers, linewidth=0, y_scale='log', x_scale='log', legends=legends, figsize=(5, 3.8), loc=3, tight_layout=[], annotate=(0.39, 0.86, "$r=%0.2f$" % (ddd.corr()['y']['x'])), font_size=11, alpha=0.8, marker=['x', '.', 'o'], colors=["#e41a1c", "#4daf4a", "#377eb8", "#fdae61"], grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))
ddd_l = np.log(ddd)

    # boxplot stagnation and communication diversity
    execfile("plotter.py")
    execfile("swarm_analyzer.py")
    df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw=None, head=2000, diff_threshold=0.01, smooth=1, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
    topologies = [2, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    runs = 30
    ds_index = range(0, len(topologies)*runs + 1, runs)
    ds = []
    for i,j in zip(ds_index, ds_index[1:]):
    #for i,j in zip(ds_index, ds_index):
        ds.append(pd.DataFrame({'x': df.sum_cd[i:j], 'y': df.s_it[i:j]}))
    topologies_legend = map(str, topologies)
    output_filename="boxplot_stagnation.pdf"
    output_filename= None
    boxes_kargs = {'color': 'black', 'linewidth': 1.3, 'zorder': 3, 'fillstyle': 'full', 'facecolor': '#a6cee3'}
    means_kargs = {'color': 'black', 'fillstyle': 'full', 'markerfacecolor': "black", 'marker': "s", 'markersize': 3, 'mew': 1, 'mec': 'black', 'zorder': 5}
    fliers_kargs = {'color': 'black', 'marker': ".", 'markersize': 3, 'mew': 1.5, 'mec': 'black'}
    whiskers_kargs = {'color': 'black', 'linewidth': 1.2, 'zorder': 2, 'linestyle': '-', 'alpha': 1.0, 'mec': 'black'}
    medians_kargs = {'color': 'black', 'linewidth': 1.6, 'zorder': 5, 'alpha': 0.3}
    caps_kargs = {'linewidth': 1.5, 'color': 'black',}
    topologies_legend[0] = 'Ring'
    topologies_legend[len(topologies)-1] = 'Global'
    legends = [range(1, len(topologies_legend) + 1), topologies_legend, {"rotation": 90}]
    Plotter.plot_boxplots([d.y for d in ds], xticks_args=legends, grid_only='y', xlim=(0.3, len(topologies)+0.6), grid=True, output=output_filename, widths=0.7, tight_layout=[-0.1, -0.10, 1.025, 1.06], bootstrap=2000, boxes_kargs=boxes_kargs, showmeans=False, ylabel="Iteration of stagnation $i_s$", xlabel="$k$-regular topologies", size=(2.26, 2.3), showfliers=False, fliers_kargs=fliers_kargs, means_kargs=means_kargs, whiskers_kargs=whiskers_kargs, medians_kargs=medians_kargs, caps_kargs=caps_kargs)
    Plotter.plot_boxplots([d.x for d in ds], topologies_legend, ylabel="Communication diversity", size=(6, 3), showfliers=False)
    output_filename="communication_diversity_.pdf"
    output_filename= None
    #xticks_args = [2, 10, 30, 20, 40, 50, 60, 70, 80, 90, 100]
    xticks_args = [2, 10, 50, 100]
    xticks_args = (xticks_args, map(str, xticks_args))
    xticks_args[1][0] = 'Ring'
    xticks_args[1][len(xticks_args[1])-1] = 'Global'
    plt.clf()
    Plotter.plot_curve({'x': topologies, 'y': [np.mean(d.x) for d in ds]}, linewidth=1.5, x_scale='log', output_filename=output_filename, xticks_args=xticks_args, tight_layout=[-0.08, -0.10, 1.08, 1.05], mew=1.5, figsize=(4.5/2, 2.3), font_size=11, xtick_rotation=90, grid=True, marker=['s'], mec="#e41a1c", colors="#e41a1c", markersize=3, x_label="$k$-regular topologies", y_label="Average communication diversity")
    plt.plot(topologies, [np.mean(d.x) for d in ds])
    plt.xticks(topologies)
    plt.show()



execfile("swarm_analyzer.py")
import numpy as np, scipy.stats as st
import re

functions = [2, 6, 14, 19]
diff_threshold = 10**-5
df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(functions=functions, tw=None, diff_threshold=diff_threshold, smooth=1, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
len_topologies = 17
runs = 30
ds_index = range(0, len(df) + 1, runs)
ds_indices = zip(ds_index, ds_index[1:])
# ds_indices = [ds_indices[i] for j in range(len(functions)) for i in range(j, len(ds_indices), len_topologies) ]

fun_top = [(i, f) for f in functions for i in range(len_topologies)]
fun_top_i = 0
jump_every = 6
jump_every_i = 0
tex = ""
for i, j  in ds_indices:
    if fun_top[fun_top_i][1] == 19:
        sum_cd = list(df.sum_cd[i:j])
        l, u = st.t.interval(0.95, len(sum_cd)-1, loc=np.mean(sum_cd), scale=st.sem(sum_cd))
        # tex += "& $%0.3f, %.3f$  " % (l, u)
        tex += re.sub("0\.", ".", "& $(%-1.3f, %.3f)$  " % (l, u))
        jump_every_i += 1
        if jump_every_i % jump_every == 0:
            print tex + "\\\\"
            tex = ""
    fun_top_i += 1
print tex + "\\\\"


# execfile("plotter.py")
execfile("swarm_analyzer.py")

topologies = ['ring', 'kregular10', 'kregular30', 'kregular40']
functions = [2, 6, 14, 19]
diff_threshold = 10**-5
df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(topologies=topologies, functions=functions, tw=None, diff_threshold=diff_threshold, smooth=1, tws=[10, 25, 50, 75, 100, 200, 300, 400, 500, 1000])
# topologies = [2, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
topologies_n = [2, 10, 30, 40]
runs = 30
ds_index = range(0, len(functions)*len(topologies)*runs + 1, runs)
ds_indices = zip(ds_index, ds_index[1:])
ds_indices = [ds_indices[i] for j in range(len(functions)) for i in range(j, len(ds_indices), len(topologies)) ]
ds = []
space_i = 0
for i,j in ds_indices:
    ds.append(pd.DataFrame({'x': df.sum_cd[i:j], 'y': df.s_it[i:j]}))
    space_i += 1
    if space_i % len(topologies) == 0:
        ds.append(None)
topologies_legend = map(str, topologies)
output_filename="boxplot_stagnation.pdf"
output_filename= None
colorblind = [(0, 0, 0), (230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66), (0, 114, 178), (213, 94, 0), (204, 121, 167)]
colorblind = ['#%02x%02x%02x' % (c[0], c[1], c[2]) for c in colorblind]
colors = [colorblind[3]]*len(functions) + [colorblind[5]]*len(functions) + [colorblind[6]]*len(functions) + [colorblind[1]]*len(functions) #, colorblind[2], colorblind[7]
colors = [colorblind[3], colorblind[5], colorblind[6], colorblind[1], "black"]*4
boxes_kargs = {'color': 'black', 'linewidth': 1.3, 'zorder': 3, 'fillstyle': 'full', 'facecolor': colors}
means_kargs = {'color': 'black', 'fillstyle': 'full', 'markerfacecolor': "black", 'marker': "s", 'markersize': 3, 'mew': 1, 'mec': 'black', 'zorder': 5}
fliers_kargs = {'color': 'black', 'marker': ".", 'markersize': 3, 'mew': 1.5, 'mec': 'black'}
whiskers_kargs = {'color': 'black', 'linewidth': 1.2, 'zorder': 2, 'linestyle': '-', 'alpha': 1.0, 'mec': 'black'}
medians_kargs = {'color': 'black', 'linewidth': 1.6, 'zorder': 5, 'alpha': 0.3}
caps_kargs = {'linewidth': 1.5, 'color': 'black',}
topologies_legend[0] = 'Ring'
topologies_legend[len(topologies)-1] = 'Global'
topologies_legend = ['%d-regular' % i for i in topologies_n]
legends = [range(1, len(topologies_legend) + 1), topologies_legend, {"rotation": 90}]
legends = [[2.5, 7.5, 12.5, 17.5] + range(5, len(topologies)*len(functions) + len(topologies) + 1, 5), topologies_legend + [""]*4, {"rotation": 0}]
legends = None
Plotter.plos_style()
Plotter.plot_boxplots([d.x if d is not None else [] for d in ds ], whis=1.5, ylim=(0.1, 0.45), size=(4.7, 3), tight_layout=[-0.03, -0.04, 1.02, 1.04], loc=1, xticks_args=legends, legends=["$F_{%d}$" % i for i in [2, 6, 14, 19]], grid_only='y', xlim=(0.3, len(topologies) + len(topologies)*len(functions)), grid=True, output=output_filename, widths=0.7,  bootstrap=2000, boxes_kargs=boxes_kargs, showmeans=False, ylabel="Communication diversity", xlabel="$k$-regular topologies", showfliers=True, fliers_kargs=fliers_kargs, means_kargs=means_kargs, whiskers_kargs=whiskers_kargs, medians_kargs=medians_kargs, caps_kargs=caps_kargs)
Plotter.plot_boxplots([d.y if d is not None else [] for d in ds ], ylim=(1, 2200), xticks_args=legends, size=(4.7, 4), tight_layout=[-0.03, -0.04, 1.02, 1.04], legends=["$F_{%d}$" % i for i in [2, 6, 14, 19]], grid_only='y', xlim=(0.3, len(topologies) + len(topologies)*len(functions)), grid=True, output=output_filename, widths=0.7,  bootstrap=2000, boxes_kargs=boxes_kargs, showmeans=False, ylabel="Iteration of stagnation $i_s$", xlabel="$k$-regular topologies", showfliers=True, fliers_kargs=fliers_kargs, means_kargs=means_kargs, whiskers_kargs=whiskers_kargs, medians_kargs=medians_kargs, caps_kargs=caps_kargs)
    """

    @staticmethod
    def create_influence_graph_graphml(filename, output_file_name, window_size=1000, calculate_on=1000):
        influence_graph_grep = 'ig\:#'
        pre_callback = Callback.to_symmetric
        #for calculate_on in calculates_on:
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
    def plot_fitness(topologies=None, function=2, functions=None, runs=1, **kargs):
        if topologies is None:
            topologies = ['ring', 'vonneumann', 'kregular5', 'kregular6', 'kregular7', 'kregular8', 'kregular9',
                          'kregular10', 'kregular20', 'kregular30', 'kregular40', 'kregular50',
                          'kregular60', 'kregular70', 'kregular80', 'kregular90', 'global']
        if functions is None:
            functions = [function]
        for function in functions:
            for t in topologies:
                for r in range(runs):
                    filename = "./data/100_particles/%s_F%02d_%02d.hdf" % (t, function, r)
                    diff_threshold, iter_threshold = 0.02, 500
                    _, stagnation_iteration = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative(filename, 10, smooth=1, diff_threshold=diff_threshold, iter_threshold=iter_threshold)
                    df = pd.read_hdf(filename, 'df')
                    Plotter.plot_curve({'x': df.index, 'y': df['it:#']}, vline_at=stagnation_iteration, figsize=(10, 4))

    @staticmethod
    def analyze_analysis_df_file(filename):
        df = pd.read_hdf(filename, 'df')
    """
    execfile("swarm_analyzer.py")

    filename = "./data/100_particles/ring_F06_00.teste.hdf"
    df = pd.read_hdf(filename, 'df')
    measures = ['aggregation_factor:#', 'average_around_center:#', 'average_of_average_around_all_particles:#', 'coherence:#', 'diameter:#', 'normalized_average_around_center:#', 'radius:#']
    titles = ['Aggregation\nfactor', 'Average\naround\ncenter', 'Average of\naverage\nall particles', 'Coherence', 'Diameter', 'Normalized\naverage\naround center', 'Radius']

    execfile("plotter.py")
    # for each
    i = 0
    for i in range(1):
        filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
        df = pd.read_hdf(filename, 'df').head(2000)
        df.to_hdf
        k = 50
        tw = 100
        x = df.index[k:len(df)]
        y = [sum(df[tw][n-k:n])/k for n in range(k, len(df))]
        curves = [{'x': x, 'y': y}]
        a = np.array(df['it:#'])


        b = a[:len(a)-1] - a[1:len(a)]
        b /= max(b)
        a = np.array(df['it:#']/max(df['it:#']))
        curves += [{'x': df.index, 'y': a}]
        measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']

        d = (a[:len(a)-1] - a[1:len(a)]) > 0.0001
        counts = []
        last = d[0]
        last_c = 1
        for c in d[1:]:
            if c == last:
                last_c += 1
            else:
                counts.append(last_c)
                last_c = 1
            last = c
        if last_c != 1:
            counts.append(last_c)

        threshold = 500

        for threshold_i in range(len(counts)):
            if counts[threshold_i] > threshold:
                break
        stagnation_iteration = sum(counts[:threshold_i])

        #curves += [{'x': df.index[1:], 'y': b}]
        for m in measures:
            ddx = df.index[k:len(df)]
            ddy = [sum(df[m][n-k:n])/k for n in range(k, len(df))]
            ddy /= max(ddy)
            curves += [{'x': ddx, 'y': ddy}]
        #Plotter.plot_curve(curves, legends=['auc', 'it']+measures, vline_at=stagnation_iteration, figsize=(15,4), markersize=3, output_filename="auc_smoothed_vonneumann_%02d.png" % i, ylim=(0, 1.0), grid=True)
        Plotter.plot_curve(curves, legends=['auc', 'it']+measures, vline_at=stagnation_iteration, figsize=(15,4), markersize=3, ylim=(0, 1.0), grid=True)

    execfile("swarm_analyzer.py")
    execfile("plotter.py")
    # for each
    for i in range(1):
        filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
        df = pd.read_hdf(filename, 'df')
        measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']
        titles = [ 'Aggregation\nfactor', 'Coherence', 'Normalized\naverage\naround center', 'Average\naround\ncenter', 'Average of\naverage\nall particles', 'Diameter', 'Radius']

        Plotter.plot_heatmap(df[measures].corr(), output_filename="correlation_measures_vonneumann_%02d.png" % i, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)

    Plotter.plot_heatmap(df[measures].corr(), titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.coolwarm, values_on=True, font_size=11)
    Plotter.plot_heatmap(df[measures].corr(), output_filename=None, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.coolwarm, values_on=True, font_size=14)

    # the average and std dev
    all = None
    each = []
    for i in range(30):
        filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
        filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
        # filename = "./data/100_particles/dynamicring_F06_%02d.teste.hdf" % i
        # filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
        df = pd.read_hdf(filename, 'df')
        measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']
        titles = [ 'Aggregation\nfactor', 'Coherence', 'Normalized\naverage\naround center', 'Average\naround\ncenter', 'Average of\naverage\nall particles', 'Diameter', 'Radius']
        each.append(df[measures].corr())
    avg = sum(each)/30
    std_dev = [np.power(i - avg, 2) for i in each]
    avg_array = np.array(avg)
    std_dev = np.array(np.power(sum(std_dev)/30.0, 0.5))
    labels = []
    for i in range(len(avg_array)):
        labels_i = []
        for j in range(len(avg_array[i])):
            labels_i.append("%0.2f\n(%0.2f)" % (avg_array[i][j], std_dev[i][j]))
        labels.append(labels_i)
    df = pd.DataFrame(labels)
    df.columns = avg.columns
    df.index = avg.index
    titles = ["%d" % i for i in range(1, len(titles) + 2)]
    Plotter.plot_heatmap(avg, output_filename="correlation_measures_ring.pdf", colorbar_on=False, values_on_text=df, titles=titles, ordered=False, tight_layout=[], figsize=(6,6), cmap=plt.cm.coolwarm, values_on=True, font_size=14)
    Plotter.plot_heatmap(avg, output_filename=None, values_on_text=df, titles=titles, ordered=False, colorbar_on=False, tight_layout=[], figsize=(6,6), cmap=plt.cm.coolwarm, values_on=True, font_size=14)

    # the average and std dev (with amortization)
    all = None
    each = []
    for i in range(30):
        print i
        filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
        filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
        filename = "./data/100_particles/dynamicring_F06_%02d.teste.hdf" % i
        filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
        df = pd.read_hdf(filename, 'df')
        measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']
        titles = [ 'Aggregation\nfactor', 'Coherence', 'Normalized\naverage\naround center', 'Average\naround\ncenter', 'Average of\naverage\nall particles', 'Diameter', 'Radius']
        curves = {}
        k = 20
        for m in measures:
            ddy = [sum(df[m][n-k:n])/k for n in range(k, len(df))]
            ddy /= max(ddy)
            curves[m] =  ddy
        df = pd.DataFrame(curves)
        each.append(df[measures].corr())
    avg = sum(each)/30
    std_dev = [np.power(i - avg, 2) for i in each]
    avg_array = np.array(avg)
    std_dev = np.array(np.power(sum(std_dev)/30.0, 0.5))
    labels = []
    for i in range(len(avg_array)):
        labels_i = []
        for j in range(len(avg_array[i])):
            labels_i.append("%0.2f\n(%0.2f)" % (avg_array[i][j], std_dev[i][j]))
        labels.append(labels_i)
    df = pd.DataFrame(labels)
    df.columns = avg.columns
    df.index = avg.index
    titles = ["%d" % i for i in range(1, len(titles) + 2)]
    #Plotter.plot_heatmap(avg, output_filename="correlation_measures_vonneumann_.pdf", values_on_text=df, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=14)
    Plotter.plot_heatmap(avg, values_on_text=df, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=14)

    # only some measures and AUC of a run
    filename = "./data/100_particles/dynamicring_F06_00.teste.hdf"
    for i in range(30):
        filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
        df = pd.read_hdf(filename, 'df')
        df.columns = map(str, df.columns)
        #df = df.tail(len(df)).head(1000)
        measures = ['coherence:#', 'aggregation_factor:#',  'average_of_average_around_all_particles:#', 'normalized_average_around_center:#', 'radius:#']
        measures += [u'10', u'25', u'50', u'75', u'100', u'200', u'300', u'400', u'500', u'1000']
        titles_y = ['Coherence','Aggregation\nfactor',  'Average of\naverage\nall particles', 'Normalized\naverage\naround center', 'Radius']
        titles_x = [u'10', u'25', u'50', u'75', u'100', u'200', u'300', u'400', u'500', u'1000']
        only_measures = df[measures].corr()
        only_measures.columns = map(str, only_measures.columns)
        for tw in tws:
            del only_measures[tw]
        only_measures = only_measures.tail(len(tws))
        Plotter.plot_heatmap(only_measures.transpose(), main_title="Global topology", output_filename="correlations_global2_%02d.png" % i, ordered=False, titles_x=titles_x, titles_y=titles_y, tight_layout=[], figsize=(9,3.5), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)

    # average/std dev some measures and AUC of runs
    for maximum in [2000]: #[1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]:
        each = []
        for i in range(30):
            filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
            df = pd.read_hdf(filename, 'df')
            df.columns = map(str, df.columns)
            df = df.tail(len(df)).head(maximum) #### remeber this!
            measures = ['coherence:#', 'aggregation_factor:#',  'average_of_average_around_all_particles:#', 'normalized_average_around_center:#', 'radius:#']
            measures += [u'10', u'25', u'50', u'75', u'100', u'200', u'300', u'400', u'500', u'1000']
            titles_y = ['Coherence','Aggregation\nfactor',  'Average of\naverage\nall particles', 'Normalized\naverage\naround center', 'Radius']
                    #1: aggregation factor, 2: coherence, 3: normalized average around center, 4: average around center, 5: average of average all particles, 6: diameter, 7: radius
            titles_y = [1, 2, 3,
            titles_x = [u'10', u'25', u'50', u'75', u'100', u'200', u'300', u'400', u'500', u'1000']
            only_measures = df[measures].corr()
            only_measures.columns = map(str, only_measures.columns)
            for tw in titles_x:
                del only_measures[tw]
            only_measures = only_measures.tail(len(titles_x))
            each.append(only_measures.transpose())
        avg = sum(each)/30
        std_dev = [np.power(i - avg, 2) for i in each]
        avg_array = np.array(avg)
        std_dev = np.array(np.power(sum(std_dev)/30.0, 0.5))
        labels = []
        for i in range(len(avg_array)):
            labels_i = []
            for j in range(len(avg_array[i])):
                labels_i.append("%0.2f\n(%0.2f)" % (avg_array[i][j], std_dev[i][j]))
            labels.append(labels_i)
        df = pd.DataFrame(labels)
        df.columns = avg.columns
        df.index = avg.index
        titles = ["%d" % i for i in range(1, len(titles) + 2)]
        Plotter.plot_heatmap(avg, values_on_text=df, main_title="Global topology", ordered=False, titles_x=titles_x, titles_y=titles_y, tight_layout=[], figsize=(9,3.5), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)
        #Plotter.plot_heatmap(avg, values_on_text=df, main_title="Global topology", ordered=False, output_filename="correlations_global_it%d.png" % maximum, titles_x=titles_x, titles_y=titles_y, tight_layout=[], figsize=(9,3.5), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)

    # average/std dev some measures and AUC of runs with amortization
    for maximum in [2000]: #[1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]:
        each = []
        for i in range(30):
            filename = "./data/100_particles/dynamicring_F06_%02d.teste.hdf" % i
            df = pd.read_hdf(filename, 'df')
            df.columns = map(str, df.columns)
            df = df.tail(len(df)).head(maximum) #### remeber this!
            measures = ['coherence:#', 'aggregation_factor:#',  'average_of_average_around_all_particles:#', 'normalized_average_around_center:#', 'radius:#']
            measures += [u'10', u'25', u'50', u'75', u'100', u'200', u'300', u'400', u'500', u'1000']
            titles_y = ['Coherence','Aggregation\nfactor',  'Average of\naverage\nall particles', 'Normalized\naverage\naround center', 'Radius']
            titles_x = [u'10', u'25', u'50', u'75', u'100', u'200', u'300', u'400', u'500', u'1000']
            curves = {}
            k = 100
            for m in measures:
                ddy = [sum(df[m][n-k:n])/k for n in range(k, len(df))]
                ddy /= max(ddy)
                curves[m] =  ddy
            df = pd.DataFrame(curves)
            only_measures = df[measures].corr()
            only_measures.columns = map(str, only_measures.columns)
            for tw in titles_x:
                del only_measures[tw]
            only_measures = only_measures.tail(len(titles_x))
            each.append(only_measures.transpose())
        avg = sum(each)/30
        std_dev = [np.power(i - avg, 2) for i in each]
        avg_array = np.array(avg)
        std_dev = np.array(np.power(sum(std_dev)/30.0, 0.5))
        labels = []
        for i in range(len(avg_array)):
            labels_i = []
            for j in range(len(avg_array[i])):
                labels_i.append("%0.2f\n(%0.2f)" % (avg_array[i][j], std_dev[i][j]))
            labels.append(labels_i)
        df = pd.DataFrame(labels)
        df.columns = avg.columns
        df.index = avg.index
        titles = ["%d" % i for i in range(1, len(titles) + 2)]
        Plotter.plot_heatmap(avg, values_on_text=df, main_title="Global topology", ordered=False, titles_x=titles_x, titles_y=titles_y, tight_layout=[], figsize=(9,3.5), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)
        #Plotter.plot_heatmap(avg, values_on_text=df, main_title="Global topology", ordered=False, output_filename="correlations_global_it%d.png" % maximum, titles_x=titles_x, titles_y=titles_y, tight_layout=[], figsize=(9,3.5), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)



    Plotter.plot_heatmap(avg, values_on_text=df, main_title="Global topology", output_filename="correlations_global.pdf", ordered=False, titles_x=titles_x, titles_y=titles_y, tight_layout=[], figsize=(9,3.5), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)


    filename = "./data/100_particles/dynamicring_F06_00.teste.hdf"
    filename = "./data/100_particles/ring_F06_00.teste.hdf"
    filename = "./data/100_particles/global_F06_00.teste.hdf"
    df = pd.read_hdf(filename, 'df')
    df.columns = map(str, df.columns)
    tws = [u'10', u'25', u'50', u'75', u'100', u'200', u'300', u'400', u'500', u'1000']
    titles = map(str, tws)
    Plotter.plot_heatmap(df[tws].corr(), titles=tws, ordered=True, tight_layout=[], figsize=(7,6), cmap=plt.cm.Spectral, values_on=True, font_size=11)

    execfile("plotter.py")
    topologies = ["vonneumann", "ring", "dynamicring", "global", "30kregular"]
    topology = "vonneumann"
    for topology in topologies:
        df = pd.read_hdf("./data/100_particles/number_of_components/heatmap_1000_%s.hdf" % topology, 'df')
        components = np.array(df)#/100
        ws = list(np.arange(0, 1.01, 0.01))
        yticks = range(10, len(ws), 10)
        titles_y = map(str, [ws[i] for i in yticks])
        tws = range(1, 1001, 10) + [1000]
        xticks = range(0, len(tws), 50)
        titles_x = map(str, [1]+[tws[i]-1 for i in xticks[1:len(xticks)-1]]+[1000])
        # yticks.reverse()
        # titles_x = None
        # titles_y = None
        # xticks = [0]
        # yticks = [0]
        output_filename = "heatmap_%s.svg" % topology
        # output_filename = None
        plt.clf()
        cmap = plt.cm.Spectral_r
        # import seaborn.apionly as sns
        # cmap = sns.cubehelix_palette(20, start=2.8, light=0.9, dark=0.1, rot=0.75, as_cmap=True)
        Plotter.plot_heatmap(components, colorbar_on=False, output_filename=output_filename, subplot_adjust=[0.17, 0.19], font_size=16, vmin=0, vmax=100, tight_layout=[-0.03, 0, 1.06, 1.06], figsize=(2.5, 2.3), ordered=False, titles_x=titles_x, set_xticks=xticks, x_label="Time window", y_label="Filter", titles_y=titles_y, set_yticks=yticks, grid=False, cmap=cmap)

    ### tests
    #for i in range(len(ws)):
    #for j in range(len(tws)):
    #components[i][j] = SwarmAnalyzer.get_number_of_components(filename, tw, w, calculate_on=1000)[0][1]
#    Plotter.plot_heatmap(np.flipud(components), vmin=components.min(), vmax=components.max(), figsize=(5,3.5), ordered=False, titles_x=titles_x, set_xticks=xticks, titles_y=titles_y, set_yticks=yticks, grid=False, y_label="Time window", x_label="Filter", cmap=plt.cm.Spectral_r)
     ws = range(10)
     tws = range(10)
     intt = 0
     components = np.zeros((len(ws), len(tws)))
     for i in range(len(ws)):
         for j in range(len(tws)):
             components[i][j] = intt
             intt += 1
    titles_y = ws
    titles_y.reverse()
    Plotter.plot_heatmap(components, vmin=components.min(), vmax=components.max())
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
        array = array[:len(array)-n] - array[n:]
        data = {'y': list(array), 'x': list(df['x'][:len(df['x']) - n])}
        return data

    @staticmethod
    def get_giant_component_destruction_area_from_files(basename, window_size, runs=30):
        filenames = [basename+"%02d" % i for i in range(1, runs+1)]
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
        df.to_hdf(basename+str(window_size)+".hdf", 'df')

    @staticmethod
    def read_hdfs_and_plot(basename):
        filenames = [basename+"%02d" % i for i in range(1, 30)]
        windows_size = 1000
        for filename in filenames:
            df = pd.read_hdf(filename+"_"+str(windows_size)+".hdf", 'df')
            Plotter.plot_curve(df, figsize=(18, 6))

    @staticmethod
    def read_files_and_plot(filenames, windows_size, calculate_on):
        influence_graph_grep = 'ig\:#'
        pre_callback = Callback.to_symmetric
        graph_matrices, _ = SwarmParser.read_files_and_measures(filenames, influence_graph_grep=influence_graph_grep,
                                                                pos_callback=pre_callback, windows_size=windows_size,
                                                                calculate_on=calculate_on)
        normalize = [2*i for i in windows_size]
        pd_datas = []
        for title, _ in filenames:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(graphs,
                                                                             weight_normalize=normalize)
            pd_datas.append((title, dict(zip(windows_size, curves_areas))))
        GiantComponentDeathPlotter.giant_component_death_curve(calculate_on, pd_datas, windows_size, xlim=(0, 1.0), figsize=(4.5, 4))
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
        normalize = [2*i for i in windows_size]
        pd_datas = []
        for title,_ in filenames:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(graphs,
                                                                             weight_normalize=normalize)
            pd_datas.append((title, dict(zip(windows_size, curves_areas))))
        GiantComponentDeathPlotter.giant_component_death_curve_with_area(pd_datas, xlim=(0, 1.0), figsize=(4.5, 1.9), mew=1.2, tight_layout=[-0.02, 0.01, 1.02, 1.06], output_filename="graph_destruction_area.pdf") #, **kargs)
    """
    plt.clf()
    execfile("swarm_analyzer.py")
    execfile("giant_component_analysis_plotter.py")
    filenames = [('Global', "./data/100_particles/global_F06_06.teste"), ('Ring', "./data/100_particles/ring_F06_06.teste")]
    SwarmAnalyzer.read_files_and_plot_with_area(filenames, windows_size=[100, 1000], calculate_on=1000, output_filename="graph_destruction_area_10002.pdf")
    SwarmAnalyzer.read_files_and_plot_with_area(filenames, windows_size=[100, 1000], calculate_on=1000)
    """
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

#if __name__ == "__main__":
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


