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
topology = "global"
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

    #Plotter.plot_curve(curves, legends=legends, marker=['.', '+', 'x'], markersize=4, font_size=10, markevery=17, alpha=0.9, tight_layout=[], figsize=(5,4), grid=True, linewidth=0.7, xlim=(window_size, 1400), y_label="Number of\ncomponents", x_label="Iterations")
    Plotter.plot_subplot_curve([curvess2, curvess1], titles=['Von Neumann', 'Ring'], legends=legends, marker=['.', '+', 'x'], markersize=4, font_size=13, markevery=17, alpha=0.9, tight_layout=[], figsize=(13, 4), grid=True, linewidth=0.7, xlim=(window_size, 3000), y_label="Number of\ncomponents", x_label="Iterations")
    Plotter.plot_subplot_curve([curvess2, curvess1], colors=["#e41a1c", "#4daf4a", "#377eb8", "#fdae61"], titles=['Ring', 'Von Neumann'], legends=legends, marker=['v', '+', '.', 'x'], markevery=2, markersize=7, font_size=13, tight_layout=[], figsize=(13, 4), grid=True, linewidth=1.2, xlim=(window_size, 3000), y_label="Number of\ncomponents", x_label="Iterations", output_filename="number_of_components_2.pdf")
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
execfile("swarm_analyzer.py")
base_name = './data/100_particles/'
topology = 'vonneumann'
function = 6
suffix = '.teste'
filenames = ['%s/%s_F%02d_%02d%s' % (base_name, topology, function, run, suffix) for run in range(14, 30)]
for filename in filenames:
    print time.strftime(">> %H:%M:%S-%d/%m/%y")
    df = SwarmAnalyzer.create_swarm_analysis_df_file_from_file(filename, filename+'.hdf')

    """

    @staticmethod
    def analyze_analysis_df_file_stagnation_derivative(filename, tw, smooth=50, diff_threshold=0.0001,
                                                       iter_threshold=500, head=None, plot=False):
        df = pd.read_hdf(filename, 'df')
        if head:
            df = df.head(head)
        x = df.index[smooth:len(df)]
        communication_diversity_smoothed = [sum(df[tw][n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
        # curves = [{'x': x, 'y': communication_diversity_smoothed}]
        fitness = np.array(df['it:#'][smooth:])
        # curves += [{'x': x, 'y': fitness/fitness.max()}]
        # measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']

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
        communication_diversity_until_stagnation = communication_diversity_smoothed[:stagnation_iteration]
        if plot:
            Plotter.plot_curve({'x': x, 'y': fitness/max(fitness)}, vline_at=stagnation_iteration, figsize=(15, 4),
                               markersize=3, ylim=(0, 1.0), grid=True)

        return sum(communication_diversity_until_stagnation)/float(stagnation_iteration), float(stagnation_iteration)
    """
    execfile("swarm_analyzer.py")
    i = 10
    filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
    filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
    sum_cd, stagnation_it = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative(filename, 100, head=2000, diff_threshold=0.02, plot=True)

    """

    @staticmethod
    def analyze_analysis_df_file_stagnation_derivative_and_measures(filename, tw, smooth=50, diff_threshold=0.0001,
                                                                    iter_threshold=500, head=None, **kargs):
        df = pd.read_hdf(filename, 'df')
        if head:
            df = df.head(head)
        x = df.index[smooth:len(df)]
        communication_diversity_smoothed = [sum(df[tw][n-smooth:n])/float(smooth) for n in range(smooth, len(df))]
        curves = [{'x': x, 'y': communication_diversity_smoothed}]
        fitness = np.array(df['it:#'][smooth:])
        curves += [{'x': x, 'y': fitness/fitness.max()}]
        measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']

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
        Plotter.plot_curve(curves, legends=legends, vline_at=stagnation_iteration, figsize=(15, 4),
                           markersize=3, ylim=(0, 1.0), grid=True, loc=1, **kargs)
    """
    execfile("swarm_analyzer.py")
    i = 10
    tw = 500
    for i in range(1):
        filename = "./data/100_particles/dynamicring_F06_%02d.teste.hdf" % i
        # filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
        # filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
        # filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
        output_filename="auc_smoothed_ring_%03d_%02d.png" % (tw, i)
        output_filename="auc_smoothed_vonneumann_%03d_%02d.png" % (tw, i)
        output_filename="auc_smoothed_global_%03d_%02d.png" % (tw, i)
        output_filename="auc_smoothed_dynamicring_%03d_%02d.png" % (tw, i)
        SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative_and_measures(filename, tw, head=3000, diff_threshold=0.02, output_filename=output_filename, smooth=1)


    """

    @staticmethod
    def analyze_stagnation_correlation_communication_diversity(tw, **kargs):
        runs = 30
        topologies = ['vonneumann', 'ring', 'global', 'dynamicring']
        #topologies = ['vonneumann', 'ring', 'dynamicring']
        #topologies = ['global']
        sums = []
        stagnation_its = []
        for t in topologies:
            for r in range(runs):
                filename = "./data/100_particles/%s_F06_%02d.teste.hdf" % (t, r)
                sum_cd, stagnation_it = SwarmAnalyzer.analyze_analysis_df_file_stagnation_derivative(filename,
                                                                                                     tw, **kargs)
                sums.append(sum_cd)
                stagnation_its.append(stagnation_it)
        df = pd.DataFrame({'sum_cd': sums, 's_it': stagnation_its})
        return df
    """
execfile("swarm_analyzer.py")
execfile("plotter.py")
for tw in [10, 25, 50, 75, 100, 200, 300, 400, 500, 1000]:
    df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw, head=3000, diff_threshold=0.02, smooth=1)
    #df.sum_cd /= max(df.sum_cd)
    #df.s_it /= max(df.s_it)
    title = ""
    xlabel = "Average communication diversity"
    ylabel = "Stagnation iteration"
    d = pd.DataFrame({'x': df.sum_cd, 'y': df.s_it})
    #Plotter.plot_curve(d, linewidth=0, figsize=(5, 3.8), tight_layout=[], output_filename=output_filename, annotate=(0.15, 0.84, "$R=%0.2f$" % (d.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors='#de2d26', grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel, xlim=(0,1), ylim=(0,1))
    d1 = pd.DataFrame({'x': df.sum_cd[:30], 'y': df.s_it[:30]})
    d2 = pd.DataFrame({'x': df.sum_cd[30:60], 'y': df.s_it[30:60]})
    d3 = pd.DataFrame({'x': df.sum_cd[60:90], 'y': df.s_it[60:90]})
    d4 = pd.DataFrame({'x': df.sum_cd[90:120], 'y': df.s_it[90:120]})
    legends = ['Neumann', 'Ring', 'Global', 'Dynamic']
    #legends = ['Neumann', 'Ring', 'Dynamic']
    output_filename = "correlation_cd_si_normalized.pdf"
    # d1.x /= max(df.sum_cd)
    # d2.x /= max(df.sum_cd)
    # d3.x /= max(df.sum_cd)
    # d4.x /= max(df.sum_cd)
    # d1.y /= max(df.s_it)
    # d2.y /= max(df.s_it)
    # d3.y /= max(df.s_it)
    # d4.y /= max(df.s_it)
    # Plotter.plot_curve([d1, d2, d3, d4], linewidth=0, legends=legends, figsize=(5, 3.8), loc=1, tight_layout=[], output_filename=output_filename, annotate=(0.39, 0.84, "$r=%0.2f$" % (d.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors=["#e41a1c", "#4daf4a", "#377eb8", "#fdae61"], grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))
    dd1 = d1
    dd2 = d2
    dd3 = d3
    dd4 = d4
    # dd1.x /= max(dd1.x)
    # dd1.y /= max(df.s_it)
    # dd2.x /= max(dd2.x)
    # dd2.y /= max(df.s_it)
    # dd3.x /= max(dd3.x)
    # dd3.y /= max(df.s_it)
    # dd4.x /= max(dd4.x)
    # dd4.y /= max(df.s_it)
    # dd = dd4[dd4.y < 1000]
    # dd = dd2[dd2.y < 1000].append(dd4[dd4.y < 1000]).append(dd3).append(dd1[dd1.y < 1000])
    # dd = dd3
    # output_filename = "correlation_cd_si_no_outliers.pdf"
    # Plotter.plot_curve(dd, linewidth=0, legends=['Neumann', 'Ring', 'Global', 'Dynamic'], figsize=(5, 3.8), loc=4, tight_layout=[], output_filename=output_filename, annotate=(0.15, 0.84, "$r=%0.2f$" % (dd.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors="#4daf4a", grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))

    ddd2 = dd2[dd2.y < 1000]
    ddd4 = dd4[dd4.y < 1000]
    ddd3 = dd3
    ddd1 = dd1[dd1.y < 1000]
    # ddd = ddd1.append(ddd2.append(ddd3.append(ddd4)))

    # ddd1 /= ddd.max()
    # ddd2 /= ddd.max()
    # ddd3 /= ddd.max()
    # ddd4 /= ddd.max()
    # Plotter.plot_curve([ddd1, ddd2, ddd3, ddd4], linewidth=0, legends=legends, figsize=(5, 3.8), loc=1, tight_layout=[], output_filename=output_filename, annotate=(0.39, 0.84, "$r=%0.2f$" % (ddd.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors=["#e41a1c", "#4daf4a", "#377eb8", "#fdae61"], grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))

    ddds = [ddd1, ddd2, ddd3, ddd4]
    for i in range(len(ddds)):
        title = legends[i] + " $t_w=%d$" %tw
        output_filename = "correlation_cd_si_no_outliers_%d_%s.pdf" %  (tw, title)
        ddd = ddds[i]
        Plotter.plot_curve(ddd, linewidth=0, legends=legends, figsize=(5, 3.8), loc=1, tight_layout=[], output_filename=output_filename, annotate=(0.39, 0.84, "$r=%0.2f$" % (ddd.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))
        output_filename = "correlation_cd_si_no_outliers_%s.png" %  title
        Plotter.plot_curve(ddd, linewidth=0, legends=legends, figsize=(5, 3.8), loc=1, tight_layout=[], output_filename=output_filename, annotate=(0.39, 0.84, "$r=%0.2f$" % (ddd.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors="#e41a1c", grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel) #, xlim=(0,1), ylim=(0,1))

execfile("swarm_analyzer.py")
execfile("plotter.py")
tws = [10, 25, 50, 75, 100, 200, 300, 400, 500, 1000]
correlations = {}
for tw in tws:
    print "heatmap_correlation_cd_si_%04d" % tw
    df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw, head=3000, diff_threshold=0.02, smooth=1)

    # outliers: (would need to use std dev and etc, but 1000 is pretty accurate)
    d1 = df[:30]
    d2 = df[30:60]
    d3 = df[60:90]
    d4 = df[90:120]

    d1 = d1[d1.s_it < 1000]
    d2 = d2[d2.s_it < 1000]
    d3 = d3[d3.s_it < 1000]
    d4 = d4[d4.s_it < 1000]
    d1 = pd.DataFrame({'x': d1.sum_cd, 'y': d1.s_it})
    d2 = pd.DataFrame({'x': d2.sum_cd, 'y': d2.s_it})
    d3 = pd.DataFrame({'x': d3.sum_cd, 'y': d3.s_it})
    d4 = pd.DataFrame({'x': d4.sum_cd, 'y': d4.s_it})
    correlations[tw] = {}
    correlations[tw]['vonneumann'] = d1.corr()['y']['x']
    correlations[tw]['ring'] = d2.corr()['y']['x']
    correlations[tw]['global'] = d3.corr()['y']['x']
    correlations[tw]['dynamicring'] = d4.corr()['y']['x']

output_filename = "heatmap_correlation_cd_si_no_outliers.png"
Plotter.plot_heatmap(pd.DataFrame(correlations), titles_y=['Dynamic', 'Global', 'Ring', 'Neumann'], output_filename=output_filename, x_label="Time Window", titles_x=map(str, tws), ordered=False, tight_layout=[], figsize=(10, 3.5), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)
    # print "%03d:  %f" % (iter, df.corr()['y']['x'])

    # boxplot stagnation and communication diversity
    execfile("plotter.py")
    tws = [10, 25, 50, 75, 100, 200, 300, 400, 500, 1000]
    correlations = {}
    for tw in tws:
        df = SwarmAnalyzer.analyze_stagnation_correlation_communication_diversity(tw, head=3000, diff_threshold=0.02, smooth=1)
        # ['vonneumann', 'ring', 'global', 'dynamicring']
        d4 = pd.DataFrame({'x': df.sum_cd[:30], 'y': df.s_it[:30]})
        d3 = pd.DataFrame({'x': df.sum_cd[30:60], 'y': df.s_it[30:60]})
        d2 = pd.DataFrame({'x': df.sum_cd[60:90], 'y': df.s_it[60:90]})
        d1 = pd.DataFrame({'x': df.sum_cd[90:120], 'y': df.s_it[90:120]})
        Plotter.plot_boxplots([d1.x, d2.x, d3.x, d4.x], ['Dynamic', 'Global', 'Ring', 'Neumann'], "Communication diversity ($tw=%d$)" % tw, ylabel="Average communication\ndiversity until stagnation", size=(3.4,5), ylim=(0,1), output="boxplot_communication_diversity_%03d.pdf" % tw)
    Plotter.plot_boxplots([d1.y, d2.y, d3.y, d4.y], ['Dynamic', 'Global', 'Ring', 'Neumann'], "Swarm stagnation", ylabel="Iteration of stagnation $i_s$", size=(3.4, 5), output="boxplot_stagnation.pdf")
    """


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
    for i in range(30):
        filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
        df = pd.read_hdf(filename, 'df')
        measures = ['aggregation_factor:#', 'coherence:#',  'normalized_average_around_center:#',  'average_around_center:#', 'average_of_average_around_all_particles:#', 'diameter:#', 'radius:#']
        titles = [ 'Aggregation\nfactor', 'Coherence', 'Normalized\naverage\naround center', 'Average\naround\ncenter', 'Average of\naverage\nall particles', 'Diameter', 'Radius']
        Plotter.plot_heatmap(df[measures].corr(), output_filename="correlation_measures_vonneumann_%02d.png" % i, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)

    Plotter.plot_heatmap(df[measures].corr(), titles=titles, ordered=False, tight_layout=[], figsize=(12,10), cmap=plt.cm.RdBu_r, values_on=True, font_size=11)

    # the average and std dev
    all = None
    each = []
    for i in range(30):
        filename = "./data/100_particles/global_F06_%02d.teste.hdf" % i
        filename = "./data/100_particles/ring_F06_%02d.teste.hdf" % i
        filename = "./data/100_particles/dynamicring_F06_%02d.teste.hdf" % i
        filename = "./data/100_particles/vonneumann_F06_%02d.teste.hdf" % i
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
    Plotter.plot_heatmap(avg, output_filename="correlation_measures_vonneumann_.pdf", values_on_text=df, titles=titles, ordered=False, tight_layout=[], figsize=(7,6), cmap=plt.cm.RdBu_r, values_on=True, font_size=14)

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
    df = pd.read_hdf("./data/100_particles/number_of_components/heatmap_vonneumann_00500.hdf", 'df')
    df = pd.read_hdf("./data/100_particles/number_of_components/heatmap_ring_01000.hdf", 'df')
    df = pd.read_hdf("./data/100_particles/number_of_components/heatmap_dynamicring_02000.hdf", 'df')
    components = np.array(df)#/100
    ws = list(np.arange(0, 1.01, 0.01))
    yticks = range(0, len(ws), 10)
    titles_y = map(str, [ws[i] for i in yticks])
    tws = range(1, 1001, 10) + [1000]
    xticks = range(0, len(tws), 10)
    titles_x = map(str, [1]+[tws[i]-1 for i in xticks[1:len(xticks)-1]]+[1000])
    yticks.reverse()

    titles_x = None
    titles_y = None
    xticks = [0]
    yticks = [0]

    #for i in range(len(ws)):
    #for j in range(len(tws)):
    #components[i][j] = SwarmAnalyzer.get_number_of_components(filename, tw, w, calculate_on=1000)[0][1]
#    Plotter.plot_heatmap(np.flipud(components), vmin=components.min(), vmax=components.max(), figsize=(5,3.5), ordered=False, titles_x=titles_x, set_xticks=xticks, titles_y=titles_y, set_yticks=yticks, grid=False, y_label="Time window", x_label="Filter", cmap=plt.cm.Spectral_r)
#     ws = range(10)
#     tws = range(10)
#     components = np.zeros((len(ws), len(tws)))
#     for i in range(len(ws)):
#         for j in range(len(tws)):
#             components[i][j] = j + i*len(tws)


    Plotter.plot_heatmap(np.fliplr(np.rot90(components,2)), subplot_adjust=[0.08, 0.135], font_size=16, vmin=0, vmax=100, tight_layout=[], figsize=(9, 6.5), ordered=False, titles_x=titles_x, set_xticks=xticks, x_label="Time window", y_label="Filter", titles_y=titles_y, set_yticks=yticks, grid=False, cmap=plt.cm.Spectral_r)


    Plotter.plot_heatmap(np.fliplr(np.rot90(components,2)), , vmin=components.min(), vmax=components.max(), figsize=(7, 5), ordered=False, titles_x=titles_x, set_xticks=xticks, x_label="Time window", y_label="Filter", titles_y=titles_y, set_yticks=yticks, grid=False, cmap=plt.cm.Spectral_r)
    Plotter.plot_heatmap(np.fliplr(np.rot90(components,2)), vmin=components.min(), vmax=components.max(), tight_layout=[], figsize=(9, 6.5), ordered=False, titles_x=titles_x, set_xticks=xticks, x_label="Time window", y_label="Filter", titles_y=titles_y, set_yticks=yticks, grid=False, cmap=plt.cm.Spectral_r)
    # Plotter.plot_heatmap(np.rot90(components,2), vmin=0, vmax=1, figsize=(6,3.5), ordered=False, titles_x=titles_x, set_xticks=xticks, titles_y=titles_y, set_yticks=yticks, grid=False, x_label="Time window", y_label="Filter", cmap=plt.cm.Spectral_r)
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
pre_callback = Callback.to_symmetric
# for iter in [999]:
for iter in [5, 50, 100, 250, 500, 750, 999]:
    graphs, _ = SwarmParser.read_file_and_measures(filename, influence_graph_grep="ig\:#", window_size=iter, calculate_on=iter, pos_callback=pre_callback)

    graph_1_ = graphs[0][1]
    df_info = SwarmAnalyzer.get_swarm_informations_from_file(filename, measures, information_map=lambda x: x, calculate_on=iter)

    graph_2_ = SwarmParser.read_matrix_from_line(list(df_info.irow(0))[1])
    graph_1 = [graph_1_[i][j] for i in range(len(graph_1_)) for j in range(i, len(graph_1_)) if i != j] # np.reshape(graph_1_, (1, 10000))
    graph_2 = [graph_2_[i][j] for i in range(len(graph_2_)) for j in range(i, len(graph_2_)) if i != j] # np.reshape(graph_2_, (1, 10000))

    graph_1 = np.array(graph_1)
    graph_2 = np.array(graph_2)

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
    df = pd.DataFrame({'x': graph_2/float(max(graph_2)), 'y': graph_3/float(max(graph_3))})
    xlabel = "Geodesic distance"
    ylabel = "Euclidean distance"
    title = "Global"
    output_filename="global_%04d.pdf" % iter
    output_filename = None
    Plotter.plot_curve(df, linewidth=0, figsize=(5, 3.8), tight_layout=[], output_filename=output_filename, annotate=(0.15, 0.84, "$R=%0.2f$" % (df.corr()['y']['x'])), font_size=11, alpha=0.6, marker='.', colors='#de2d26', grid=True, markersize=6, title=title, x_label=xlabel, y_label=ylabel, xlim=(0,1), ylim=(0,1))
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
        pd_datas = {}
        for title in graph_matrices:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(graphs,
                                                                             weight_normalize=normalize)
            pd_datas[title] = dict(zip(windows_size, curves_areas))
        GiantComponentDeathPlotter.giant_component_death_curve(calculate_on, pd_datas, windows_size, xlim=(0, 1.0))
    """
    execfile("swarm_analyzer.py")
    filename = [('global', "/mnt/pso_100_particles/global_F06_00"), ('g1', "/mnt/pso_100_particles/global_F06_01"), \
    ('g2', "/mnt/pso_100_particles/global_F06_02"), ('g3', "/mnt/pso_100_particles/global_F06_03")]
    df = SwarmAnalyzer.read_files_and_plot(filename, windows_size=[50, 100, 1000], calculate_on=999)
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


