__author__ = 'marcos'
from swarm_analyzer import SwarmAnalyzer
from giant_component_analysis import GiantComponentDeath
from opt.giant_component_analysis_plotter import GiantComponentDeathPlotter
from swarm_parser import SwarmParser
from opt.plotter import Plotter


class SwarmPlotter():
    def __init__(self):
        pass

    @staticmethod
    def read_files_and_plot_destruction_curves(filenames, windows_size, calculate_on, count='components'):
        graph_matrices = SwarmAnalyzer.get_graph_matrices_from_files(
            filenames,  windows_size=windows_size, calculate_on=calculate_on)
        normalize = [2 * i for i in windows_size]
        pd_datas = []
        for title, _ in filenames:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(
                graphs, weight_normalize=normalize, count=count)
            pd_datas.append((title, dict(zip(windows_size, curves_areas))))
        GiantComponentDeathPlotter.giant_component_death_curve(
            calculate_on, pd_datas, xlim=(0, 1.0), figsize=(4.5, 4), count=count)
    """
    from opt.plotter import Plotter
    Plotter.plos_style()
    execfile("swarm_analyzer.py")
    execfile("swarm_plotter.py")
    filenames = [('Global', "./data/global_F06_15"), ('Ring', "./data/ring_F06_15"), ('Von Neumann', "./data/vonneumann_F06_15"), ('Dynamic', "./data/dynamicring_F06_15")]
    df = SwarmPlotter.read_files_and_plot_destruction_curves(filenames, windows_size=[10, 50, 100], calculate_on=100, count='size')
    """

    @staticmethod
    def plot_boxplot_fitness(
            sets_of_filenames, output_filename=None, info_grep="it\:#", yticks_args=None, legends=None, showmeans=False,
            showfliers=True, whis=1., widths=0.7, bootstrap=2000, **kargs):
        """
        :param sets_of_filenames:
        :param output_filename:
        :param info_grep:
        :return:
        """
        get_infos = lambda x: [
            SwarmParser.read_file_and_measures(filename, informations_grep=info_grep)[1][info_grep][-1][1]
            for filename in x]
        values = map(get_infos, sets_of_filenames)
        print values
        boxes_kargs = {'color': 'black', 'linewidth': 1.3, 'zorder': 3, 'fillstyle': 'full', 'facecolor': '#a6bddb'}
        means_kargs = {'color': 'black', 'fillstyle': 'full', 'markerfacecolor': "black", 'marker': "s",
                       'markersize': 3, 'mew': 1, 'mec': 'black', 'zorder': 5}
        fliers_kargs = {'color': 'black', 'marker': "s", 'markersize': 1, 'mew': 1.2, 'mec': 'black'}
        whiskers_kargs = {'color': 'black', 'linewidth': 1.2, 'zorder': 2, 'linestyle': '-', 'alpha': 1.0,
                          'mec': 'black'}
        medians_kargs = {'color': 'black', 'linewidth': 1.6, 'zorder': 5, 'alpha': 0.3}
        caps_kargs = {'linewidth': 1.5, 'color': 'black'}
        xticks_args = None
        if legends is not None:
            xticks_args = (range(1, len(legends)+1), legends)
        Plotter.plot_boxplots(
            values, grid_only='y', xlim=(0.3, len(sets_of_filenames) + 0.6), yticks_args=yticks_args, whis=whis,
            xticks_args=xticks_args, grid=False, widths=widths, bootstrap=bootstrap, boxes_kargs=boxes_kargs,
            showmeans=showmeans, showfliers=showfliers, fliers_kargs=fliers_kargs, means_kargs=means_kargs,
            whiskers_kargs=whiskers_kargs, medians_kargs=medians_kargs, caps_kargs=caps_kargs, output=output_filename,
            **kargs)