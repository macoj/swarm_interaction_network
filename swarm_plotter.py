__author__ = 'marcos'
from swarm_analyzer import SwarmAnalyzer
from giant_component_analysis import GiantComponentDeath
from opt.giant_component_analysis_plotter import GiantComponentDeathPlotter


class SwarmPlotter():
    def __init__(self):
        pass

    @staticmethod
    def read_files_and_plot_destruction_curves(filenames, windows_size, calculate_on):
        graph_matrices = SwarmAnalyzer.get_graph_matrices_from_files(
            filenames,  windows_size=windows_size, calculate_on=calculate_on)
        normalize = [2 * i for i in windows_size]
        pd_datas = []
        for title, _ in filenames:
            graphs = [graph_matrices[title][i] for i in windows_size]
            graphs = map(lambda x: x[0], graphs)  # this was a calculate_on call
            curves_areas = GiantComponentDeath.create_giant_component_curves(graphs, weight_normalize=normalize)
            pd_datas.append((title, dict(zip(windows_size, curves_areas))))
        GiantComponentDeathPlotter.giant_component_death_curve(
            calculate_on, pd_datas, windows_size, xlim=(0, 1.0), figsize=(4.5, 4))
    """
    execfile("swarm_analyzer.py")
    execfile("swarm_plotter.py")
    filenames = [('Global', "./data/global_F06_15"), ('Ring', "./data/ring_F06_15"), ('Von Neumann', "./data/vonneumann_F06_15"), ('Dynamic', "./data/dynamicring_F06_15")]
    df = SwarmPlotter.read_files_and_plot_destruction_curves(filenames, windows_size=[100, 1000], calculate_on=1000)
    """
