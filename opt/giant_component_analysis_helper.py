__author__ = 'marcos'
import pandas as pd
from swarm_analyzer import SwarmAnalyzer


class GiantComponentDeathHelper:
    def __init__(self):
        pass

    @staticmethod
    def get_number_of_component_different_time_windows(filename, calculate_on=1000, tws=None):
        curves = []
        if tws is None:
            tws = range(1, calculate_on, 10) + [calculate_on]
        for tw in tws:
            col_x = "x_%04d" % tw
            col_y = "y_%04d" % tw
            curve = SwarmAnalyzer.get_giant_component_destruction_curves(
                filename, tw, calculate_on=calculate_on, count='components', adjusted=False)
            curve = curve[0].append([{'x': 1, 'y': curve[0]['y'].max()}])
            curve.index = range(len(curve))
            curve.columns = [col_x, col_y]
            curves.append(curve)
        df = pd.concat(curves, axis=1)
        return df
    """
execfile("opt/giant_component_analysis_helper.py")
for t in ['ring', 'vonneumann', 'global', 'dynamicring']:
    filename = './data/%s_F06_15' % t
    df = GiantComponentDeathHelper.get_number_of_component_different_time_windows(filename)
    df.to_hdf(filename + "_components_.hdf", 'df')
    """
