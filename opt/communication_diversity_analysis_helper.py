from swarm_analyzer import SwarmAnalyzer

__author__ = 'marcos'


class CommunicationDiversityHelper():
    def __init__(self):
        pass

    @staticmethod
    def calculate_communication_diversity(from_=None, to_=None):
        file_names = []
        runs = 1
        pks = [2, 25, 50, 75, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        for topology in [("dynamic", 1)]:  # initial is ring
            for function in [26, 25, 24, 23]:  # sphere, schwefel, rosenbrock, rastrigin
                for pk in pks:
                    file_names += ["%s_%04d_F%02d_%02d" % (topology[0], pk, function, r) for r in range(runs)]
        if from_ is not None and to_ is not None:
            file_names = file_names[from_:to_]
        for f in file_names:
            filename = './analysis_on_a_dynamic_topology/%s' % f
            cd = SwarmAnalyzer.communication_diversity(filename, until=2000)
            cd.to_hdf(filename+"_cd.hdf", 'df')
            del cd