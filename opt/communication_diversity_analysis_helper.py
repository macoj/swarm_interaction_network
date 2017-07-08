from swarm_analyzer import SwarmAnalyzer

__author__ = 'marcos'


class CommunicationDiversityHelper():
    def __init__(self):
        pass

    @staticmethod
    def calculate_communication_diversity(from_=None, to_=None):
        file_names = []
        runs = 30
        pks = [2, 25, 50, 75, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        for topology in [("dynamic", 1)]:  # initial is ring
            for function in [21, 22]:
                for pk in pks:
                    file_names += ["%s_%04d_F%02d_%02d" % (topology[0], pk, function, r) for r in range(runs)]
        for function in [21, 22]:
            for topology in ["ring", "global"]:
                file_names += ["%s_F%02d_%02d" % (topology, function, r) for r in range(runs)]

        if from_ is not None and to_ is not None:
            file_names = file_names[from_:to_]
        for f in file_names:
            filename = './data/analysis_on_a_dynamic_topology/%s' % f
            cd = SwarmAnalyzer.communication_diversity(filename, until=2000)
            cd.to_hdf(filename+"_cd.hdf", 'df')
            del cd
    """
execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=0, to_=106)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=106, to_=212)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=212, to_=318)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=318, to_=424)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=424, to_=530)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=530, to_=636)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=636, to_=742)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=742, to_=848)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=848, to_=960)
    """