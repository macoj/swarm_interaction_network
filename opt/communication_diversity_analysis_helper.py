from swarm_analyzer import SwarmAnalyzer

__author__ = 'marcos'


class CommunicationDiversityHelper():
    def __init__(self):
        pass

    @staticmethod
    def calculate_communication_diversity(from_=None, to_=None):
        file_names = []
        # runs = 1
        # pks = [2, 25, 50, 75, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        # for topology in [("dynamic", 1)]:  # initial is ring
        #     # for function in [26, 25, 24, 23]:  # sphere, schwefel, rosenbrock, rastrigin
        #     for function in [26, 24]:  # sphere, schwefel, rosenbrock, rastrigin
        #         for pk in pks:
        #             file_names += ["%s_%04d_F%02d_%02d" % (topology[0], pk, function, r) for r in range(10, 30)]
        #     for function in [25, 23]:  # sphere, schwefel, rosenbrock, rastrigin
        #         for pk in pks:
        #             file_names += ["%s_%04d_F%02d_%02d" % (topology[0], pk, function, r) for r in range(1, 30)]
        for function in [26, 25, 24, 23]:
            for topology in ["ring", "global"]:
                file_names += ["%s_F%02d_%02d" % (topology, function, r) for r in range(1)]
        for function in [25, 24, 23]:
            for topology in ["global"]:
                file_names += ["%s_F%02d_%02d" % (topology, function, r) for r in range(1, 30)]

        if from_ is not None and to_ is not None:
            file_names = file_names[from_:to_]
        for f in file_names:
            filename = './data/analysis_on_a_dynamic_topology/%s' % f
            cd = SwarmAnalyzer.communication_diversity(filename, until=2000)
            cd.to_hdf(filename+"_cd.hdf", 'df')
            del cd
    """
execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=0, to_=11)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=11, to_=22)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=22, to_=33)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=33, to_=44)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=44, to_=55)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=55, to_=66)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=66, to_=77)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=77, to_=88)

execfile("opt/communication_diversity_analysis_helper.py")
CommunicationDiversityHelper.calculate_communication_diversity(from_=88, to_=95)
    """