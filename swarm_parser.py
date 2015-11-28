__author__ = 'marcos'
import re
import sys
import math
import numpy as np


class SwarmParser:
    def __init__(self):
        pass

    @staticmethod
    def read_matrix_from_line(line, sep=' '):
        """ Reads a line a creates a square matrix.
        :param line: The string line.
        :param sep: The separator between the numbers. Default = ' '.
        :return: A numpy array.
        """
        matrix_result = None
        if line is not None and type(line) == str:
            matrix_read = line.strip().split(sep)
            matrix_read = list(map(float, matrix_read))
            number_of_elements = len(matrix_read)
            dimension_float = math.sqrt(number_of_elements)
            dimension = math.trunc(dimension_float)
            if dimension == dimension_float:
                matrix_slices = [matrix_read[part:part+dimension] for part in range(0, number_of_elements, dimension)]
                matrix_result = np.array(matrix_slices)
            else:
                # this matrix cannot be assembled as a square matrix
                pass
            del matrix_read
        return matrix_result

    @staticmethod
    def read_vector_from_line(line, sep=' '):
        """ Reads a string and converts to a vector.
        :param line: The string line.
        :param sep: The separator between the numbers. Default = ' '.
        :return: A numpy array.
        """
        matrix_result = None
        if line is not None and type(line) == str:
            matrix_read = line.strip().split(sep)
            matrix_read = list(map(float, matrix_read))
            matrix_result = np.array(matrix_read)
        return matrix_result

    @staticmethod
    def read_file_and_measures(filename,
                               influence_graph_grep=None,
                               informations_grep=None,
                               window_size=-1,
                               pre_callback=None,
                               pos_callback=None):
        """
        """
        input_file = open(filename, 'r')
        windowed = window_size >= 1
        window = {}
        matrix_count = 0
        graphs = []
        informations = None
        if informations_grep:
            if type(informations_grep) != list:
                informations_grep = [informations_grep]
            informations = dict([(i, []) for i in informations_grep])
        accumulated_matrix = None
        for line in input_file:
            # we add fitnesses and influence graphs in ig_line and fitnesses
            matrix_line, information, iteration = SwarmParser.grep_line_infos(line,
                                                                              influence_graph_grep,
                                                                              informations_grep)
            information_index, information = information
            if information:
                informations[information_index].append((iteration, information))
            if matrix_line:
                current_matrix = SwarmParser.read_matrix_from_line(matrix_line)
                if pre_callback:
                    current_matrix = pre_callback(current_matrix)
                if not windowed:
                    if accumulated_matrix is None:
                        accumulated_matrix = current_matrix
                    else:
                        accumulated_matrix = accumulated_matrix + current_matrix
                    current_accumulated = accumulated_matrix
                else:
                    window[matrix_count % window_size] = current_matrix
                    current_accumulated = SwarmParser.sum_matrices(window)
                if pos_callback:
                    current_accumulated = pos_callback(current_accumulated)
                graphs.append((iteration, current_accumulated))
                matrix_count += 1
        input_file.close()
        return dict(graphs), informations
    """
execfile("swarm_parser.py")
t = SwarmParser.read_file_and_measures("/mnt/pso_100_particles/global_F06_00", influence_graph_grep="ig\:#", window_size=500, informations_grep=["radius\:#", "it\:#"])
    """

    @staticmethod
    def read_files_and_measures(title_filenames, influence_graph_grep=None, informations_grep=None,
                                window_size=-1, pre_callback=None, pos_callback=None):
        all_graph_matrices = {}
        all_informations = {}
        for filename in title_filenames:
            print filename
            title, filename = filename
            graphs, informations = SwarmParser.read_file_and_measures(filename,
                                                                      influence_graph_grep=influence_graph_grep,
                                                                      informations_grep=informations_grep,
                                                                      window_size=window_size,
                                                                      pre_callback=pre_callback,
                                                                      pos_callback=pos_callback)
            all_graph_matrices[title] = graphs
            all_informations[title] = informations
        return all_graph_matrices, all_informations

    @staticmethod
    def grep_line_infos(line, influence_graph_grep, informations_grep):
        iteration_grep = "^[0-9]* "
        graph_line, information, information_grep = None, None, None
        iteration = -1
        line, times = re.subn(influence_graph_grep, "", line)
        if times != 0:
            iteration_find = re.findall(iteration_grep, line)
            if iteration_find:
                iteration = int(iteration_find[0].strip())
            line, _ = re.subn(iteration_grep, "", line)
            graph_line = line
        elif informations_grep:
            for information_grep in informations_grep:
                line, times = re.subn(information_grep, "", line)
                if times != 0:
                    break
            if times != 0:
                iteration_find = re.findall(iteration_grep, line)
                if iteration_find:
                    iteration = int(iteration_find[0].strip())
                line, _ = re.subn(iteration_grep, "", line)
                information = float(line.strip())
            else:
                information_grep, information = None, None
        return graph_line, (information_grep, information), iteration

    @staticmethod
    def sum_matrices(window):
        """

        :rtype : object
        """
        sum_matrix = None
        if window is not None:
            sum_matrix = np.zeros(window[0].shape)
        for w in window:
            sum_matrix = sum_matrix + window[w]
        return sum_matrix

    # @staticmethod
    # def grep_line(line, ig_greped, fitness_greped, ig_line, fitnesses, influence_graph_grep=None, fitness_grep=None):
    #     if influence_graph_grep is not None and not ig_greped:
    #         line, times = re.subn(influence_graph_grep, "", line)
    #         ig_greped = (times != 0)
    #         if ig_greped:
    #             ig_line.append(line)
    #     if fitness_grep is not None and not fitness_greped:
    #         fitness, times = re.subn(fitness_grep, "", line)
    #         fitness_greped = fitness_greped or (times != 0)
    #         if fitness_greped:
    #             #print str(fitness)
    #             fitnesses.append(float(fitness.strip()))
    #     return ig_greped, fitness_greped, fitnesses, ig_line
#     @staticmethod
#     def measure(matrix, pre_callback, calculate, pos_callback, matrix_count, fitnesses=None):
#         # matrix here can be symmetric or not, pre_callback should work on it...
#         sum_matrix_measured = matrix
#         if pre_callback is not None:
#             sum_matrix_measured = pre_callback(matrix)
#         if calculate is not None:
#             sum_matrix_measured = calculate(sum_matrix_measured)
#         # adds info about line read
#         matrix_out = matrix_count, sum_matrix_measured
#         ##
# #             if (matrix_count == 100):
# #                 return (sum_matrix_measured)
# #             continue
#         ##
#         if pos_callback == sys.stdout.write:
#             matrix_out = str(matrix_out).strip() + "\n"
#
#         if pos_callback is not None:
#             if fitnesses:
#                 pos_callback(matrix_out, fitnesses)
#             else:
#                 pos_callback(matrix_out)
#
#     @staticmethod
#     def measures(matrix, pre_callback, calculate, pos_callback):
#         # matrix here can be symmetric or not, pre_callback should work on it...
#         matrix_measured = matrix
#         if pre_callback is not None:
#             matrix_measured = pre_callback(matrix)
#         if calculate is not None:
#             matrix_measured = calculate(matrix_measured)
#         if pos_callback is not None:
#             if pos_callback == sys.stdout.write:
#                 matrix_measured = str(matrix_measured).strip() + "\n"
#             pos_callback(matrix_measured)
#         return matrix_measured

# @staticmethod
    # def get_influence_graph_and_particles_position(filename,
    #                                                position_grep=None,
    #                                                influence_graph_grep=None,
    #                                                window_size=-1,
    #                                                calculate_on=-1):
    #     """ Gets the influence graph and the particles positions in an iteration given.
    #     :param filename:
    #     :param position_grep:
    #     :param influence_graph_grep:
    #     :param window_size:
    #     :param calculate_on:
    #     :return:
    #     """
    #     input_file = open(filename, 'r')
    #     #window_size = 10
    #     windowed = window_size >= 1
    #     window = {}
    #     window_current = 0
    #
    #     matrix_count = 0
    #     accumulated_matrix = None
    #     particles_position = {}
    #     ig_pp_return = None
    #
    #     # let's get just the dimension of the influence graph
    #     for line in input_file:
    #         line, times = re.subn(influence_graph_grep, "", line)
    #         ig_greped = (times != 0)
    #         if ig_greped:
    #             accumulated_matrix = SwarmParser.read_matrix_from_line(line)
    #             accumulated_matrix = np.zeros(accumulated_matrix.shape)
    #             break
    #
    #     # rewind
    #     input_file.seek(0L)
    #
    #     graph_done = False
    #     if accumulated_matrix is not None:
    #         for line in input_file:
    #             # tries to get influence graph
    #             line, times = re.subn(influence_graph_grep, "", line)
    #             ig_greped = (times != 0)
    #             if ig_greped and graph_done is False:
    #                 matrix_count += 1
    #                 matrix = SwarmParser.read_matrix_from_line(line)
    #                 accumulated_matrix = accumulated_matrix + matrix
    #                 # let's keep the window history
    #                 if windowed:
    #                     window_current += 1
    #                     window[window_current % window_size] = matrix
    #             else:
    #                 # let's try to get a position vector,
    #                 # the pattern is something like: position:#0#i x y z
    #                 # what means that the position_grep = position:#
    #                 line, times = re.subn(position_grep, "", line)
    #                 pos_greped = (times != 0)
    #                 if pos_greped:
    #                     pos_match = re.match('^[0-9]*#[0-9]*', line)
    #                     if pos_match is not None:
    #                         pos_match = pos_match.group(0)
    #                         iteration, particle_id = pos_match.split('#')
    #                         if int(iteration) == calculate_on:
    #                             line_mod, _ = re.subn('^[0-9]*#[0-9]*', "", line)
    #                             particles_position[particle_id] = SwarmParser.read_vector_from_line(line_mod)
    #             if matrix_count == calculate_on:
    #                 graph_done = True   # and we do not care about the graph anymore
    #                 if len(particles_position) != accumulated_matrix.shape[0]:
    #                     continue
    #                 if windowed:
    #                     graph_return = SwarmParser.sum_matrices(window)
    #                 else:
    #                     graph_return = accumulated_matrix
    #                 ig_pp_return = (graph_return, particles_position)
    #                 break
    #     input_file.close()
    #     return ig_pp_return

    # @staticmethod
    # def read_file_line_and_measure(filename, calculate, grep=None, pre_callback=None, pos_callback=sys.stdout.write):
    #     """ Reads each line in a file, creates a graph for each one based on the line content and measure the graph.
    #     :param filename: The filename of the file to be processed.
    #     :param calculate: The function used to measure the graph.
    #     :param grep: A regexp to use just specific lines of the file. The match is removed from the line.
    #     :param pre_callback: A function used to pre-process the graph.
    #     :param pos_callback: A function used to receive the mesasurement. If none, sys.stdout.write is used.
    #     :return: None.
    #     """
    #     input_file = open(filename, 'r')
    #     for line in input_file:
    #         if grep is not None:
    #             line, times = re.subn(grep, "", line)
    #             if times == 0:
    #                 continue
    #         matrix = SwarmParser.read_matrix_from_line(line)
    #         if pre_callback is not None:
    #             matrix = pre_callback(matrix)
    #             #print matrix
    #         matrix_measured = calculate(matrix)
    #         if pos_callback == sys.stdout.write:
    #             matrix_measured = str(matrix_measured) + "\n"
    #         pos_callback(matrix_measured)
    #     input_file.close()
#     '''
#  SwarmParser.read_file_and_measure(
#  '/home/marcos/PhD/research/pso_influence_graph_communities/pso_dynamic_initial_ring_F6_30',
#  SpectraUndirect.calculate, grep = 'ig\:#[0-9]*')
#  SwarmParser.read_file_and_measure(
#  '/home/marcos/PhD/research/pso_influence_graph_communities/pso_dynamic_initial_ring_F6_30',
# igraph.Graph.degree, grep = 'ig\:#[0-9]*', pre_callback = create_igraph_from_matrix)
#
# SwarmParser.read_file_and_measure_no_window(
# '/home/marcos/PhD/research/pso_influence_graph_communities/pso_dynamic_initial_ring_F6_30',
# igraph.Graph.degree, grep = 'ig\:#[0-9]*', pre_callback = create_igraph_from_matrix)
# SwarmParser.read_file_and_measure_no_window(
# '/home/marcos/PhD/research/pso_influence_graph_communities/50_particles/pso_dynamic_initial_ring_F6_30',
# calculate = None, grep = 'ig\:#[0-9]*', pre_callback = None, pos_callback = create_and_save_plot)
#     '''
    # @staticmethod
    # def read_file_and_measure(filename,
    #                           calculate=None,
    #                           influence_graph_grep=None,
    #                           fitness_grep=None,
    #                           window_size=-1,
    #                           pre_callback=None,
    #                           pos_callback=sys.stdout.write,
    #                           calculate_on=-1):
    #     """ Measures each pair (fitnesses, graph) from the file.
    #
    #     This function uses a regexp for fitness (or any other kind of information related to the graph) and a regexp
    #     for the graph. So, after finding a pair (fitness, graph), the graph is measured and this measurement is passed
    #     to the pos_callback function with the fitness values found so far. Thus, the pos_callback function must handle
    #     two arguments, where the first one is the result of the calculate function and the second one is a list of
    #     fitness values.
    #
    #     :param filename:
    #     :param calculate:
    #     :param influence_graph_grep:
    #     :param fitness_grep:
    #     :param window_size:
    #     :param pre_callback:
    #     :param pos_callback:
    #     :return:
    #     """
    #     input_file = open(filename, 'r')
    #     #window_size = 10
    #     windowed = window_size >= 1
    #     window = {}
    #     window_current = 0
    #
    #     # gets first line in order to create the sum_matrix
    #     matrix_count = 0
    #     fitnesses = None
    #     if fitness_grep:
    #         fitnesses = []
    #     ig_greped = False
    #     fitness_greped = False
    #     accumulated_matrix = None
    #     ig_line = []
    #     for line in input_file:
    #         # we add fitnesses and influence graphs in ig_line and fitnesses
    #         ig_greped, fitness_greped, fitnesses, ig_line = SwarmParser.grep_line(line,
    #                                                                               ig_greped, fitness_greped,
    #                                                                               ig_line, fitnesses,
    #                                                                               influence_graph_grep,
    #                                                                               fitness_grep)
    #         # we will keep reading until we find at least a pair ig/fitness
    #         if not ig_greped or (not fitness_greped and fitness_grep is not None):
    #             continue
    #         accumulated_matrix = SwarmParser.read_matrix_from_line(ig_line[len(ig_line)-1])
    #         matrix_count += 1
    #         sum_matrix_measured = accumulated_matrix
    #         # let's keep the window history
    #         if windowed:
    #             window[window_current % window_size] = sum_matrix_measured
    #             window_current += 1
    #             #print str(window)
    #             sum_matrix_measured = SwarmParser.sum_matrices(window)
    #         # does it calculate all the time? or only one shot?
    #         is_to_calculate = calculate_on == -1 or (calculate_on != -1 and matrix_count == calculate_on)
    #         # let's calculate
    #         if (matrix_count >= window_size or not windowed) and is_to_calculate:
    #             SwarmParser.measure(sum_matrix_measured, pre_callback, calculate,
    #                                 pos_callback, matrix_count, fitnesses)
    #         break
    #
    #     # now we go to the other lines
    #     ig_greped = False
    #     fitness_greped = False
    #     ig_line = []
    #     # did we get the first line? did we already do what we needed?
    #     if accumulated_matrix is not None and not (calculate_on != -1 and matrix_count == calculate_on):
    #         for line in input_file:
    #             ig_greped, fitness_greped, fitnesses, line = SwarmParser.grep_line(line,
    #                                                                                ig_greped, fitness_greped,
    #                                                                                ig_line, fitnesses,
    #                                                                                influence_graph_grep,
    #                                                                                fitness_grep)
    #             # we will go until find ig and fitness pair
    #             if not ig_greped or (not fitness_greped and fitness_grep is not None):
    #                 continue
    #             matrix_count += 1
    #             matrix = SwarmParser.read_matrix_from_line(ig_line[len(ig_line)-1])
    #             accumulated_matrix = accumulated_matrix + matrix
    #             sum_matrix_measured = accumulated_matrix
    #             # reset:
    #             ig_greped = False
    #             fitness_greped = False
    #             ig_line = []
    #             # let's keep the window history
    #             if windowed:
    #                 window_current += 1
    #                 window[window_current % window_size] = matrix
    #                 sum_matrix_measured = SwarmParser.sum_matrices(window)
    #             # does it calculate all the time? or only one shot?
    #             is_to_calculate = calculate_on == -1 or (calculate_on != -1 and matrix_count == calculate_on)
    #             # let's calculate
    #             if (matrix_count >= window_size or not windowed) and is_to_calculate:
    #                 SwarmParser.measure(sum_matrix_measured, pre_callback, calculate,
    #                                     pos_callback, matrix_count, fitnesses)
    #             # already done?
    #             if calculate_on != -1 and matrix_count >= calculate_on:
    #                 break
    #     input_file.close()