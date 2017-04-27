__author__ = 'marcos'
import re
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
    def read_file_and_measures(
            filename, influence_graph_grep=None, informations_grep=None, information_map=float, window_size=-1,
            pre_callback=None, pos_callback=None, calculate_on=-1, until=-1):
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
            matrix_line, information, iteration = SwarmParser.grep_line_infos(
                line, influence_graph_grep, informations_grep, information_map=information_map)
            # let's check if we already have useful information
            if iteration == -1:
                continue
            # OK, let's start!
            information_index, information = information
            if information is not None:
                if calculate_on == -1 or iteration == calculate_on:
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
                    window[matrix_count % window_size] = current_matrix  # we can do that because order does not matter
                    current_accumulated = SwarmParser.sum_matrices(window)
                    matrix_count += 1
                if calculate_on == -1 or iteration == calculate_on:
                    if pos_callback:
                        current_accumulated = pos_callback(current_accumulated)
                    graphs.append((iteration, current_accumulated))
            if iteration == calculate_on:
                break
            if until != -1 and iteration >= until:
                break
        input_file.close()
        return graphs, informations
    """
execfile("swarm_parser.py")
t = SwarmParser.read_file_and_measures("/mnt/pso_100_particles/global_F06_00", influence_graph_grep="ig\:#", window_size=500, informations_grep=["radius\:#", "it\:#"])
    """

    @staticmethod
    def read_files_and_measures(
            title_filenames, influence_graph_grep=None, informations_grep=None, information_map=float, windows_size=-1,
            pre_callback=None, pos_callback=None, calculate_on=-1, until=-1):
        all_graph_matrices = {}
        all_informations = {}
        if type(windows_size) != list:
            windows_size = [windows_size]
        for filename in title_filenames:
            title, filename = filename
            all_graph_matrices[title] = {}
            all_informations[title] = {}
            for window_size in windows_size:
                graphs, informations = SwarmParser.read_file_and_measures(
                    filename, influence_graph_grep=influence_graph_grep, informations_grep=informations_grep,
                    information_map=information_map, window_size=window_size, pre_callback=pre_callback,
                    pos_callback=pos_callback, calculate_on=calculate_on, until=until)
                if graphs:
                    all_graph_matrices[title][window_size] = graphs
                if informations:
                    all_informations[title][window_size] = informations
        return all_graph_matrices, all_informations

    @staticmethod
    def grep_line_infos(line, influence_graph_grep, informations_grep, information_map=float):
        iteration_grep = "^[0-9]* "
        graph_line, information, information_grep = None, None, None
        iteration = -1
        times = 0
        if influence_graph_grep:
            line, times = re.subn("^"+influence_graph_grep, "", line)
        if times != 0:
            iteration_find = re.findall(iteration_grep, line)
            if iteration_find:
                iteration = int(iteration_find[0].strip())
            line, _ = re.subn(iteration_grep, "", line)
            graph_line = line
        elif informations_grep:
            for information_grep in informations_grep:
                line, times = re.subn("^"+information_grep, "", line)
                if times != 0:
                    break
            if times != 0:
                iteration_find = re.findall(iteration_grep, line)
                if iteration_find:
                    iteration = int(iteration_find[0].strip())
                line, _ = re.subn(iteration_grep, "", line)
                information = information_map(line.strip())
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
