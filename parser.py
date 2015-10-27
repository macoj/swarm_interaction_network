__author__ = 'marcos'
import math
import numpy as np


class Parser:
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