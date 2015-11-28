__author__ = 'marcos'
import numpy as np


class Callback:
    def __init__(self):
        pass

    @staticmethod
    def to_symmetric(matrix):
        is_symmetric = np.triu(matrix) == np.transpose(np.tril(matrix))
        is_symmetric = is_symmetric.min()
        # here we sum
        if not is_symmetric:
            matrix = matrix + np.transpose(np.triu(matrix)) + np.transpose(np.tril(matrix))
        return matrix

